# The ATLS thermal-damage layer (traitecoevo/plant#566), behind
# `use_thermal_damage_`.
#
# The leaf half of a lasting-damage model, functional <=> impaired -> damaged: a
# caller integrates the acclimation state A_ and the damaged fractions I_r_ and
# I_p_, seats them, and reads back phi_i_ (the impaired fraction), k_i_eff_ (the
# leak from impaired into damage) and k_rec_eff_ (the repair rate) at the leaf's
# temperature. The equations are a port of plant's July layer; the checks
# ported from its test-leaf-thermal.R are the first half of this file.
#
# ⚠️ TWO KINDS OF ASSERTION, and both are needed. Off, the layer must be INERT,
# which is asserted with identical() against inputs set to values that would
# matter if they were read. On, it must DO what it says, which is asserted
# against the equations at non-default values -- otherwise an inert gate would
# pass the first half by doing nothing at all.

# The gas constant and Kelvin offset the leaf uses. A reference built on R = 8.314
# would disagree with phi_i_ by more than rounding.
.R_gas <- 8.314462618153240
.T0 <- 273.15

thermal_leaf <- function(supply = leaf_supply_singlelayer()) {
  leaf_model(traits = leaf_traits(), control = leaf_control(), supply = supply)
}

drive_at <- function(l, leaf_temp, psi_soil = 0.5) {
  if (identical(l$supply_kind, "single")) {
    set_drivers(l, psi_soil = psi_soil, PPFD = 1200,
                root_network = series_resistance(50),
                leaf_specific_conductance_max = 1e-4, atm_vpd = 1.0, ca = 40,
                leaf_temp = leaf_temp, atm_o2_kpa = 21, atm_kpa = 101.3)
  } else {
    set_drivers(l, psi_soil = psi_soil + c(0, 0.4, 0.8), PPFD = 1200,
                leaf_specific_conductance_max = 1e-4, atm_vpd = 1.0, ca = 40,
                leaf_temp = leaf_temp, atm_o2_kpa = 21, atm_kpa = 101.3)
  }
  invisible(l)
}

# K/(1+K), K = exp[(H/R)(1/Tmid - 1/T)], Tmid in deg C. The defaults are the
# leaf's: H = 200 kJ/mol and Tmid = 200000/650 K, which is where the port's
# Jmax-tied equilibrium put its midpoint.
phi_i_ref <- function(leaf_temp, H = 200000, Tmid = 200000 / 650 - .T0) {
  K <- exp((H / .R_gas) * (1 / (Tmid + .T0) - 1 / (leaf_temp + .T0)))
  K / (1 + K)
}

# The Jmax curve's own deactivation equilibrium, K/(1+K) with K = exp((d_S*T -
# H_d)/(R*T)): what the impaired fraction USED to be, and what its defaults
# still reproduce.
jmax_deactivation_ref <- function(leaf_temp, d_S = 650, H_d = 200000) {
  T_K <- leaf_temp + .T0
  K <- exp((d_S * T_K - H_d) / (.R_gas * T_K))
  K / (1 + K)
}

# Every input of the layer, at a value other than its default.
atls_inputs <- list(topt_offset_ = 2.5, impair_H_ = 2.6e5, impair_Tmid_ = 37,
                    m_rep_ = 0.7, t_rep_cut_ = 41,
                    k_i_ = 0.2, k_i_Ea_ = 8e4, k_i_Tref_ = 30,
                    k_i_baseline_ = 0.6, k_rec_ = 0.4, dTopt_max_ = 4,
                    dTimpair_max_ = 2.5, K_A_ = 2, A_ = 1.5, I_r_ = 0.2,
                    I_p_ = 0.15,
                    c_acclim_maint_ = 0.05, c_repair_maint_ = 0.3,
                    c_protect_maint_ = 0.4, c_repair_flux_ = 0.8)

set_inputs <- function(l, inputs = atls_inputs) {
  for (nm in names(inputs)) l[[nm]] <- inputs[[nm]]
  invisible(l)
}

# --- ported from plant's July test-leaf-thermal.R ---------------------------

test_that("thermal damage is off by default and writes zero outputs", {
  l <- thermal_leaf()
  expect_false(l$use_thermal_damage_)
  drive_at(l, 30)
  expect_identical(l$phi_i_, 0)
  expect_identical(l$k_i_eff_, 0)
  expect_identical(l$k_rec_eff_, 0)

  # On but undamaged and unshifted: the shift is an exact zero in exact
  # arithmetic, and d_S_shifted's round trip H_d/(H_d/(d_S - corr)) + corr
  # returns d_S to rounding only, hence a tolerance rather than identity.
  off <- drive_at(thermal_leaf(), 45)
  on <- thermal_leaf()
  on$use_thermal_damage_ <- TRUE
  drive_at(on, 45)
  expect_equal(on$vcmax_, off$vcmax_, tolerance = 1e-9)
  expect_equal(on$jmax_, off$jmax_, tolerance = 1e-9)
})

test_that("phi_i is a centred Arrhenius equilibrium with its own midpoint", {
  for (temp in c(20, 34.54, 45, 55)) {
    l <- thermal_leaf()
    l$use_thermal_damage_ <- TRUE
    drive_at(l, temp)
    expect_equal(l$phi_i_, phi_i_ref(temp), tolerance = 1e-12,
                 info = paste("leaf_temp", temp))
    # The defaults reproduce the port, where the onset was the Jmax curve's
    # deactivation equilibrium, to rounding.
    expect_equal(l$phi_i_, jmax_deactivation_ref(temp), tolerance = 1e-12,
                 info = paste("leaf_temp", temp))
  }
  # Half-impaired (K = 1) exactly at the midpoint, which is a parameter now.
  l <- thermal_leaf()
  l$use_thermal_damage_ <- TRUE
  drive_at(l, l$impair_Tmid_)
  expect_equal(l$phi_i_, 0.5, tolerance = 1e-12)
  l$impair_Tmid_ <- 41
  l$impair_H_ <- 3e5
  drive_at(l, 41)
  expect_equal(l$phi_i_, 0.5, tolerance = 1e-12)
  drive_at(l, 44)
  expect_equal(l$phi_i_, phi_i_ref(44, H = 3e5, Tmid = 41), tolerance = 1e-12)
})

test_that("a non-positive impairment enthalpy is refused", {
  l <- thermal_leaf()
  l$impair_H_ <- -1
  expect_silent(drive_at(l, 30))   # gate off: not read
  l$use_thermal_damage_ <- TRUE
  expect_error(drive_at(l, 30), "needs impair_H_ > 0")
})

# ⚠️ THE DECISION THIS LAYER ENCODES. The Vcmax and Jmax curves are the
# instantaneous, reversible response and impairment is a separate process, so
# they share no parameter: moving one must not move the other.
test_that("the impairment onset is independent of the Jmax curve", {
  base <- thermal_leaf()
  base$use_thermal_damage_ <- TRUE
  drive_at(base, 38)

  jm <- thermal_leaf()
  jm$use_thermal_damage_ <- TRUE
  jm$jmax_H_d_ <- 210000
  jm$jmax_d_S_ <- 670
  drive_at(jm, 38)
  expect_false(isTRUE(all.equal(jm$jmax_, base$jmax_)))  # the curve moved
  expect_identical(jm$phi_i_, base$phi_i_)                # the onset did not

  im <- thermal_leaf()
  im$use_thermal_damage_ <- TRUE
  im$impair_Tmid_ <- 40
  im$impair_H_ <- 3e5
  drive_at(im, 38)
  expect_false(isTRUE(all.equal(im$phi_i_, base$phi_i_)))  # the onset moved
  expect_identical(c(im$vcmax_, im$jmax_), c(base$vcmax_, base$jmax_))
})

test_that("the leak into damage is a centred Arrhenius in leaf temperature", {
  # Default k_i_Ea_ = 0: exactly k_i_ at every temperature.
  for (temp in c(5, 25, 45)) {
    l <- thermal_leaf()
    l$use_thermal_damage_ <- TRUE
    drive_at(l, temp)
    expect_identical(l$k_i_eff_, l$k_i_, info = paste("leaf_temp", temp))
  }
  l <- thermal_leaf()
  l$use_thermal_damage_ <- TRUE
  l$k_i_ <- 0.3
  l$k_i_Ea_ <- 1.5e5
  l$k_i_Tref_ <- 35
  drive_at(l, 35)
  expect_identical(l$k_i_eff_, 0.3)   # exact at the reference
  drive_at(l, 42)
  expect_equal(l$k_i_eff_,
               0.3 * exp((1.5e5 / .R_gas) * (1 / (35 + .T0) - 1 / (42 + .T0))),
               tolerance = 1e-14)
  # The sharpness Ea sets: tenfold between T1 and the T2 with
  # 1/T2 = 1/T1 - R*ln(10)/Ea, about ln(10)*R*T^2/Ea degrees above it.
  T2 <- 1 / (1 / (35 + .T0) - .R_gas * log(10) / 1.5e5) - .T0
  drive_at(l, T2)
  expect_equal(l$k_i_eff_ / 0.3, 10, tolerance = 1e-12)
  # Bounded: the exponent is clamped at 30, so an absurd Ea stays finite.
  l$k_i_Ea_ <- 1e8
  drive_at(l, 90)
  expect_identical(l$k_i_eff_, 0.3 * exp(30))
})

test_that("acclimation moves the impairment onset by its own amount", {
  # Acclimated, with the onset's own maximum shift switched off: the capacity
  # optimum moves, the onset does not.
  base <- thermal_leaf()
  base$use_thermal_damage_ <- TRUE
  drive_at(base, 38)
  acc <- thermal_leaf()
  acc$use_thermal_damage_ <- TRUE
  acc$A_ <- 3
  acc$dTimpair_max_ <- 0
  drive_at(acc, 38)
  expect_false(isTRUE(all.equal(acc$vcmax_, base$vcmax_)))
  expect_equal(acc$phi_i_, base$phi_i_, tolerance = 1e-15)

  # And the other way round.
  only <- thermal_leaf()
  only$use_thermal_damage_ <- TRUE
  only$A_ <- 3
  only$dTopt_max_ <- 0
  drive_at(only, 38)
  expect_lt(only$phi_i_, base$phi_i_)
  expect_equal(only$vcmax_, base$vcmax_, tolerance = 1e-12)

  # At the defaults the two maxima are equal, so one shift moves both by the
  # same number of degrees, as in the port.
  expect_identical(base$dTimpair_max_, base$dTopt_max_)
})

test_that("lasting damage discounts BOTH capacities by 1 - I_r - I_p", {
  undamaged <- thermal_leaf()
  undamaged$use_thermal_damage_ <- TRUE
  drive_at(undamaged, 30)
  damaged <- thermal_leaf()
  damaged$use_thermal_damage_ <- TRUE
  damaged$I_r_ <- 0.3
  damaged$I_p_ <- 0.1
  drive_at(damaged, 30)
  expect_equal(damaged$vcmax_ / undamaged$vcmax_, 0.6, tolerance = 1e-14)
  expect_equal(damaged$jmax_ / undamaged$jmax_, 0.6, tolerance = 1e-14)

  # Past total damage the survival factor sits on its floor, so the capacities
  # stay strictly positive and the solve stays well-conditioned.
  crushed <- thermal_leaf()
  crushed$use_thermal_damage_ <- TRUE
  crushed$I_r_ <- 0.8
  crushed$I_p_ <- 0.5
  drive_at(crushed, 30)
  expect_equal(crushed$vcmax_ / undamaged$vcmax_, 1e-6, tolerance = 1e-12)
  expect_equal(crushed$jmax_ / undamaged$jmax_, 1e-6, tolerance = 1e-12)
  # A leaf this damaged cannot cover its respiration, so the solve takes the
  # no-flow exit a dark leaf takes -- and completes there, with finite carbon.
  # (lambda_emergent is NaN at that exit by design, as it is for a dark leaf.)
  crushed$find_root_collar_psi()
  expect_true(is.finite(crushed$profit_))
  expect_identical(crushed$assim_colimited_, -crushed$R_d_)
})

test_that("a thermostability offset moves the whole response up in temperature", {
  base <- thermal_leaf()
  base$use_thermal_damage_ <- TRUE
  drive_at(base, 40)
  tolerant <- thermal_leaf()
  tolerant$use_thermal_damage_ <- TRUE
  tolerant$topt_offset_ <- 8
  drive_at(tolerant, 40)
  expect_gt(tolerant$vcmax_, base$vcmax_)  # the capacity optimum moved up
  expect_lt(tolerant$phi_i_, base$phi_i_)  # and so did the impairment onset
})

test_that("the repair rate is gated off above t_rep_cut", {
  cold <- thermal_leaf()
  cold$use_thermal_damage_ <- TRUE
  drive_at(cold, 20)
  hot <- thermal_leaf()
  hot$use_thermal_damage_ <- TRUE
  drive_at(hot, 60)
  expect_equal(cold$k_rec_eff_, cold$k_rec_, tolerance = 1e-3)
  expect_lt(hot$k_rec_eff_, 1e-2 * hot$k_rec_)
  # At exactly t_rep_cut the logistic's argument is zero, so this is exact.
  cut <- thermal_leaf()
  cut$use_thermal_damage_ <- TRUE
  drive_at(cut, cut$t_rep_cut_)
  expect_identical(cut$k_rec_eff_, cut$k_rec_ / 2)
})

test_that("phi_i is finite, in [0, 1] and non-decreasing across a wide sweep", {
  temps <- seq(-30, 90, by = 1)
  phis <- vapply(temps, function(temp) {
    l <- thermal_leaf()
    l$use_thermal_damage_ <- TRUE
    drive_at(l, temp)
    l$phi_i_
  }, numeric(1))
  expect_true(all(is.finite(phis)))
  expect_true(all(phis >= 0 & phis <= 1))
  expect_true(all(diff(phis) >= 0))
})

test_that("outlandish inputs still give finite capacities and outputs", {
  l <- thermal_leaf()
  l$use_thermal_damage_ <- TRUE
  l$topt_offset_ <- 50
  l$I_r_ <- 0.4
  l$I_p_ <- 0.4
  drive_at(l, 99)
  l$k_i_Ea_ <- 5e5
  l$impair_Tmid_ <- 200
  drive_at(l, 99)
  expect_true(all(is.finite(c(l$vcmax_, l$jmax_, l$phi_i_, l$k_i_eff_,
                              l$k_rec_eff_, l$electron_transport_))))
  expect_true(l$phi_i_ >= 0 && l$phi_i_ <= 1)
})

# On the energy-balance path the block is re-seated at each candidate's leaf
# temperature, so phi_i_ follows the operating point: more transpiration, a
# cooler leaf, less impairment. That is the avoidance axis.
test_that("on the energy-balance path phi_i is taken at the operating Tleaf", {
  pm_leaf <- function(gate) {
    l <- thermal_leaf()
    l$use_energy_balance_ <- TRUE
    l$use_thermal_damage_ <- gate
    drive_at(l, 42, psi_soil = 0.2)
    l
  }
  lo <- pm_leaf(TRUE)
  lo$set_leaf_states_rates_from_psi_stem(1.0, 0.2)
  hi <- pm_leaf(TRUE)
  hi$set_leaf_states_rates_from_psi_stem(3.5, 0.2)
  expect_gt(hi$transpiration_, lo$transpiration_)
  expect_lt(hi$Tleaf_, lo$Tleaf_)
  expect_lt(hi$phi_i_, lo$phi_i_)
  # At THAT temperature, not at the air temperature driver.
  expect_equal(hi$phi_i_, phi_i_ref(hi$Tleaf_), tolerance = 1e-12)
  expect_equal(lo$phi_i_, phi_i_ref(lo$Tleaf_), tolerance = 1e-12)

  off <- pm_leaf(FALSE)
  off$set_leaf_states_rates_from_psi_stem(3.5, 0.2)
  expect_identical(off$phi_i_, 0)
})

# --- respiration: phylloptim's R_d, plus the costs --------------------------
#
# ⚠️ THE ONE INTENDED DIVERGENCE FROM THE JULY LAYER. There R_d was 0.015 *
# Vcmax after the damage discount, so a damaged leaf respired less. Here R_d is
# the R_d_25 trait on its Q10, for every model, and the layer only ADDS its
# costs: damage does not lower respiration.
test_that("R_d is not scaled by survival, and the costs add to it", {
  off <- drive_at(thermal_leaf(), 30)
  damaged <- thermal_leaf()
  damaged$use_thermal_damage_ <- TRUE
  damaged$I_r_ <- 0.3
  damaged$I_p_ <- 0.1
  drive_at(damaged, 30)
  # Coefficients zero: the four terms add an exact zero.
  expect_identical(damaged$R_d_, off$R_d_)

  cost_of <- function(...) {
    l <- thermal_leaf()
    l$use_thermal_damage_ <- TRUE
    set_inputs(l, list(...))
    drive_at(l, 30)
    l
  }
  acc <- cost_of(A_ = 3, c_acclim_maint_ = 0.1)
  expect_equal(acc$R_d_ - off$R_d_, 0.1 * 3, tolerance = 1e-12)
  # The standing repair cost is paid on the CAPACITY, even when cold.
  rep <- cost_of(c_repair_maint_ = 1e-2)
  expect_equal(rep$R_d_ - off$R_d_, 1e-2 * rep$k_rec_, tolerance = 1e-12)
  prot <- cost_of(k_i_ = 0.2, c_protect_maint_ = 0.5)
  expect_equal(prot$R_d_ - off$R_d_, 0.5 * (prot$k_i_baseline_ - 0.2),
               tolerance = 1e-12)
  # Protection above the baseline costs nothing extra: max(0, k_i_baseline - k_i).
  none <- cost_of(k_i_ = 0.9, c_protect_maint_ = 0.5)
  expect_identical(none$R_d_, off$R_d_)
  flux <- cost_of(I_r_ = 0.4, c_repair_flux_ = 1.0)
  expect_equal(flux$R_d_ - off$R_d_, flux$k_rec_eff_ * 0.4, tolerance = 1e-12)
  expect_gt(flux$R_d_, off$R_d_)

  # All four at once, with every input at a non-default value.
  all4 <- do.call(cost_of, atls_inputs)
  expected <- with(atls_inputs,
                   c_acclim_maint_ * A_ + c_repair_maint_ * k_rec_ +
                     c_protect_maint_ * (k_i_baseline_ - k_i_) +
                     c_repair_flux_ * all4$k_rec_eff_ * I_r_)
  expect_equal(all4$R_d_ - off$R_d_, expected, tolerance = 1e-12)
})

# --- off: inert, with every input set -----------------------------------------

test_that("with the gate off, no thermal-damage input is read", {
  # Every input at a value that would move the answer if it were read, against a
  # leaf that never had them set. identical(), on every reported output, across
  # both supply paths, both temperature paths and two temperatures.
  checked <- 0L
  for (supply in list(leaf_supply_singlelayer(), leaf_supply_multilayer())) {
    for (eb in c(FALSE, TRUE)) {
      for (temp in c(25, 40)) {
        clean <- thermal_leaf(supply)
        clean$use_energy_balance_ <- eb
        drive_at(clean, temp)
        clean$find_root_collar_psi()

        set <- set_inputs(thermal_leaf(supply))
        set$use_energy_balance_ <- eb
        drive_at(set, temp)
        set$find_root_collar_psi()

        at <- paste(supply$kind, "eb", eb, "T", temp)
        expect_identical(set$operating_point_values(),
                         clean$operating_point_values(), info = at)
        expect_identical(c(set$vcmax_, set$jmax_, set$R_d_),
                         c(clean$vcmax_, clean$jmax_, clean$R_d_), info = at)
        psi <- clean$opt_root_psi_
        expect_identical(set$dprofit_droot_collar_psi(psi),
                         clean$dprofit_droot_collar_psi(psi), info = at)
        expect_identical(c(set$phi_i_, set$k_i_eff_, set$k_rec_eff_),
                         c(0, 0, 0), info = at)

        # ...and the same inputs DO move it with the gate on, or the identity
        # above could be an inert gate rather than an inert layer.
        set$use_thermal_damage_ <- TRUE
        drive_at(set, temp)
        set$find_root_collar_psi()
        expect_false(identical(set$operating_point_values(),
                               clean$operating_point_values()), info = at)
        checked <- checked + 1L
      }
    }
  }
  expect_identical(checked, 8L)
})

# --- on: what it does, at non-default values ----------------------------------

test_that("the entropy shift moves each capacity optimum by the requested amount", {
  # Topt = H_d/(d_S - R*ln(Ea/(H_d - Ea))), so the argmax of vcmax_(T) and of
  # jmax_(T) is known in closed form; locate it numerically through the model.
  argmax_T <- function(l, field) {
    stats::optimize(function(temp) {
      drive_at(l, temp)
      l[[field]]
    }, c(15, 55), maximum = TRUE, tol = 1e-9)$maximum
  }
  topt <- function(ha, H_d, d_S) H_d / (d_S - .R_gas * log(ha / (H_d - ha))) - .T0

  base <- thermal_leaf()
  base$use_thermal_damage_ <- TRUE
  shifted <- set_inputs(thermal_leaf())
  shifted$use_thermal_damage_ <- TRUE
  # 2.5 + 4 * 1.5 / (1.5 + 2) = 4.214...
  delta <- with(atls_inputs, topt_offset_ + dTopt_max_ * A_ / (A_ + K_A_))

  for (field in c("vcmax_", "jmax_")) {
    pre <- if (field == "vcmax_") "vcmax" else "jmax"
    T0 <- topt(base[[paste0(pre, "_ha_")]], base[[paste0(pre, "_H_d_")]],
               base[[paste0(pre, "_d_S_")]])
    expect_equal(argmax_T(base, field), T0, tolerance = 1e-5, info = field)
    expect_equal(argmax_T(shifted, field) - argmax_T(base, field), delta,
                 tolerance = 1e-5, info = field)
  }

  # The survival factor is applied on top of the shifted curve.
  surv <- 1 - atls_inputs$I_r_ - atls_inputs$I_p_
  undamaged <- set_inputs(thermal_leaf(),
                          utils::modifyList(atls_inputs,
                                            list(I_r_ = 0, I_p_ = 0)))
  undamaged$use_thermal_damage_ <- TRUE
  drive_at(undamaged, 33)
  drive_at(shifted, 33)
  expect_equal(shifted$vcmax_ / undamaged$vcmax_, surv, tolerance = 1e-14)
  expect_equal(shifted$jmax_ / undamaged$jmax_, surv, tolerance = 1e-14)
})

test_that("the onset moves with its own shift, and repair halves at t_rep_cut", {
  l <- set_inputs(thermal_leaf())
  l$use_thermal_damage_ <- TRUE
  delta <- with(atls_inputs, topt_offset_ + dTimpair_max_ * A_ / (A_ + K_A_))
  # The midpoint, moved by the constitutive offset and the onset's own
  # acclimated shift.
  drive_at(l, l$impair_Tmid_ + delta)
  expect_equal(l$phi_i_, 0.5, tolerance = 1e-12)

  drive_at(l, l$t_rep_cut_)
  expect_identical(l$k_rec_eff_, l$k_rec_ / 2)
  drive_at(l, l$t_rep_cut_ + 3)
  expect_equal(l$k_rec_eff_,
               l$k_rec_ / (1 + exp(l$m_rep_ * 3)), tolerance = 1e-14)
})

test_that("the two switches are the guarded logistic and softplus", {
  l <- thermal_leaf()
  expect_identical(l$logistic_(0), 0.5)
  expect_equal(l$logistic_(1.3), 1 / (1 + exp(-1.3)), tolerance = 1e-15)
  # Clamped at |x| = 30, so neither end overflows.
  expect_identical(l$logistic_(1e6), 1 / (1 + exp(-30)))
  expect_identical(l$logistic_(-1e6), 1 / (1 + exp(30)))
  expect_equal(l$softplus_(0, 1), log(2), tolerance = 1e-15)
  expect_equal(l$softplus_(0.7, 2.5), log1p(exp(2.5 * 0.7)) / 2.5,
               tolerance = 1e-15)
  expect_true(is.finite(l$softplus_(1e6, 1)))
})

# --- the temperature cache ----------------------------------------------------
#
# ⚠️ THE INVARIANT A LATER EDIT COULD BREAK. The cache key does not hold the
# layer's inputs, because set_physiology bypasses the cache whenever the gate is
# on: A_, I_r_ and I_p_ are state a caller moves every step. Remove the bypass and
# the second drive below takes a cache hit and keeps the first damage level.
test_that("with the gate on, every re-drive reads the current state", {
  l <- thermal_leaf()
  l$use_thermal_damage_ <- TRUE
  drive_at(l, 30)
  l$find_root_collar_psi()
  vcmax_first <- l$vcmax_
  A_first <- l$assim_colimited_

  l$I_r_ <- 0.25                      # same drivers, new damage state
  drive_at(l, 30)
  l$find_root_collar_psi()
  expect_equal(l$vcmax_ / vcmax_first, 0.75, tolerance = 1e-14)
  expect_lt(l$assim_colimited_, A_first)

  l$A_ <- 2                           # and new acclimation state
  vcmax_damaged <- l$vcmax_
  drive_at(l, 30)
  expect_false(identical(l$vcmax_, vcmax_damaged))

  # Off again at the same drivers: the gate is in the key, so this recomputes
  # rather than serving the gate-on block, and matches a leaf never switched on.
  l$use_thermal_damage_ <- FALSE
  drive_at(l, 30)
  l$find_root_collar_psi()
  never <- drive_at(thermal_leaf(), 30)
  never$find_root_collar_psi()
  expect_identical(l$operating_point_values(), never$operating_point_values())
  expect_identical(c(l$phi_i_, l$k_i_eff_, l$k_rec_eff_), c(0, 0, 0))
})

# --- refusals -----------------------------------------------------------------

test_that("the trait gradient refuses a leaf with the layer on", {
  l <- thermal_leaf(leaf_supply_multilayer())
  l$use_thermal_damage_ <- TRUE
  expect_error(leaf_gradient(psi_soil = 2, x = l, traits = leaf_traits(),
                             pars = "vcmax_25"),
               "leaf_gradient\\(\\): refused with use_thermal_damage_ on")
  # The refusal comes before the leaf is touched.
  expect_true(l$use_thermal_damage_)

  b <- leaf_batch(psi_soil = c(1, 2), PPFD = 900)
  b$leaf$use_thermal_damage_ <- TRUE
  expect_error(leaf_gradient_batch(b, pars = "vcmax_25"),
               "leaf_gradient_batch\\(\\): refused with use_thermal_damage_ on")
  # And the batch runs again once the gate is off.
  b$leaf$use_thermal_damage_ <- FALSE
  g <- leaf_gradient_batch(b, pars = "vcmax_25")
  expect_true(all(g$status %in% c("interior", "pinned")))
})

test_that("a temperature response with no optimum to move is refused", {
  # The activation and deactivation energies are settable, and the peaked
  # Arrhenius has an optimum only for Ea < H_d. With the gate off nothing
  # moves the optimum, so nothing needs one.
  l <- thermal_leaf()
  l$vcmax_ha_ <- 2.5e5
  expect_silent(drive_at(l, 30))
  l$use_thermal_damage_ <- TRUE
  expect_error(drive_at(l, 30), "exists only for 0 < Ea < H_d")
})

test_that("the closed form refuses a leaf with the layer on", {
  l <- thermal_leaf()
  drive_at(l, 25)
  l$set_model("TF24", "stem", "closed")
  l$use_thermal_damage_ <- TRUE
  expect_error(l$optimise(), "closed form is refused with use_thermal_damage_")
  # The exact method on the same leaf solves.
  l$set_model("TF24", "stem")
  l$optimise()
  expect_true(is.finite(l$profit_))
})
