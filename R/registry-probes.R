# How registry_check() reads each equation's quantity off the package.
#
# One probe per equation id in inst/model-registry.yml. A probe takes the
# replicas (an environment from registry_replicas(), so `rep$arrhenius(...)`
# is the registry's own function) and the temperature grid, drives a Leaf, and
# returns a data.frame with columns `pkg` (what the package reports) and `mine`
# (what the replica gives for the same inputs). registry_check() does the
# comparing; a probe only pairs numbers.
#
# Two things a probe must do: read the package through a binding, never
# re-derive the quantity it is checking (that is the replica's job), and reach
# the branches the equation has (a clamp, a floor, a gate off and on), or the
# check says nothing about them. Where a quantity has no binding of its own it
# is read through the one that exposes it, and the entry's note says so.
#
# ⚠️ The temperature block is cached on every scalar it reads (photo_temp_key),
# so re-driving after writing a temperature parameter does recompute; with the
# thermal-damage gate on the cache is bypassed altogether. A trait still goes
# through set_traits() (hazard 10 of the developer guide).

.probe_pairs <- function(pkg, mine) {
  data.frame(pkg = as.numeric(pkg), mine = as.numeric(mine))
}

.probe_fields <- function(l, nm) stats::setNames(lapply(nm, function(n) l[[n]]), nm)

.probe_drive <- function(l, T, PPFD = 900, atm_vpd = 2) {
  set_drivers(l, psi_soil = 1, PPFD = PPFD, atm_vpd = atm_vpd, leaf_temp = T)
  invisible(l)
}

# A single-potential leaf with the energy balance on, driven at daylight NIR.
.probe_eb <- function(PPFD = 900, atm_vpd = 2, Tair = 25, d = 0.05, U = 2,
                      NIR = PPFD / 4.57, alpha_NIR = NULL) {
  e <- leaf_model(supply = leaf_supply_singlelayer())
  e$use_energy_balance_ <- TRUE
  e$d_ <- d
  e$wind_speed_ <- U
  e$NIR_int_ <- NIR
  if (!is.null(alpha_NIR)) e$alpha_NIR_ <- alpha_NIR
  set_drivers(e, psi_soil = 1, PPFD = PPFD, atm_vpd = atm_vpd, leaf_temp = Tair,
              root_network = series_resistance(1e3))
  e
}

.probe_eb_grid <- expand.grid(d = c(0.02, 0.05, 0.2), U = c(0.5, 2, 5),
                              Q = c(300, 900, 1800), Ta = c(5, 10, 25, 40, 50),
                              D = c(1, 2, 4))

.probe_temp_fields <- c("vcmax_25", "vcmax_ha_", "vcmax_H_d_", "vcmax_d_S_",
                        "jmax_25", "jmax_ha_", "jmax_H_d_", "jmax_d_S_",
                        "gamma_25_", "gamma_ha_", "kc_25_", "kc_ha_", "ko_25_", "ko_ha_")

.probe_damage_fields <- c("topt_offset_", "impair_H_", "impair_Tmid_", "k_i_", "k_i_Ea_",
                          "k_i_Tref_", "k_rec_", "m_rep_", "t_rep_cut_", "k_mat_",
                          "acclim_alpha_", "acclim_beta_", "acclim_T_", "acclim_s_",
                          "dTimpair_max_", "K_A_", "dTopt_max_", "k_i_baseline_",
                          "c_acclim_maint_", "c_repair_maint_", "c_protect_maint_",
                          "c_repair_flux_")

# The layer's inputs, moved cumulatively: each set is applied on top of the
# last, so the fourth reaches the functional-fraction floor at the third's
# kinetics. Negative A, every clamp and the floor are all reached.
.probe_damage_sets <- list(
  list(A_ = 0, I_ = 0, D_ = 0, topt_offset_ = 0),
  list(A_ = 1, I_ = 0.1, D_ = 0.05, topt_offset_ = 2, k_i_Ea_ = 1.5e5, k_i_Tref_ = 35),
  list(A_ = 4, I_ = 0.3, D_ = 0.2, topt_offset_ = 0, dTopt_max_ = 10, dTimpair_max_ = 3,
       K_A_ = 2, impair_H_ = 3e5, impair_Tmid_ = 41, t_rep_cut_ = 40, m_rep_ = 0.8,
       k_i_Ea_ = 5e4, k_mat_ = 0.05, acclim_alpha_ = 0.1, acclim_beta_ = 0.3,
       acclim_T_ = 28, acclim_s_ = 0.7, c_acclim_maint_ = 0.05, c_repair_maint_ = 0.3,
       c_protect_maint_ = 0.4, c_repair_flux_ = 0.8, k_i_ = 0.2),
  list(A_ = -0.5, I_ = 0.7, D_ = 0.4))

# Run `f(d, T, p)` on a damage leaf at each cumulative set and temperature; `p`
# is the layer's fields as a named list. `f` returns a pairs data.frame.
.probe_damage <- function(temps, f) {
  d <- leaf_model()
  d$use_thermal_damage_ <- TRUE
  out <- list()
  for (st in .probe_damage_sets) {
    for (nm in names(st)) d[[nm]] <- st[[nm]]
    p <- .probe_fields(d, .probe_damage_fields)
    for (T in temps) {
      .probe_drive(d, T)
      out[[length(out) + 1L]] <- f(d, T, p)
    }
  }
  do.call(rbind, out)
}

# The Rn_ inversion: leaf_temp_from_E(0) is affine in Rn_, so two calls give
# the slope and Rn_ can be set to put the zero-transpiration temperature at a
# chosen Tleaf; the zero-flow branch of set_leaf_states_rates_from_psi_stem
# then seats the block and the deficit there. Override Rn_ AFTER set_drivers,
# which rewrites it.
.probe_place_tleaf <- function(e, Tl) {
  e$Rn_ <- 0; T0 <- e$leaf_temp_from_E(0)
  e$Rn_ <- 1; T1 <- e$leaf_temp_from_E(0)
  e$Rn_ <- (Tl - T0) / (T1 - T0)
  e$set_leaf_states_rates_from_psi_stem(1, 1)
  invisible(e)
}

.registry_probes <- list(

  # --- biochemistry ----------------------------------------------------------
  arrhenius = function(rep, temps) {
    l <- leaf_model()
    g <- expand.grid(T = temps, Ea = c(30000, 60000, 90000))
    .probe_pairs(mapply(function(T, Ea) l$arrh_curve(Ea, 42.75, T), g$T, g$Ea),
                 rep$arrhenius(g$T, g$Ea, 42.75))
  },
  peaked_arrhenius = function(rep, temps) {
    l <- leaf_model()
    g <- expand.grid(T = temps, H_d = c(150000, 200000, 250000), d_S = c(600, 650, 700))
    .probe_pairs(mapply(function(T, H_d, d_S) l$peak_arrh_curve(60000, 96, T, H_d, d_S),
                        g$T, g$H_d, g$d_S),
                 rep$peaked_arrhenius(g$T, 60000, 96, g$H_d, g$d_S))
  },
  # The identity: the package's own peaked curve is stationary at the replica's
  # optimum. The slope is a central difference scaled by the peak.
  thermal_optimum = function(rep, temps) {
    l <- leaf_model()
    sets <- list(c(60000, 200000, 650), c(30000, 200000, 650),
                 c(80000, 200000, 640), c(40000, 150000, 600))
    h <- 1e-3
    slopes <- vapply(sets, function(s) {
      Topt <- rep$thermal_optimum(s[1], s[2], s[3])
      f <- function(T) l$peak_arrh_curve(s[1], 1, T, s[2], s[3])
      (f(Topt + h) - f(Topt - h)) / (2 * h) / f(Topt)
    }, 0)
    .probe_pairs(slopes, 0)
  },
  temperature_block = function(rep, temps) {
    l <- leaf_model()
    p <- .probe_fields(l, .probe_temp_fields)
    do.call(rbind, lapply(temps, function(T) {
      .probe_drive(l, T)
      r <- rep$temperature_block(T, p)
      .probe_pairs(c(l$vcmax_, l$jmax_, l$gamma_, l$kc_, l$ko_),
                   r[c("vcmax", "jmax", "gamma", "kc", "ko")])
    }))
  },
  michaelis_constant = function(rep, temps) {
    l <- leaf_model()
    do.call(rbind, lapply(temps, function(T) {
      .probe_drive(l, T)
      .probe_pairs(l$km_, rep$michaelis_constant(l$kc_, l$ko_, l$atm_o2_kpa_, l$atm_kpa_))
    }))
  },
  electron_transport = function(rep, temps) {
    l <- leaf_model()
    g <- expand.grid(T = temps, Q = c(100, 900, 1800))
    do.call(rbind, mapply(function(T, Q) {
      .probe_drive(l, T, PPFD = Q)
      .probe_pairs(l$electron_transport_,
                   rep$electron_transport(l$PPFD_, l$phi_J, l$alpha_PAR_, l$jmax_,
                                          l$curv_fact_elec_trans))
    }, g$T, g$Q, SIMPLIFY = FALSE))
  },
  assimilation = function(rep, temps) {
    l <- leaf_model()
    g <- expand.grid(T = temps, ci = c(10, 20, 30, 38))
    do.call(rbind, mapply(function(T, ci) {
      .probe_drive(l, T)
      .probe_pairs(l$assim_colimited(ci),
                   rep$assimilation(ci, l$vcmax_, l$electron_transport_, l$gamma_, l$km_,
                                    l$curv_fact_colim, l$R_d_, l$atm_kpa_))
    }, g$T, g$ci, SIMPLIFY = FALSE))
  },
  thermal_cost = function(rep, temps) {
    on <- leaf_model(); on$use_thermal_cost_ <- TRUE
    moved <- leaf_model(); moved$use_thermal_cost_ <- TRUE; moved$Tcrit_ <- 40
    set_traits(moved, leaf_traits())
    off <- leaf_model()
    do.call(rbind, lapply(temps, function(T) {
      .probe_drive(on, T); .probe_drive(moved, T); .probe_drive(off, T)
      jmax <- function(l) rep$peaked_arrhenius(T, l$jmax_ha_, l$jmax_25, l$jmax_H_d_, l$jmax_d_S_)
      .probe_pairs(c(on$thermal_cost_at(T), on$jmax_, moved$thermal_cost_at(T), moved$jmax_,
                     off$thermal_cost_at(T), off$jmax_),
                   c(rep$thermal_cost(T, on$T50_, on$Tcrit_),
                     jmax(on) * (1 - rep$thermal_cost(T, on$T50_, on$Tcrit_)),
                     rep$thermal_cost(T, moved$T50_, moved$Tcrit_),
                     jmax(moved) * (1 - rep$thermal_cost(T, moved$T50_, moved$Tcrit_)),
                     0, jmax(off)))
    }))
  },

  # --- respiration ------------------------------------------------------------
  # At the driver temperature with the gate off, and at a placed leaf
  # temperature with it on, since the block is re-seated there.
  q10_respiration = function(rep, temps) {
    l <- leaf_model()
    e <- .probe_eb(Tair = 25, atm_vpd = 1)
    do.call(rbind, lapply(temps, function(T) {
      .probe_drive(l, T)
      .probe_place_tleaf(e, T)
      .probe_pairs(c(l$R_d_, e$R_d_),
                   c(rep$q10_respiration(T, l$R_d_25, l$rd_q10_intercept_, l$rd_q10_slope_),
                     rep$q10_respiration(e$Tleaf_, e$R_d_25, e$rd_q10_intercept_, e$rd_q10_slope_)))
    }))
  },

  # --- water vapour -----------------------------------------------------------
  saturation_vapour_pressure = function(rep, temps) {
    l <- leaf_model()
    .probe_pairs(vapply(temps, l$saturation_vapour_pressure, 0),
                 rep$saturation_vapour_pressure(temps))
  },
  saturation_vapour_pressure_slope = function(rep, temps) {
    l <- leaf_model()
    .probe_pairs(vapply(temps, l$saturation_vapour_pressure_slope, 0),
                 rep$saturation_vapour_pressure_slope(temps))
  },
  # Placed at a chosen leaf temperature with the gate on (the floor binds for a
  # cold leaf under a dry air); exactly the air's deficit with the gate off.
  leaf_vpd = function(rep, temps) {
    l <- leaf_model()
    g <- expand.grid(Ta = c(15, 25, 35), D = c(0.5, 1, 2, 4))
    out <- lapply(seq_len(nrow(g)), function(i) {
      e <- .probe_eb(Tair = g$Ta[i], atm_vpd = g$D[i])
      do.call(rbind, lapply(temps, function(Tl) {
        .probe_place_tleaf(e, Tl)
        .probe_pairs(c(e$Tleaf_, e$vpd_leaf_),
                     c(Tl, rep$leaf_vpd(e$Tleaf_, g$Ta[i], g$D[i])))
      }))
    })
    off <- do.call(rbind, lapply(temps, function(T) {
      .probe_drive(l, T, atm_vpd = 1.5)
      .probe_pairs(l$vpd_leaf_, rep$leaf_vpd(T, T, 1.5, energy_balance = FALSE))
    }))
    rbind(do.call(rbind, out), off)
  },
  # The transpiring branch at a prescribed stem potential: the reported
  # conductance against Fick's law on the reported transpiration and deficit,
  # with the gate off (the air's deficit) and on (the leaf's, at the leaf
  # temperature that transpiration gives).
  stomatal_conductance = function(rep, temps) {
    off <- leaf_model(supply = leaf_supply_singlelayer())
    set_drivers(off, psi_soil = 1, PPFD = 900, atm_vpd = 2, leaf_temp = 25,
                root_network = series_resistance(1e3))
    g <- expand.grid(psi = c(1.2, 1.5, 2, 3), Ta = c(15, 25, 35))
    on <- lapply(seq_len(nrow(g)), function(i) {
      e <- .probe_eb(Tair = g$Ta[i])
      e$set_leaf_states_rates_from_psi_stem(g$psi[i], 1)
      .probe_pairs(c(e$stom_cond_CO2_, e$Tleaf_, e$vpd_leaf_),
                   c(rep$stomatal_conductance(e$transpiration_, e$atm_kpa_, e$vpd_leaf_,
                                              e$H2O_CO2_stom_diff_ratio_),
                     e$leaf_temp_from_E(e$transpiration_),
                     rep$leaf_vpd(e$Tleaf_, g$Ta[i], 2)))
    })
    gate_off <- do.call(rbind, lapply(c(1.2, 1.5, 2, 3), function(psi) {
      off$set_leaf_states_rates_from_psi_stem(psi, 1)
      .probe_pairs(off$stom_cond_CO2_,
                   rep$stomatal_conductance(off$transpiration_, off$atm_kpa_, off$vpd_leaf_,
                                            off$H2O_CO2_stom_diff_ratio_))
    }))
    rbind(do.call(rbind, on), gate_off)
  },

  # --- energy balance ----------------------------------------------------------
  # Not reported on its own: read back from Rn_, which is the absorbed shortwave
  # less (1 - emissivity) times the blackbody term.
  atmospheric_emissivity = function(rep, temps) {
    g <- expand.grid(Ta = c(5, 10, 15, 20, 25, 35, 40, 45, 50), D = c(0.5, 1, 2, 4))
    do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
      e <- .probe_eb(PPFD = 900, Tair = g$Ta[i], atm_vpd = g$D[i])
      TK <- g$Ta[i] + 273.15
      sw <- e$alpha_PAR_ * e$PPFD_ / 4.57 + e$alpha_NIR_ * e$NIR_int_
      .probe_pairs(1 - (sw - e$Rn_) / (5.670374419e-8 * TK^4),
                   rep$atmospheric_emissivity(g$Ta[i], g$D[i]))
    }))
  },
  net_radiation = function(rep, temps) {
    g <- .probe_eb_grid
    main <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
      e <- .probe_eb(PPFD = g$Q[i], Tair = g$Ta[i], atm_vpd = g$D[i], d = g$d[i], U = g$U[i])
      .probe_pairs(e$Rn_, rep$net_radiation(g$Q[i], g$Q[i] / 4.57, g$Ta[i], g$D[i],
                                            e$alpha_PAR_, e$alpha_NIR_))
    }))
    # The absorptances and NIR are read, not folded in.
    led <- .probe_eb(PPFD = 900, Tair = 25, NIR = 0, alpha_NIR = 0.35)
    extra <- .probe_pairs(led$Rn_, rep$net_radiation(900, 0, 25, 2, led$alpha_PAR_, led$alpha_NIR_))
    rbind(main, extra)
  },
  radiative_conductance = function(rep, temps) {
    do.call(rbind, lapply(c(5, 10, 25, 40, 50), function(Ta) {
      e <- .probe_eb(Tair = Ta)
      .probe_pairs(e$g_rad_, rep$radiative_conductance(Ta))
    }))
  },
  aerodynamic_resistance = function(rep, temps) {
    g <- expand.grid(d = c(0, 0.02, 0.05, 0.2), U = c(0, 0.5, 2, 5))
    do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
      e <- .probe_eb(d = g$d[i], U = g$U[i])
      .probe_pairs(e$ra_, rep$aerodynamic_resistance(g$d[i], g$U[i]))
    }))
  },
  # Over the driver grid at five transpirations, the last of which reaches the
  # lower clamp.
  leaf_temperature = function(rep, temps) {
    g <- .probe_eb_grid
    E <- c(0, 2e-5, 1e-4, 5e-4, 5e-3)
    do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
      e <- .probe_eb(PPFD = g$Q[i], Tair = g$Ta[i], atm_vpd = g$D[i], d = g$d[i], U = g$U[i])
      .probe_pairs(vapply(E, e$leaf_temp_from_E, 0),
                   rep$leaf_temperature(g$Ta[i], e$Rn_, E, e$ra_, e$g_rad_))
    }))
  },
  # The bound stem-route derivative, which carries the term, against a central
  # difference of the profit it differentiates. Only feasible points compare.
  energy_balance_foc = function(rep, temps) {
    out <- list()
    for (Ta in c(20, 30, 40)) {
      e <- .probe_eb(PPFD = 900, Tair = Ta, atm_vpd = 2)
      e$set_model("TF24", "stem", "exact")
      for (psi in c(1.3, 1.8, 2.5, 3.5)) {
        v <- e$dprofit_dpsi_stem_checked(psi)
        if (v[2] != 1) next
        out[[length(out) + 1L]] <- .probe_pairs(
          v[1], rep$energy_balance_foc(function(x) e$evaluate_psi_stem_at(x), psi, 1e-4))
      }
    }
    do.call(rbind, out)
  },

  # --- damage, repair and acclimation ---------------------------------------------
  # The shift and the shifted entropy term reach the package only through the
  # capacities, so both are checked as the capacity the leaf reports, divided
  # by the functional fraction it reports: J_max here, V_cmax in entropy_shift.
  acclimation_shift = function(rep, temps) {
    .probe_damage(temps, function(d, T, p) {
      delta <- rep$acclimation_shift(d$A_, p$topt_offset_, p$dTopt_max_, p$K_A_)
      jdS <- rep$entropy_shift(d$jmax_ha_, d$jmax_H_d_, d$jmax_d_S_, delta)
      .probe_pairs(d$jmax_ / d$functional_,
                   rep$peaked_arrhenius(T, d$jmax_ha_, d$jmax_25, d$jmax_H_d_, jdS))
    })
  },
  entropy_shift = function(rep, temps) {
    .probe_damage(temps, function(d, T, p) {
      delta <- rep$acclimation_shift(d$A_, p$topt_offset_, p$dTopt_max_, p$K_A_)
      vdS <- rep$entropy_shift(d$vcmax_ha_, d$vcmax_H_d_, d$vcmax_d_S_, delta)
      .probe_pairs(d$vcmax_ / d$functional_,
                   rep$peaked_arrhenius(T, d$vcmax_ha_, d$vcmax_25, d$vcmax_H_d_, vdS))
    })
  },
  functional_fraction = function(rep, temps) {
    .probe_damage(temps, function(d, T, p) {
      F <- rep$functional_fraction(d$I_, d$D_)
      .probe_pairs(c(d$functional_, d$vcmax_, d$jmax_),
                   c(F,
                     rep$peaked_arrhenius(T, d$vcmax_ha_, d$vcmax_25, d$vcmax_H_d_,
                       rep$entropy_shift(d$vcmax_ha_, d$vcmax_H_d_, d$vcmax_d_S_,
                         rep$acclimation_shift(d$A_, p$topt_offset_, p$dTopt_max_, p$K_A_))) * F,
                     rep$peaked_arrhenius(T, d$jmax_ha_, d$jmax_25, d$jmax_H_d_,
                       rep$entropy_shift(d$jmax_ha_, d$jmax_H_d_, d$jmax_d_S_,
                         rep$acclimation_shift(d$A_, p$topt_offset_, p$dTopt_max_, p$K_A_))) * F))
    })
  },
  respiration_costs = function(rep, temps) {
    .probe_damage(temps, function(d, T, p) {
      .probe_pairs(d$R_d_,
                   rep$q10_respiration(T, d$R_d_25, d$rd_q10_intercept_, d$rd_q10_slope_) +
                     rep$respiration_costs(d$A_, d$I_, d$k_rec_eff_, p))
    })
  },
  logistic = function(rep, temps) {
    d <- leaf_model()
    x <- seq(-40, 40, by = 5)
    .probe_pairs(vapply(x, d$logistic_, 0), rep$logistic(x))
  },
  softplus = function(rep, temps) {
    d <- leaf_model()
    g <- expand.grid(x = seq(-40, 40, by = 5), s = c(0.5, 1, 2))
    .probe_pairs(mapply(function(x, s) d$softplus_(x, s), g$x, g$s), rep$softplus(g$x, g$s))
  },
  # The three pieces of the impairment rate are reported only as their product,
  # k_imp_eff_. Each is checked through it: the shift at A > 0 (sets 2 and 3),
  # the switch alone at the set where the Arrhenius scale is exactly k_i, and
  # the scale as the product divided by the replica's switch.
  delta_impair = function(rep, temps) {
    .probe_damage(temps, function(d, T, p) {
      if (d$A_ <= 0) return(NULL)
      .probe_pairs(d$k_imp_eff_, rep$impairment_rate(T, d$A_, p))
    })
  },
  impair_switch = function(rep, temps) {
    .probe_damage(temps, function(d, T, p) {
      if (p$k_i_Ea_ != 0) return(NULL)
      delta <- rep$delta_impair(d$A_, p$topt_offset_, p$dTimpair_max_, p$K_A_)
      .probe_pairs(d$k_imp_eff_ / d$k_i_, rep$impair_switch(T, p$impair_H_, p$impair_Tmid_, delta))
    })
  },
  k_i_eff = function(rep, temps) {
    .probe_damage(temps, function(d, T, p) {
      delta <- rep$delta_impair(d$A_, p$topt_offset_, p$dTimpair_max_, p$K_A_)
      sw <- rep$impair_switch(T, p$impair_H_, p$impair_Tmid_, delta)
      .probe_pairs(d$k_imp_eff_ / sw, rep$k_i_eff(T, p$k_i_, p$k_i_Ea_, p$k_i_Tref_))
    })
  },
  impairment_rate = function(rep, temps) {
    .probe_damage(temps, function(d, T, p) {
      .probe_pairs(d$k_imp_eff_, rep$impairment_rate(T, d$A_, p))
    })
  },
  k_rec_eff = function(rep, temps) {
    .probe_damage(temps, function(d, T, p) {
      .probe_pairs(d$k_rec_eff_, rep$k_rec_eff(T, p$k_rec_, p$m_rep_, p$t_rep_cut_))
    })
  },
  damage_rates = function(rep, temps) {
    g <- expand.grid(A = c(-0.5, 0, 2), I = c(0, 0.25), D = c(0, 0.1))
    .probe_damage(temps, function(d, T, p) {
      do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
        .probe_pairs(d$damage_rates(T, g$A[i], g$I[i], g$D[i]),
                     rep$damage_rates(T, g$A[i], g$I[i], g$D[i], p))
      }))
    })
  },
  # A heatwave and its recovery, at each cumulative parameter set.
  damage_trajectory = function(rep, temps) {
    times <- c(0, 0.3, 0.5, 1, 2, 5, 10, 20)
    series <- c(25, 46, 46, 30, 25, 40, 25, 25)
    d <- leaf_model()
    out <- list()
    for (st in .probe_damage_sets) {
      for (nm in names(st)) d[[nm]] <- st[[nm]]
      p <- .probe_fields(d, .probe_damage_fields)
      out[[length(out) + 1L]] <- .probe_pairs(
        d$damage_trajectory(times, series, 0.5, 0.05, 0.01),
        as.vector(rep$damage_trajectory(times, series, c(0.5, 0.05, 0.01), p)))
    }
    do.call(rbind, out)
  }
)
