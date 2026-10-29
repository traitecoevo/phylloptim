# registry_check()'s probes for the Medlyn route, the soil-moisture price and the closed form.
#
# Same contract as R/registry-probes.R: one probe per equation id in
# inst/model-registry.yml, each driving a Leaf and pairing what the package
# reports with the registry's replica. registry_probes() collects every
# `.registry_probes_*` list in the namespace, so this file needs no other wiring.

# Probes for the medlyn fragment: the empirical Medlyn route, the CF77 soil-moisture
# factor and the closed-form method. Sourced with the phylloptim namespace as
# parent, so .probe_pairs, .probe_drive, leaf_model, set_drivers ... are in scope.
#
# The closed form's internals are not bound. Each closed-form probe runs
# optimise() on a leaf seated with set_model(curve, "stem", "closed") and reads
# what it reports: `opt_psi_stem_` (the closed form's psi* on a served row) and
# `last_solve_fell_back_` (the guard's verdict). The replicas take the leaf's
# state fields and the bound supply-side functions (lambda, transpiration and its
# inverse, the conductivity fraction) as arguments; nothing checked is
# re-derived on the package side.

# Four temperatures out of the check's grid, for the probes that solve.
.mf_T <- function(temps) temps[unique(round(seq(1, length(temps), length.out = 4)))]

# The leaf's state as the replicas read it.
.mf_p <- function(l) {
  list(ca = l$ca_, D = l$vpd_leaf_, ratio = l$H2O_CO2_stom_diff_ratio_, atm_kpa = l$atm_kpa_,
       vcmax = l$vcmax_, J = l$electron_transport_, gamma = l$gamma_, km = l$km_,
       theta_c = l$curv_fact_colim, R_d = l$R_d_, scale = l$TF24_cost_scale,
       beta2 = l$TF24_beta2, c = l$stem_c, b = l$stem_b, kmax = l$leaf_specific_conductance_max_,
       g0 = l$g0, g1 = l$g1, theta = l$theta_, theta_fc = l$theta_fc_, theta_w = l$theta_w_)
}

# The bound supply side, measured from psi_soil (the stem route's upstream).
.mf_h <- function(l, lambda = NULL) {
  up <- l$single_psi_soil_
  list(psi_up = up, psi_crit = l$psi_crit, E_crit = l$transpiration(l$psi_crit, up),
       transpiration = function(psi) l$transpiration(psi, up),
       inverse = function(E) l$transpiration_to_psi_stem(E, up),
       lambda = lambda, f = function(psi) l$proportion_of_conductivity(psi))
}

# A single-layer stem-route leaf seated closed; `set` is applied before the drivers.
.mf_leaf <- function(curve, psi_soil = 1, PPFD = 1500, D = 1, T = 25, lambda = 1.5e5,
                     lambda_o = 5e4, set = list(), traits = NULL) {
  l <- leaf_model(supply = leaf_supply_singlelayer())
  if (!is.null(traits)) set_traits(l, traits)
  l$CF77_lambda_ <- lambda
  l$TF24_floor_lambda_o <- lambda_o
  for (nm in names(set)) l[[nm]] <- set[[nm]]
  set_drivers(l, psi_soil = psi_soil, PPFD = PPFD, atm_vpd = D, leaf_temp = T,
              root_network = series_resistance(1500))
  l$set_model(curve, "stem", "closed")
  l
}

# Solve closed and run the replica arm. Returns psi (pkg, mine), the fall-back flag
# and the replica's own solution.
.mf_arm <- function(rep, curve, ...) {
  l <- .mf_leaf(curve, ...)
  l$optimise()
  p <- .mf_p(l)
  s <- switch(curve,
    CF77 = rep$closed_form_CF77(rep$cf77_soil_price(l$CF77_lambda_, l$CF77_soil_beta_, l$theta_,
                                                    l$theta_fc_, l$theta_w_), p, .mf_h(l)),
    TF24 = rep$closed_form_TF24(p, .mf_h(l, function(x) l$lambda_TF24(x))),
    TF24_floor = rep$closed_form_TF24_floor(l$TF24_floor_lambda_o, p,
                                            .mf_h(l, function(x) l$lambda_TF24_floor(x))))
  list(pkg = l$opt_psi_stem_, fb = l$last_solve_fell_back_, s = s, l = l,
       guard = rep$within_guard(s$psi, s$ci, s$A, l$ca_, l$single_psi_soil_, l$psi_crit))
}

# Served rows' psi*; a grid row that falls back contributes nothing here.
.mf_served <- function(rep, curve, g) {
  out <- lapply(seq_len(nrow(g)), function(i) {
    a <- lapply(as.list(g[i, , drop = FALSE]), function(x) if (is.list(x)) x[[1]] else x)
    r <- do.call(.mf_arm, c(list(rep = rep, curve = curve), a))
    if (r$fb) NULL else .probe_pairs(r$pkg, r$s$psi)
  })
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) < 5) stop("closed-form probe served fewer than 5 rows")
  out
}

# The residual of the Medlyn numerical solve, at the replica.
.mf_medlyn_leaf <- function(T, PPFD = 900, D = 2, theta = 0.3, g0 = 0.022, g1 = 2.57, R_d_25 = NULL) {
  l <- leaf_model(control = leaf_control(ci_abs_tol = 1e-12))
  if (!is.null(R_d_25)) set_traits(l, leaf_traits(R_d_25 = R_d_25))
  l$theta <- theta; l$g0 <- g0; l$g1 <- g1
  set_drivers(l, psi_soil = 1, PPFD = PPFD, atm_vpd = D, leaf_temp = T)
  l
}

.registry_probes_medlyn <- list(

  # --- medlyn -----------------------------------------------------------------
  # Called at a grid of A through the binding, with beta below 0, at 0, inside,
  # at 1 and above 1 (unclamped), and g0 from 0.
  medlyn_gs = function(rep, temps) {
    g <- expand.grid(A = c(-2, 0, 5, 15, 30), D = c(0.25, 1, 2, 4), theta = c(0.1, 0.2, 0.3, 0.5, 0.8),
                     g0 = c(0, 0.022, 0.1), g1 = c(1, 2.57, 6))
    l <- leaf_model()
    out <- lapply(split(g, list(g$D, g$theta, g$g0, g$g1), drop = TRUE), function(s) {
      l$theta <- s$theta[1]; l$g0 <- s$g0[1]; l$g1 <- s$g1[1]
      set_drivers(l, psi_soil = 1, atm_vpd = s$D[1])
      pkg <- vapply(s$A, l$medlyn_model_gs, 0)
      .probe_pairs(c(pkg, l$medlyn_model_gs_),
                   c(rep$medlyn_gs(s$A, l$g0, l$g1, l$atm_vpd_, l$ca_, l$atm_kpa_, l$theta_,
                                   l$theta_fc_, l$theta_w_), pkg[length(pkg)]))
    })
    do.call(rbind, out)
  },
  # The root branch: the package's ci against one Newton step on the replica
  # residual from it, then A and gs at it. The dark branch (no sign change, a
  # stop) and the R_d = 0, beta < 0 branch (the peak is <= 0, ci is the left end)
  # as indicators.
  medlyn_ci_numerical = function(rep, temps) {
    g <- expand.grid(T = .mf_T(temps), D = c(0.3, 1, 2, 4), theta = c(0.2001, 0.21, 0.3, 0.5, 0.7),
                     g0 = c(0, 0.022), Q = c(50, 900))
    root <- lapply(seq_len(nrow(g)), function(i) {
      l <- .mf_medlyn_leaf(g$T[i], g$Q[i], g$D[i], g$theta[i], g$g0[i])
      ok <- tryCatch({ l$solve_medlyn_ci_numerical(); TRUE }, error = function(e) FALSE)
      p <- .mf_p(l); p$D <- l$atm_vpd_
      r <- function(x) rep$medlyn_ci_numerical(x, p)
      x <- l$ci_
      if (!ok) return(.probe_pairs(1, as.numeric(r(p$ca) > 0)))
      h <- 1e-4
      A <- rep$assimilation(x, p$vcmax, p$J, p$gamma, p$km, p$theta_c, p$R_d, p$atm_kpa)
      .probe_pairs(c(x, l$assim_colimited_, l$stom_cond_CO2_),
                   c(x - r(x) / ((r(x + h) - r(x - h)) / (2 * h)), A,
                     rep$medlyn_gs(A, p$g0, p$g1, p$D, p$ca, p$atm_kpa, p$theta, p$theta_fc, p$theta_w)))
    })
    dark <- lapply(.mf_T(temps), function(T) {
      l <- .mf_medlyn_leaf(T, PPFD = 0)
      ok <- tryCatch({ l$solve_medlyn_ci_numerical(); TRUE }, error = function(e) FALSE)
      p <- .mf_p(l); p$D <- l$atm_vpd_
      .probe_pairs(as.numeric(!ok), as.numeric(rep$medlyn_ci_numerical(p$ca, p) > 0))
    })
    flat <- lapply(.mf_T(temps), function(T) {
      l <- .mf_medlyn_leaf(T, theta = 0.15, g0 = 0, R_d_25 = 0)
      l$solve_medlyn_ci_numerical()
      p <- .mf_p(l); p$D <- l$atm_vpd_
      lo <- p$gamma * p$atm_kpa * 1e-3
      peak <- max(vapply(seq(lo, p$ca, length.out = 50), function(x) rep$medlyn_ci_numerical(x, p), 0))
      .probe_pairs(c(as.numeric(abs(l$ci_ - lo) < 1e-6 * lo), 1),
                   c(1, as.numeric(peak <= 1e-12)))
    })
    do.call(rbind, c(root, dark, flat))
  },
  medlyn_ci_analytical = function(rep, temps) {
    g <- expand.grid(T = .mf_T(temps), D = c(0.3, 1, 2, 4), g1 = c(1, 2.57, 6))
    out <- lapply(seq_len(nrow(g)), function(i) {
      l <- .mf_medlyn_leaf(g$T[i], D = g$D[i], g1 = g$g1[i])
      l$solve_medlyn_ci_analytical()
      p <- .mf_p(l); p$D <- l$atm_vpd_
      ci <- rep$medlyn_ci_analytical(l$ca_, l$g1, l$atm_vpd_)
      A <- rep$assimilation(l$ci_, p$vcmax, p$J, p$gamma, p$km, p$theta_c, p$R_d, p$atm_kpa)
      # The root of the residual only at g0 = 0, beta = 1.
      w <- p; w$g0 <- 0; w$theta <- w$theta_fc
      .probe_pairs(c(l$ci_, l$assim_colimited_, l$stom_cond_CO2_, 0),
                   c(ci, A, rep$medlyn_gs(A, p$g0, p$g1, p$D, p$ca, p$atm_kpa, p$theta, p$theta_fc, p$theta_w),
                     rep$medlyn_ci_numerical(ci, w)))
    })
    do.call(rbind, out)
  },

  # --- soil moisture ----------------------------------------------------------
  soil_moisture_state = function(rep, temps) {
    l <- leaf_model()
    g <- expand.grid(theta = c(0.1, 0.25, 0.45), theta_fc = c(0.4, 0.6), theta_w = c(0.05, 0.2))
    out <- lapply(seq_len(nrow(g)), function(i) {
      l$theta <- g$theta[i]; l$theta_fc <- g$theta_fc[i]; l$theta_w <- g$theta_w[i]
      l$theta_ <- 99; l$theta_fc_ <- 99; l$theta_w_ <- 99  # undone by the next set_drivers
      .probe_drive(l, 25)
      .probe_pairs(c(l$theta_, l$theta_fc_, l$theta_w_),
                   rep$soil_moisture_state(l$theta, l$theta_fc, l$theta_w))
    })
    do.call(rbind, out)
  },
  # lambda_emergent after an exact CF77 stem solve: gate off, beta inside (0, 1),
  # at 1 and clamped above it; the two refusals as indicators. Divided by the
  # prescribed lambda so the comparison is O(1).
  cf77_soil_price = function(rep, temps) {
    g <- expand.grid(theta = c(0.15, 0.2, 0.21, 0.3, 0.45, 0.5, 0.8), on = c(FALSE, TRUE),
                     lambda = c(1e5, 2e5), theta_fc = c(0.5, 0.2))
    out <- lapply(seq_len(nrow(g)), function(i) {
      l <- leaf_model(supply = leaf_supply_singlelayer())
      l$CF77_lambda_ <- g$lambda[i]
      l$CF77_soil_beta_ <- g$on[i]
      l$theta <- g$theta[i]; l$theta_fc <- g$theta_fc[i]
      set_drivers(l, psi_soil = 0.5, root_network = series_resistance(1500))
      l$set_model("CF77", "stem")
      ok <- tryCatch({ l$optimise(); TRUE }, error = function(e) FALSE)
      mine <- rep$cf77_soil_price(l$CF77_lambda_, l$CF77_soil_beta_, l$theta_, l$theta_fc_, l$theta_w_)
      if (!ok || is.na(mine)) return(.probe_pairs(as.numeric(!ok), as.numeric(is.na(mine))))
      .probe_pairs(l$lambda_emergent / g$lambda[i], mine / g$lambda[i])
    })
    do.call(rbind, out)
  },

  # --- closed form ------------------------------------------------------------
  uso_group = function(rep, temps) {
    g <- expand.grid(T = .mf_T(temps), set = list(list(H2O_CO2_stom_diff_ratio_ = 1.6),
                                                  list(H2O_CO2_stom_diff_ratio_ = 1.67)),
                     psi_soil = c(0.5, 1.5))
    .mf_served(rep, "CF77", g)
  },
  transpiration_from_assim = function(rep, temps) {
    .mf_served(rep, "CF77", expand.grid(D = c(0.5, 1, 2, 3), PPFD = c(300, 1500), psi_soil = c(0.5, 2)))
  },
  stom_cond_from_E = function(rep, temps) {
    g <- expand.grid(psi = c(1.2, 1.5, 2, 3), T = .mf_T(temps), D = c(0.5, 2))
    out <- lapply(seq_len(nrow(g)), function(i) {
      l <- leaf_model(supply = leaf_supply_singlelayer())
      set_drivers(l, psi_soil = 1, PPFD = 900, atm_vpd = g$D[i], leaf_temp = g$T[i],
                  root_network = series_resistance(1e3))
      l$set_leaf_states_rates_from_psi_stem(g$psi[i], 1)
      # x 1e3: the conductance is ~0.1-0.3 mol, so judge it on a scale of order one.
      1e3 * .probe_pairs(l$stom_cond_CO2_,
                         rep$stom_cond_from_E(l$transpiration_, l$atm_kpa_, l$vpd_leaf_,
                                              l$H2O_CO2_stom_diff_ratio_))
    })
    do.call(rbind, out)
  },
  dassim_dci = function(rep, temps) {
    l <- leaf_model()
    g <- expand.grid(T = temps, ci = c(6, 10, 20, 30, 38))
    h <- 1e-4
    do.call(rbind, mapply(function(T, ci) {
      .probe_drive(l, T)
      .probe_pairs((l$assim_colimited(ci + h) - l$assim_colimited(ci - h)) / (2 * h),
                   rep$dassim_dci(ci, l$vcmax_, l$electron_transport_, l$gamma_, l$km_,
                                  l$curv_fact_colim, l$atm_kpa_))
    }, g$T, g$ci, SIMPLIFY = FALSE))
  },
  # Relative to lambda itself, so the scale of the marginal cost (~1e5) drops out.
  dlambda_TF24 = function(rep, temps) {
    h <- 1e-5
    out <- lapply(c(0.8, 1.5, 3), function(b2) {
      l <- leaf_model(supply = leaf_supply_singlelayer())
      set_traits(l, leaf_traits(TF24_beta2 = b2))
      set_drivers(l, psi_soil = 1, root_network = series_resistance(1500))
      psi <- c(0.3, 0.7, 1, 2, 3, 4, 5, 5.5)
      fd <- vapply(psi, function(x) (l$lambda_TF24(x + h) - l$lambda_TF24(x - h)) / (2 * h), 0)
      lam <- vapply(psi, l$lambda_TF24, 0)
      .probe_pairs(fd / lam, rep$dlambda_TF24(psi, l$TF24_cost_scale, l$TF24_beta2, l$stem_c,
                                              l$stem_b, l$leaf_specific_conductance_max_) / lam)
    })
    do.call(rbind, out)
  },
  # Interior rows served; the dark (E < 0) and low-price (E >= E(psi_crit))
  # branches return a bound, which the guard always rejects.
  psi_from_supply = function(rep, temps) {
    served <- .mf_served(rep, "CF77", expand.grid(psi_soil = c(0.5, 1, 2, 3), lambda = c(1e5, 2e5, 4e5)))
    bound <- lapply(list(list(PPFD = 0), list(PPFD = 10, T = 40), list(lambda = 1e4),
                         list(lambda = 3e4, psi_soil = 3)), function(a) {
      r <- do.call(.mf_arm, c(list(rep = rep, curve = "CF77"), a))
      at_bound <- r$s$psi == r$l$single_psi_soil_ || r$s$psi == r$l$psi_crit
      .probe_pairs(as.numeric(r$fb), as.numeric(at_bound))
    })
    rbind(served, do.call(rbind, bound))
  },
  evaluate_at_lambda = function(rep, temps) {
    g <- expand.grid(psi_soil = c(0.3, 1, 2, 3), D = c(0.5, 2, 4), PPFD = c(300, 1500),
                     curve = c("TF24", "TF24_floor"), stringsAsFactors = FALSE)
    out <- lapply(seq_len(nrow(g)), function(i) {
      r <- .mf_arm(rep, g$curve[i], psi_soil = g$psi_soil[i], D = g$D[i], PPFD = g$PPFD[i])
      if (r$fb) .probe_pairs(0, r$guard) else .probe_pairs(c(r$pkg, 1), c(r$s$psi, r$guard))
    })
    do.call(rbind, out)
  },
  wet_end_p = function(rep, temps) {
    .mf_served(rep, "TF24", expand.grid(psi_soil = c(0, 0.3, 1, 2, 3), D = c(0.5, 1, 2),
                                        PPFD = c(500, 1500)))
  },
  closed_form_TF24 = function(rep, temps) {
    g <- expand.grid(psi_soil = c(0.3, 1, 2, 3), D = c(0.5, 1, 2, 4), T = .mf_T(temps))
    exact <- expand.grid(psi_soil = c(0.5, 1, 2), D = c(0.5, 1, 2))
    exact$traits <- list(leaf_traits(TF24_beta2 = 1 / 2.680147))
    rbind(.mf_served(rep, "TF24", g), .mf_served(rep, "TF24", exact))
  },
  closed_form_CF77 = function(rep, temps) {
    g <- expand.grid(psi_soil = c(0.3, 1, 2, 3), lambda = c(1e5, 1.5e5, 3e5), T = .mf_T(temps),
                     set = list(list(), list(CF77_soil_beta_ = TRUE, theta = 0.35),
                                list(CF77_soil_beta_ = TRUE, theta = 0.8)))
    .mf_served(rep, "CF77", g)
  },
  closed_form_TF24_floor = function(rep, temps) {
    g <- expand.grid(psi_soil = c(0.3, 1, 2, 3), lambda_o = c(1e4, 5e4, 1.5e5), D = c(0.5, 2),
                     T = .mf_T(temps))
    .mf_served(rep, "TF24_floor", g)
  },
  # The fall-back flag against the replica guard, on all three arms over a grid
  # that reaches every clause.
  within_guard = function(rep, temps) {
    g <- expand.grid(psi_soil = c(0.3, 0.5, 1, 2, 3, 4.5), D = c(0.5, 1, 2, 4), PPFD = c(10, 300, 1500),
                     T = c(10, 25, 40), curve = c("CF77", "TF24", "TF24_floor"), stringsAsFactors = FALSE)
    out <- lapply(seq_len(nrow(g)), function(i) {
      r <- .mf_arm(rep, g$curve[i], psi_soil = g$psi_soil[i], D = g$D[i], PPFD = g$PPFD[i], T = g$T[i])
      .probe_pairs(as.numeric(!r$fb), r$guard)
    })
    out <- do.call(rbind, out)
    if (length(unique(out$pkg)) < 2) stop("within_guard probe reached one outcome only")
    out
  },
  optimise_closed = function(rep, temps) {
    cfg <- list(list(curve = "TF24"), list(curve = "CF77"), list(curve = "TF24_floor"),
                list(curve = "TF24", collar = TRUE), list(curve = "TF24", layers = 3),
                list(curve = "TF24", eb = TRUE), list(curve = "TF24", damage = TRUE),
                list(curve = "TF24_floor", lambda_o = 0),
                list(curve = "JS22"), list(curve = "CMax"), list(curve = "ProfitMax"),
                list(curve = "SOX"), list(curve = "JW26"), list(curve = "LeastCost"))
    refusals <- lapply(cfg, function(a) {
      layers <- if (is.null(a$layers)) 1 else a$layers
      l <- leaf_model(supply = if (layers > 1) leaf_supply_multilayer() else leaf_supply_singlelayer())
      l$CF77_lambda_ <- 1.5e5
      l$TF24_floor_lambda_o <- if (is.null(a$lambda_o)) 5e4 else a$lambda_o
      if (isTRUE(a$eb)) { l$use_energy_balance_ <- TRUE; l$NIR_int_ <- 1500 / 4.57 }
      if (isTRUE(a$damage)) l$use_thermal_damage_ <- TRUE
      if (layers > 1) {
        set_drivers(l, psi_soil = c(0.5, 0.75, 1.0), soil_depth = c(1, 2, 3),
                    root_network = root_network_from_carbon(rep(1 / 3, 3), c(1, 2, 3)))
      } else {
        set_drivers(l, psi_soil = 1, root_network = series_resistance(1500))
      }
      l$set_model(a$curve, if (isTRUE(a$collar)) "collar" else "stem", "closed")
      refused <- tryCatch({ l$optimise(); FALSE }, error = function(e) TRUE)
      .probe_pairs(as.numeric(refused),
                   rep$optimise_closed(a$curve, isTRUE(a$collar), layers, isTRUE(a$eb),
                                       isTRUE(a$damage), l$TF24_floor_lambda_o))
    })
    # Soil drier than psi_crit: one feasible potential, no fall-back.
    dry <- lapply(c("TF24", "CF77", "TF24_floor"), function(cv) {
      l <- .mf_leaf(cv, psi_soil = 6)
      l$optimise()
      .probe_pairs(c(l$opt_psi_stem_, as.numeric(l$last_solve_fell_back_)), c(l$single_psi_soil_, 0))
    })
    # The counters: three served or fallen-back solves on one leaf.
    l <- .mf_leaf("TF24", psi_soil = 1)
    fbs <- 0
    for (ps in c(1, 4.5, 0.5)) {
      set_drivers(l, psi_soil = ps, PPFD = 1500, atm_vpd = 1, root_network = series_resistance(1500))
      l$optimise()
      fbs <- fbs + l$last_solve_fell_back_
    }
    counters <- .probe_pairs(c(l$closed_form_calls_, l$closed_form_fallbacks_), c(3, fbs))
    rbind(do.call(rbind, refusals), do.call(rbind, dry), counters)
  }
)
