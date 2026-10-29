# registry_check()'s probes for the cost curves and the marginal cost of water.
#
# Same contract as R/registry-probes.R: one probe per equation id in
# inst/model-registry.yml, each driving a Leaf and pairing what the package
# reports with the registry's replica. registry_probes() collects every
# `.registry_probes_*` list in the namespace, so this file needs no other wiring.

# Probes for the `costs` and `marginal_cost` sections of the model registry.
#
# Sourced in an environment whose parent is the phylloptim namespace, so the
# package's own helpers (.probe_pairs, .probe_drive, .probe_eb) are visible.
#
# The derivative checks rest on one decomposition. CMax at a = b = 0 is a curve
# with an identity link and an exactly zero cost, so its dprofit IS the benefit
# term h'A' at h' = 1. Every other curve's dprofit at the same point is then
# h'_K(A) * d0 - C'_K (stem route) or h'_K(A) * d0 - C'_K * s - C'_up,K (collar
# route, s = dpsi_stem/dpsi_collar), with every number on the package side read
# off a binding.

.c_cfg <- list(CF77_lambda_ = 1e5, TF24_floor_lambda_o = 5e4, LeastCost_beta = 146)

# A single-potential leaf seated on `curve`. The collar route needs a series
# resistance; the stem route needs none.
.c_leaf <- function(curve = "TF24", route = "stem", psi_soil = 1, PPFD = 900,
                    atm_vpd = 2, leaf_temp = 25, traits = leaf_traits()) {
  l <- leaf_model(traits = traits, supply = leaf_supply_singlelayer())
  if (route == "collar") {
    set_drivers(l, psi_soil = psi_soil, PPFD = PPFD, atm_vpd = atm_vpd,
                leaf_temp = leaf_temp, root_network = series_resistance(1e3))
  } else {
    set_drivers(l, psi_soil = psi_soil, PPFD = PPFD, atm_vpd = atm_vpd,
                leaf_temp = leaf_temp)
  }
  for (nm in names(.c_cfg)) l[[nm]] <- .c_cfg[[nm]]
  l$set_model(curve, route)
  l
}

# A three-layer leaf on a carbon-built root network: the multilayer case.
.c_multi <- function(curve = "TF24", traits = leaf_traits()) {
  m <- leaf_model(traits = traits)
  set_drivers(m, psi_soil = c(0.5, 0.75, 1.0), soil_depth = c(1, 2, 3), PPFD = 900,
              root_network = root_network_from_carbon(rep(1 / 3, 3), c(1, 2, 3)))
  for (nm in names(.c_cfg)) m[[nm]] <- .c_cfg[[nm]]
  m$set_model(curve, "collar")
  m
}

# Every parameter a cost replica reads, off the leaf.
.c_pars <- function(l) {
  list(gamma = l$TF24_cost_scale, beta = l$TF24_beta2, b = l$stem_b, c = l$stem_c,
       kmax = l$leaf_specific_conductance_max_, psi_crit = l$psi_crit,
       lambda = l$CF77_lambda_, gJ = l$JS22_gamma, a = l$CMax_a, bm = l$CMax_b,
       lo = l$TF24_floor_lambda_o, beta_lc = l$LeastCost_beta)
}

# Relative pairs: both sides divided by the replica's magnitude, so a quantity
# far below 1 is judged relatively rather than absolutely.
.c_rel <- function(pkg, mine) {
  s <- ifelse(mine == 0, 1, abs(mine))
  .probe_pairs(pkg / s, mine / s)
}

.c_curves <- c("TF24", "CF77", "JS22", "CMax", "SOX", "JW26", "ProfitMax",
               "TF24_floor", "LeastCost")

# Trait sets that move every cost parameter: a concave (beta < 1) and a convex
# TF24 cost, a signed CMax_b, a narrower vulnerability curve.
.c_trait_sets <- list(
  leaf_traits(),
  leaf_traits(TF24_beta2 = 0.8, TF24_cost_scale = 3, JS22_gamma = 0.3, CMax_a = 0.9,
              CMax_b = -0.4),
  leaf_traits(TF24_beta2 = 2.5, stem_P50 = 2.2, stem_c = 4, JS22_gamma = 1.7,
              CMax_b = 0.3))

# Grid of stem potentials from just above psi_soil to just below psi_crit.
.c_psi <- function(l, psi_soil, n = 7) {
  psi_soil + (l$psi_crit - psi_soil) * c(1e-3, seq(0.05, 0.95, length.out = n - 2), 0.999)
}

# The benefit term d0 = A' at h' = 1, on either route, from CMax at a = b = 0.
.c_zero_cost <- function(route, psi_soil = 1, traits = leaf_traits(), ...) {
  traits$CMax_a <- 0; traits$CMax_b <- 0
  .c_leaf("CMax", route, psi_soil = psi_soil, traits = traits, ...)
}

# Run `f(l, p, ts)` over the trait sets on a stem-route leaf of `curve`.
.c_over_traits <- function(curve, f, psi_soil = 1, route = "stem") {
  do.call(rbind, lapply(.c_trait_sets, function(ts) {
    l <- .c_leaf(curve, route, psi_soil = psi_soil, traits = ts)
    f(l, .c_pars(l), ts)
  }))
}

.registry_probes_costs <- list(

  # --- the cost curves --------------------------------------------------------
  cost_TF24 = function(rep, temps) {
    .c_over_traits("TF24", function(l, p, ts) {
      psi <- c(1e-3, 0.3, .c_psi(l, 1), 1.2 * l$psi_crit)
      pkg <- vapply(psi, l$hydraulic_cost_TF, 0)
      field <- vapply(psi, function(x) { l$hydraulic_cost_TF(x); l$hydraulic_cost_ }, 0)
      .c_rel(c(pkg, field), rep(rep$cost_TF24(psi, p$gamma, p$beta, p$b, p$c), 2))
    })
  },
  # E is the package's (the integral is the hydraulics section's); the price
  # with the soil-beta gate off, on below field capacity, and on above it (clamp).
  cost_CF77 = function(rep, temps) {
    out <- list()
    for (gate in list(NULL, c(0.2, 0.1, 0.3), c(0.35, 0.1, 0.3), c(0.12, 0.05, 0.4))) {
      l <- .c_leaf("CF77")
      if (!is.null(gate)) {
        l$CF77_soil_beta_ <- TRUE
        l$theta <- gate[1]; l$theta_w <- gate[2]; l$theta_fc <- gate[3]
        set_drivers(l, psi_soil = 1, PPFD = 900, atm_vpd = 2, leaf_temp = 25)
      }
      psi <- c(1, .c_psi(l, 1))
      E <- vapply(psi, function(x) l$transpiration(x, 1), 0)
      pkg <- vapply(psi, function(x) l$hydraulic_cost_CF77(x, 1), 0)
      mine <- rep$cost_CF77(E, l$CF77_lambda_, isTRUE(l$CF77_soil_beta_),
                            l$theta_, l$theta_w_, l$theta_fc_)
      out[[length(out) + 1L]] <- .c_rel(c(pkg, l$hydraulic_cost_), c(mine, mine[length(mine)]))
    }
    do.call(rbind, out)
  },
  cost_JS22 = function(rep, temps) {
    .c_over_traits("JS22", function(l, p, ts) {
      g <- expand.grid(psi = c(1, 1 + 1e-8, 1.5, 3, 5), up = c(0.5, 1))
      .c_rel(mapply(l$hydraulic_cost_JS22, g$psi, g$up), rep$cost_JS22(g$psi, g$up, p$gJ))
    })
  },
  # The factored form, including a 1e-10 drop where the difference of squares
  # would cancel.
  cost_CMax = function(rep, temps) {
    .c_over_traits("CMax", function(l, p, ts) {
      g <- expand.grid(psi = c(1, 1 + 1e-10, 1.5, 3, 5), up = c(0.5, 1))
      .c_rel(mapply(l$hydraulic_cost_CMax, g$psi, g$up), rep$cost_CMax(g$psi, g$up, p$a, p$bm))
    })
  },
  # At the probe's lambda_o, at lambda_o = 0 (TF24 exactly) and with the TF24
  # half zeroed (CF77 exactly).
  cost_TF24_floor = function(rep, temps) {
    out <- list()
    for (ts in .c_trait_sets) {
      for (lo in c(5e4, 0)) {
        l <- .c_leaf("TF24_floor", traits = ts); l$TF24_floor_lambda_o <- lo
        p <- .c_pars(l)
        psi <- c(1, .c_psi(l, 1))
        E <- vapply(psi, function(x) l$transpiration(x, 1), 0)
        out[[length(out) + 1L]] <- .c_rel(
          vapply(psi, function(x) l$hydraulic_cost_TF24_floor(x, 1), 0),
          rep$cost_TF24_floor(psi, E, p$gamma, p$beta, p$b, p$c, lo))
      }
    }
    z <- .c_leaf("TF24_floor", traits = leaf_traits(TF24_cost_scale = 0))
    psi <- .c_psi(z, 1)
    E <- vapply(psi, function(x) z$transpiration(x, 1), 0)
    rbind(do.call(rbind, out),
          .c_rel(vapply(psi, function(x) z$hydraulic_cost_TF24_floor(x, 1), 0),
                 rep$cost_TF24_floor(psi, E, 0, z$TF24_beta2, z$stem_b, z$stem_c, 5e4)))
  },
  # After a solve on each curve: lambda_o * E on TF24_floor, zero elsewhere.
  shadow_cost = function(rep, temps) {
    do.call(rbind, lapply(c("TF24_floor", "TF24", "CF77"), function(cv) {
      do.call(rbind, lapply(c("stem", "collar"), function(route) {
        do.call(rbind, lapply(c(0.5, 1, 2), function(ps) {
          l <- .c_leaf(cv, route, psi_soil = ps)
          l$optimise()
          .c_rel(l$shadow_cost, rep$shadow_cost(l$transpiration_, l$TF24_floor_lambda_o, cv))
        }))
      }))
    }))
  },
  # Past psi_crit the clamp binds.
  sox_reduction = function(rep, temps) {
    .c_over_traits("SOX", function(l, p, ts) {
      psi <- c(0, 1e-3, .c_psi(l, 0.5), l$psi_crit, 1.1 * l$psi_crit, 2 * l$psi_crit)
      .probe_pairs(vapply(psi, l$sox_reduction, 0), rep$sox_reduction(psi, p$b, p$c))
    })
  },
  # No binding: read through lambda_SOX = -A (g'/g) / (kmax f), inverted with
  # the package's own g, f and A at the same point.
  sox_reduction_deriv = function(rep, temps) {
    .c_over_traits("SOX", function(l, p, ts) {
      psi <- .c_psi(l, 1)
      pkg <- vapply(psi, function(x) {
        l$set_leaf_states_rates_from_psi_stem(x, 1)
        -l$lambda_SOX(x, 1) * l$sox_reduction(x) * p$kmax * l$proportion_of_conductivity(x) /
          l$assim_colimited_
      }, 0)
      .c_rel(pkg, rep$sox_reduction_deriv(psi, p$b, p$c))
    })
  },
  jw26_reduction = function(rep, temps) {
    .c_over_traits("JW26", function(l, p, ts) {
      psi <- c(0, 1e-3, .c_psi(l, 0.5), l$psi_crit, 1.1 * l$psi_crit, 2 * l$psi_crit)
      .probe_pairs(vapply(psi, l$jw26_reduction, 0), rep$jw26_reduction(psi, p$psi_crit))
    })
  },
  # Unreachable from R as itself (only tests/cpp calls it). Checked as the
  # -g'/g it collapses to in cost_deriv<JW26>, read off the stem derivative:
  # C' = d0/A - dprofit_JW26.
  jw26_reduction_deriv = function(rep, temps) {
    do.call(rbind, lapply(.c_trait_sets, function(ts) {
      l <- .c_leaf("JW26", traits = ts); z <- .c_zero_cost("stem", traits = ts)
      psi <- .c_psi(l, 1)
      pkg <- vapply(psi, function(x) {
        d <- l$dprofit_dpsi_stem_checked(x); d0 <- z$dprofit_dpsi_stem_checked(x)
        l$set_leaf_states_rates_from_psi_stem(x, 1)
        if (d[2] != 1 || d0[2] != 1) return(NA_real_)
        d0[1] / l$assim_colimited_ - d[1]
      }, 0)
      keep <- is.finite(pkg)
      .c_rel(pkg[keep], -rep$jw26_reduction_deriv(psi[keep], l$psi_crit) /
                         rep$jw26_reduction(psi[keep], l$psi_crit))
    }))
  },
  # W read back from the LeastCost objective, A/W.
  leastcost_flux = function(rep, temps) {
    out <- list()
    for (beta in c(20, 146, 1000)) {
      for (Q in c(300, 1500)) {
        l <- .c_leaf("LeastCost", PPFD = Q); l$LeastCost_beta <- beta
        for (x in .c_psi(l, 1)) {
          prof <- l$evaluate_psi_stem_at(x)
          out[[length(out) + 1L]] <- .c_rel(l$assim_colimited_ / prof,
            rep$leastcost_flux(l$transpiration_, l$vcmax_, beta))
        }
      }
    }
    do.call(rbind, out)
  },
  # |A|max against the largest |A| the package reports on a grid of the
  # stream, refined three times about the grid argmax. The 48 C row has A < 0
  # throughout, so |A|max is max(-A) there.
  profitmax_A_max = function(rep, temps) {
    # (psi_soil, PPFD, T, energy balance, deficit). With the energy balance off
    # A rises to the dry bound; with it on, in cold air, cooling by transpiration
    # puts the maximum inside the stream and the root-find branch is reached.
    cases <- list(c(1, 900, 25, 0, 2), c(2.5, 900, 25, 0, 2), c(1, 150, 15, 0, 2),
                  c(0.5, 1800, 35, 0, 2), c(1, 900, 48, 0, 2), c(1, 300, 5, 1, 1),
                  c(1, 1800, 5, 1, 3), c(1, 300, 15, 1, 1), c(1, 900, 25, 1, 2))
    do.call(rbind, lapply(cases, function(cs) {
      l <- if (cs[4] == 1) {
        e <- .probe_eb(Tair = cs[3], PPFD = cs[2], atm_vpd = cs[5])
        e$set_model("ProfitMax", "stem"); e
      } else .c_leaf("ProfitMax", psi_soil = cs[1], PPFD = cs[2], leaf_temp = cs[3])
      l$prepare_profitmax()
      pkg <- l$profitmax_A_max
      A_of <- function(x) vapply(x, function(y) {
        l$set_leaf_states_rates_from_psi_stem(y, cs[1]); l$assim_colimited_ }, 0)
      lo <- cs[1] + 1e-9 * (l$psi_crit - cs[1])
      x <- seq(lo, l$psi_crit, length.out = 401)
      A <- A_of(x)
      for (k in 1:3) {
        i <- if (max(A) > 0) which.max(A) else which.max(-A)
        h <- diff(x[1:2])
        x <- seq(max(lo, x[i] - h), min(l$psi_crit, x[i] + h), length.out = 201)
        A <- c(A[i], A_of(x))
        x <- c(x[1], x)
      }
      .c_rel(pkg, rep$profitmax_A_max(A))
    }))
  },
  profitmax_normalisers = function(rep, temps) {
    do.call(rbind, lapply(.c_trait_sets, function(ts) {
      do.call(rbind, lapply(c(0.3, 1, 2.5), function(ps) {
        l <- .c_leaf("ProfitMax", psi_soil = ps, traits = ts)
        l$prepare_profitmax()
        .c_rel(l$profitmax_k_span,
               rep$profitmax_normalisers(ps, l$psi_crit, l$leaf_specific_conductance_max_,
                                         l$stem_b, l$stem_c)[["k_span"]])
      }))
    }))
  },
  # The three components and the objective, gates off, thermal cost on with the
  # energy balance off (TC constant at the driver temperature) and both on (TC
  # at the leaf temperature the candidate's transpiration gives).
  profitmax_objective = function(rep, temps) {
    out <- list()
    for (cfg in list(c(FALSE, FALSE, 25), c(TRUE, FALSE, 40), c(TRUE, TRUE, 30),
                     c(TRUE, TRUE, 42), c(FALSE, TRUE, 30))) {
      tc <- as.logical(cfg[1]); eb <- as.logical(cfg[2]); Ta <- cfg[3]
      l <- if (eb) .probe_eb(Tair = Ta) else .c_leaf("ProfitMax", leaf_temp = Ta)
      l$use_thermal_cost_ <- tc
      if (eb) set_drivers(l, psi_soil = 1, PPFD = 900, atm_vpd = 2, leaf_temp = Ta,
                          root_network = series_resistance(1e3))
      l$prepare_profitmax()
      for (x in .c_psi(l, 1)) {
        prof <- l$profit_psi_stem_ProfitMax(x, 1)
        T <- if (eb) l$Tleaf_ else Ta
        mine <- rep$profitmax_objective(l$assim_colimited_, l$profitmax_A_max, x, 1, l$psi_crit,
                                        l$leaf_specific_conductance_max_, l$stem_b, l$stem_c,
                                        T, l$T50_, l$Tcrit_, tc)
        out[[length(out) + 1L]] <- .probe_pairs(
          c(l$carbon_gain_, l$hydraulic_cost_norm_, l$thermal_cost_, prof), mine)
      }
    }
    do.call(rbind, out)
  },
  # The seated curve's objective at a prescribed stem potential, every curve:
  # the costs are the replicas above, A and E the package's at that point.
  objective = function(rep, temps) {
    out <- list()
    for (cv in .c_curves) {
      for (ts in .c_trait_sets[1:2]) {
        l <- .c_leaf(cv, traits = ts); p <- .c_pars(l)
        for (x in .c_psi(l, 1)) {
          prof <- l$evaluate_psi_stem_at(x)
          A <- l$assim_colimited_; E <- l$transpiration_
          C <- switch(cv,
            TF24 = rep$cost_TF24(x, p$gamma, p$beta, p$b, p$c),
            CF77 = rep$cost_CF77(E, p$lambda),
            JS22 = rep$cost_JS22(x, 1, p$gJ),
            CMax = rep$cost_CMax(x, 1, p$a, p$bm),
            TF24_floor = rep$cost_TF24_floor(x, E, p$gamma, p$beta, p$b, p$c, p$lo),
            NA_real_)
          g <- switch(cv, SOX = rep$sox_reduction(x, p$b, p$c),
                      JW26 = rep$jw26_reduction(x, p$psi_crit), NA_real_)
          W <- if (cv == "LeastCost") rep$leastcost_flux(E, l$vcmax_, p$beta_lc) else NA_real_
          pm <- if (cv == "ProfitMax") {
            rep$profitmax_objective(A, l$profitmax_A_max, x, 1, p$psi_crit, p$kmax, p$b, p$c,
                                    25, l$T50_, l$Tcrit_, FALSE)[4]
          } else NA_real_
          out[[length(out) + 1L]] <- .c_rel(prof, rep$objective(cv, A, C, g, W, pm))
        }
      }
    }
    do.call(rbind, out)
  },
  # h'(A) as the package applies it: (dprofit_K + C'_K) / d0 on the stem route.
  benefit_link_deriv = function(rep, temps) {
    z <- .c_zero_cost("stem")
    out <- list()
    for (cv in .c_curves) {
      l <- .c_leaf(cv)
      if (cv == "ProfitMax") l$prepare_profitmax()
      p <- .c_pars(l); p$k_span <- l$profitmax_k_span
      for (x in .c_psi(l, 1)) {
        d <- l$dprofit_dpsi_stem_checked(x); d0 <- z$dprofit_dpsi_stem_checked(x)
        if (d[2] != 1 || d0[2] != 1) next
        l$set_leaf_states_rates_from_psi_stem(x, 1)
        p$E <- l$transpiration_; p$vcmax <- l$vcmax_
        out[[length(out) + 1L]] <- .c_rel(
          (d[1] + rep$cost_deriv(cv, x, 1, p)) / d0[1],
          rep$benefit_link_deriv(cv, l$assim_colimited_, l$profitmax_A_max))
      }
    }
    do.call(rbind, out)
  },
  # C'_K = h'_K(A) d0 - dprofit_K, every curve, every trait set, two soils.
  cost_deriv = function(rep, temps) {
    out <- list()
    for (ts in .c_trait_sets) {
      for (ps in c(0.5, 1.5)) {
        z <- .c_zero_cost("stem", psi_soil = ps, traits = ts)
        for (cv in .c_curves) {
          l <- .c_leaf(cv, psi_soil = ps, traits = ts)
          if (cv == "ProfitMax") l$prepare_profitmax()
          p <- .c_pars(l); p$k_span <- l$profitmax_k_span
          for (x in .c_psi(l, ps)) {
            d <- l$dprofit_dpsi_stem_checked(x); d0 <- z$dprofit_dpsi_stem_checked(x)
            if (d[2] != 1 || d0[2] != 1) next
            l$set_leaf_states_rates_from_psi_stem(x, ps)
            p$E <- l$transpiration_; p$vcmax <- l$vcmax_
            hp <- rep$benefit_link_deriv(cv, l$assim_colimited_, l$profitmax_A_max)
            out[[length(out) + 1L]] <- .c_rel(hp * d0[1] - d[1], rep$cost_deriv(cv, x, ps, p))
          }
        }
      }
    }
    do.call(rbind, out)
  },
  # The collar route, single-layer and three-layer: s from TF24 (no upstream
  # term), then C'_up,K = h'_K d0 - C'_K s - dprofit_K for every curve.
  cost_deriv_upstream = function(rep, temps) {
    out <- list()
    setups <- list(
      list(mk = function(cv) .c_leaf(cv, "collar"),
           z = .c_zero_cost("collar"), r = c(1.05, 1.2, 1.5, 1.9)),
      list(mk = function(cv) .c_multi(cv),
           z = .c_multi("CMax", leaf_traits(CMax_a = 0, CMax_b = 0)),
           r = c(1.0, 1.1, 1.2, 1.35)))
    for (su in setups) {
      tf <- su$mk("TF24")
      for (cv in .c_curves) {
        l <- su$mk(cv)
        if (cv == "ProfitMax") l$prepare_profitmax()
        for (r in su$r) {
          d <- l$dprofit_droot_collar_psi_checked(r)
          d0 <- su$z$dprofit_droot_collar_psi_checked(r)
          dT <- tf$dprofit_droot_collar_psi_checked(r)
          if (d[2] != 1 || d0[2] != 1 || dT[2] != 1) next
          l$evaluate_root_collar_psi(r)
          ps <- l$opt_psi_stem_
          p <- .c_pars(l); p$k_span <- l$profitmax_k_span
          p$E <- l$transpiration_; p$vcmax <- l$vcmax_
          s <- (d0[1] - dT[1]) / rep$cost_deriv("TF24", ps, r, p)
          hp <- rep$benefit_link_deriv(cv, l$assim_colimited_, l$profitmax_A_max)
          pkg <- hp * d0[1] - rep$cost_deriv(cv, ps, r, p) * s - d[1]
          out[[length(out) + 1L]] <- .probe_pairs(pkg, rep$cost_deriv_upstream(cv, ps, r, p))
        }
      }
    }
    do.call(rbind, out)
  },

  # --- the marginal cost of water ---------------------------------------------
  lambda_TF24 = function(rep, temps) {
    .c_over_traits("TF24", function(l, p, ts) {
      psi <- c(0.05, 0.3, .c_psi(l, 1))
      .probe_pairs(vapply(psi, l$lambda_TF24, 0),
                   rep$lambda_TF24(psi, p$gamma, p$beta, p$b, p$c, p$kmax))
    })
  },
  lambda_JS22 = function(rep, temps) {
    .c_over_traits("JS22", function(l, p, ts) {
      g <- expand.grid(psi = .c_psi(l, 1), up = c(0.5, 1))
      .probe_pairs(mapply(l$lambda_JS22, g$psi, g$up),
                   rep$lambda_JS22(g$psi, g$up, p$gJ, p$b, p$c, p$kmax))
    })
  },
  lambda_CMax = function(rep, temps) {
    .c_over_traits("CMax", function(l, p, ts) {
      psi <- .c_psi(l, 1)
      .probe_pairs(vapply(psi, function(x) l$lambda_CMax(x, 1), 0),
                   rep$lambda_CMax(psi, p$a, p$bm, p$b, p$c, p$kmax))
    })
  },
  lambda_TF24_floor = function(rep, temps) {
    .c_over_traits("TF24_floor", function(l, p, ts) {
      do.call(rbind, lapply(c(0, 5e4, 2e5), function(lo) {
        l$TF24_floor_lambda_o <- lo
        psi <- c(0.05, .c_psi(l, 1))
        .probe_pairs(vapply(psi, l$lambda_TF24_floor, 0),
                     rep$lambda_TF24_floor(psi, p$gamma, p$beta, p$b, p$c, p$kmax, lo))
      }))
    })
  },
  # A is the package's at the point, which lambda_SOX reads off the leaf.
  lambda_SOX = function(rep, temps) {
    .c_over_traits("SOX", function(l, p, ts) {
      psi <- .c_psi(l, 1)
      pkg <- mine <- numeric(length(psi))
      for (i in seq_along(psi)) {
        l$set_leaf_states_rates_from_psi_stem(psi[i], 1)
        pkg[i] <- l$lambda_SOX(psi[i], 1)
        mine[i] <- rep$lambda_SOX(psi[i], l$assim_colimited_, p$b, p$c, p$kmax)
      }
      .probe_pairs(pkg, mine)
    })
  },
  lambda_JW26 = function(rep, temps) {
    .c_over_traits("JW26", function(l, p, ts) {
      psi <- .c_psi(l, 1)
      pkg <- mine <- numeric(length(psi))
      for (i in seq_along(psi)) {
        l$set_leaf_states_rates_from_psi_stem(psi[i], 1)
        pkg[i] <- l$lambda_JW26(psi[i], 1)
        mine[i] <- rep$lambda_JW26(psi[i], l$assim_colimited_, p$psi_crit, p$b, p$c, p$kmax)
      }
      .probe_pairs(pkg, mine)
    })
  },
  # No binding: read through lambda_emergent after a stem-route ProfitMax solve,
  # which is lambda_ProfitMax at the operating point. Gates off, and the energy
  # balance with the thermal cost on, where the thermal term is non-zero.
  lambda_ProfitMax = function(rep, temps) {
    out <- list()
    for (ps in c(0.5, 1, 2)) {
      for (Q in c(300, 1500)) {
        l <- .c_leaf("ProfitMax", psi_soil = ps, PPFD = Q)
        l$optimise()
        out[[length(out) + 1L]] <- .probe_pairs(l$lambda_emergent,
          rep$lambda_ProfitMax(l$opt_psi_stem_, l$profitmax_A_max, l$profitmax_k_span,
                               l$stem_b, l$stem_c))
      }
    }
    for (Ta in c(25, 35, 42)) {
      e <- .probe_eb(Tair = Ta)
      e$use_thermal_cost_ <- TRUE
      set_drivers(e, psi_soil = 1, PPFD = 900, atm_vpd = 2, leaf_temp = Ta,
                  root_network = series_resistance(1e3))
      e$set_model("ProfitMax", "stem")
      e$optimise()
      Tl <- e$leaf_temp_from_E(e$transpiration_)
      out[[length(out) + 1L]] <- .probe_pairs(e$lambda_emergent,
        rep$lambda_ProfitMax(e$opt_psi_stem_, e$profitmax_A_max, e$profitmax_k_span,
                             e$stem_b, e$stem_c, Tl, e$T50_, e$Tcrit_, e$ra_, e$g_rad_, TRUE))
    }
    do.call(rbind, out)
  },
  # No binding: read through lambda_emergent after a stem-route LeastCost solve.
  lambda_LeastCost = function(rep, temps) {
    out <- list()
    for (beta in c(50, 146, 600)) {
      for (Q in c(300, 1500)) {
        for (D in c(1, 3)) {
          l <- .c_leaf("LeastCost", PPFD = Q, atm_vpd = D); l$LeastCost_beta <- beta
          l$optimise()
          out[[length(out) + 1L]] <- .probe_pairs(l$lambda_emergent,
            rep$lambda_LeastCost(l$assim_colimited_, l$transpiration_, l$vcmax_, beta))
        }
      }
    }
    do.call(rbind, out)
  },
  # Each bound lambda against (dC/dpsi)/(dE/dpsi) of the package's own cost and
  # transpiration by central difference, and the product curves' against
  # A (d(-log g)/dpsi)/(dE/dpsi). Step 1e-5 MPa from a sweep.
  lambda_ratio = function(rep, temps) {
    h <- 1e-5
    out <- list()
    for (ts in .c_trait_sets) {
      l <- .c_leaf("TF24", traits = ts)
      Ef <- function(x) l$transpiration(x, 1)
      for (x in .c_psi(l, 1)[2:6]) {
        l$set_leaf_states_rates_from_psi_stem(x, 1); A <- l$assim_colimited_
        out[[length(out) + 1L]] <- .probe_pairs(
          c(l$lambda_TF24(x), l$lambda_JS22(x, 1), l$lambda_CMax(x, 1),
            l$lambda_TF24_floor(x), l$lambda_SOX(x, 1), l$lambda_JW26(x, 1)),
          c(rep$lambda_ratio(l$hydraulic_cost_TF, Ef, x, h),
            rep$lambda_ratio(function(y) l$hydraulic_cost_JS22(y, 1), Ef, x, h),
            rep$lambda_ratio(function(y) l$hydraulic_cost_CMax(y, 1), Ef, x, h),
            rep$lambda_ratio(function(y) l$hydraulic_cost_TF24_floor(y, 1), Ef, x, h),
            rep$lambda_ratio(function(y) -log(l$sox_reduction(y)), Ef, x, h, 1 / A),
            rep$lambda_ratio(function(y) -log(l$jw26_reduction(y)), Ef, x, h, 1 / A)))
      }
    }
    do.call(rbind, out)
  },
  # TF24's lambda at whatever optimum the last solve found, on every curve and
  # both routes.
  marginal_cost_water = function(rep, temps) {
    do.call(rbind, lapply(.c_curves, function(cv) {
      do.call(rbind, lapply(c("stem", "collar"), function(route) {
        l <- .c_leaf(cv, route); l$optimise()
        .probe_pairs(l$lambda, rep$marginal_cost_water(l$opt_psi_stem_, l$TF24_cost_scale,
                                                       l$TF24_beta2, l$stem_b, l$stem_c,
                                                       l$leaf_specific_conductance_max_))
      }))
    }))
  },
  marginal_cost_water_molar = function(rep, temps) {
    do.call(rbind, lapply(c("TF24", "CF77", "SOX"), function(cv) {
      do.call(rbind, lapply(c(0.5, 1, 2), function(ps) {
        l <- .c_leaf(cv, psi_soil = ps); l$optimise()
        .c_rel(l$lambda_molar, rep$marginal_cost_water_molar(l$lambda))
      }))
    }))
  },
  # S = dE_up/dpsi_collar by central difference of E_from_Soil_to_Root_Collar,
  # single-layer and three-layer, every curve; lambda_multilayer and the
  # collar-route lambda_emergent (which is the same number) both checked.
  marginal_cost_water_multilayer = function(rep, temps) {
    h <- 1e-4
    out <- list()
    for (cv in .c_curves) {
      for (multi in c(FALSE, TRUE)) {
        l <- if (multi) .c_multi(cv) else .c_leaf(cv, "collar")
        l$optimise()
        soil <- if (multi) l$psi_soil_ else 1
        Eup <- function(x) { l$E_from_Soil_to_Root_Collar(x, soil); l$E_up_ }
        r <- l$opt_root_psi_; ps <- l$opt_psi_stem_
        lam <- rep$lambda_TF24(ps, l$TF24_cost_scale, l$TF24_beta2, l$stem_b, l$stem_c,
                               l$leaf_specific_conductance_max_)
        k_r <- l$leaf_specific_conductance_max_ * l$proportion_of_conductivity(r)
        # read at the optimum, before the difference moves E_up_
        pkg <- l$lambda_multilayer
        mine <- rep$marginal_cost_water_multilayer(lam, k_r, Eup, r, h)
        out[[length(out) + 1L]] <- .probe_pairs(pkg, mine)
      }
    }
    do.call(rbind, out)
  },
  # Collar route, every curve, at its interior optima: lambda_emergent against
  # dA/dE along the collar path. An interior optimum is one where the bound
  # first-order condition is feasible and vanishes; a pinned one carries no
  # such identity and is skipped. CF77 at 3e5, where it is interior; at the
  # shared 1e5 it pins at psi_crit on these supplies.
  lambda_collar = function(rep, temps) {
    h <- 1e-4
    out <- list()
    for (cv in .c_curves) {
      for (setup in list(list(multi = FALSE, ps = 0.5), list(multi = FALSE, ps = 1),
                         list(multi = TRUE))) {
        l <- if (setup$multi) .c_multi(cv) else
          .c_leaf(cv, "collar", psi_soil = setup$ps, PPFD = 1500)
        if (cv == "CF77") l$CF77_lambda_ <- 3e5
        l$optimise()
        r <- l$opt_root_psi_
        g <- l$dprofit_droot_collar_psi_checked(r)
        if (!(g[2] == 1) || !is.finite(g[1]) || abs(g[1]) > 1e-8) next
        lam <- l$lambda_emergent
        Af <- function(x) { l$evaluate_root_collar_psi(x); l$assim_colimited_ }
        Ef <- function(x) { l$evaluate_root_collar_psi(x); l$transpiration_ }
        out[[length(out) + 1L]] <- .c_rel(lam, rep$lambda_collar(Af, Ef, r, h))
      }
    }
    do.call(rbind, out)
  },
  # Stem route, every curve: the seated curve's lambda at its own optimum.
  lambda_emergent = function(rep, temps) {
    out <- list()
    for (cv in .c_curves) {
      for (ps in c(0.5, 1.5)) {
        for (Q in c(400, 1500)) {
          l <- .c_leaf(cv, psi_soil = ps, PPFD = Q); l$optimise()
          p <- .c_pars(l)
          p$A <- l$assim_colimited_; p$E <- l$transpiration_; p$vcmax <- l$vcmax_
          p$A_max <- l$profitmax_A_max; p$k_span <- l$profitmax_k_span
          out[[length(out) + 1L]] <- .probe_pairs(l$lambda_emergent,
            rep$lambda_emergent(cv, l$opt_psi_stem_, ps, p))
        }
      }
    }
    do.call(rbind, out)
  },
  # After solves on both routes and with the energy balance on, where the leaf
  # deficit differs from the air's.
  g1_eff = function(rep, temps) {
    out <- list()
    for (cv in c("TF24", "CF77", "SOX")) {
      for (route in c("stem", "collar")) {
        for (D in c(0.8, 2, 4)) {
          l <- .c_leaf(cv, route, atm_vpd = D); l$optimise()
          out[[length(out) + 1L]] <- .probe_pairs(l$g1_eff, rep$g1_eff(l$ci_, l$ca_, l$vpd_leaf_))
        }
      }
    }
    for (Ta in c(20, 35)) {
      e <- .probe_eb(Tair = Ta); e$optimise()
      out[[length(out) + 1L]] <- .probe_pairs(e$g1_eff, rep$g1_eff(e$ci_, e$ca_, e$vpd_leaf_))
    }
    do.call(rbind, out)
  },
  stom_cond_H2O = function(rep, temps) {
    out <- list()
    for (ratio in c(1.67, 1.6)) {
      for (ps in c(0.5, 1, 2)) {
        l <- .c_leaf("TF24", psi_soil = ps); l$H2O_CO2_stom_diff_ratio_ <- ratio; l$optimise()
        out[[length(out) + 1L]] <- .c_rel(l$gs_H2O, rep$stom_cond_H2O(l$stom_cond_CO2_, ratio))
      }
    }
    do.call(rbind, out)
  }
)
