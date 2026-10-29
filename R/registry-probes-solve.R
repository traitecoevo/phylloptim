# registry_check()'s probes for the solve: the ci root-find, the objective, the first-order conditions and the collar bracket.
#
# Same contract as R/registry-probes.R: one probe per equation id in
# inst/model-registry.yml, each driving a Leaf and pairing what the package
# reports with the registry's replica. registry_probes() collects every
# `.registry_probes_*` list in the namespace, so this file needs no other wiring.

# Probes for the solve fragment. Sourced in an environment whose parent is the
# phylloptim namespace, so .probe_pairs, .probe_eb and the exports are in scope.
# Helpers are prefixed .sv_ so they cannot collide with another fragment's.

# A single-potential leaf behind a series resistance, energy balance off.
.sv_single <- function(psi_soil = 1, PPFD = 900, T = 25, r = 1e3, atm_vpd = 2) {
  s <- leaf_model(supply = leaf_supply_singlelayer())
  set_drivers(s, psi_soil = psi_soil, PPFD = PPFD, atm_vpd = atm_vpd, leaf_temp = T,
              root_network = series_resistance(r))
  s
}

# A multi-layer leaf at the default root network.
.sv_multi <- function(psi_soil, PPFD = 900, T = 25, atm_vpd = 2, traits = leaf_traits()) {
  m <- leaf_model(traits = traits)
  set_drivers(m, psi_soil = psi_soil, PPFD = PPFD, atm_vpd = atm_vpd, leaf_temp = T)
  m
}

# A single-potential leaf with the energy balance on (as .probe_eb, any psi_soil).
.sv_eb <- function(psi_soil = 1, PPFD = 900, Tair = 25, atm_vpd = 2) {
  e <- leaf_model(supply = leaf_supply_singlelayer())
  e$use_energy_balance_ <- TRUE
  e$d_ <- 0.05
  e$wind_speed_ <- 2
  e$NIR_int_ <- PPFD / 4.57
  set_drivers(e, psi_soil = psi_soil, PPFD = PPFD, atm_vpd = atm_vpd, leaf_temp = Tair,
              root_network = series_resistance(1e3))
  e
}

# Net assimilation at the leaf's CURRENT temperature block, as a function of ci,
# by the registry's own replica.
.sv_A <- function(rep, l) {
  vcmax <- l$vcmax_; J <- l$electron_transport_; gs <- l$gamma_; km <- l$km_
  th <- l$curv_fact_colim; Rd <- l$R_d_; P <- l$atm_kpa_
  function(ci) rep$assimilation(ci, vcmax, J, gs, km, th, Rd, P)
}

# The supply's critical potential: the root's on the multi-layer path.
.sv_supply_crit <- function(l) if (identical(l$supply_kind, "single")) l$psi_crit else l$root_psi_crit

# The quantities the bracket is built from, read through the bindings: the
# wettest layer, the column flux at psi_crit, and the two continuity roots (NA
# where the root-find has no bracket, which only happens past an earlier exit).
.sv_bracket_inputs <- function(l, ps) {
  wettest <- min(ps)
  l$E_from_Soil_to_Root_Collar(l$psi_crit, ps)
  E_at_crit <- l$E_up_ - l$transpiration(l$psi_crit, l$psi_crit)
  root <- function(crit) tryCatch(l$find_root_psi(wettest, ps, crit), error = function(e) NA_real_)
  list(wettest = wettest, E_at_crit = E_at_crit, root_crit = root(1), root_zero_E = root(0),
       supply_crit = .sv_supply_crit(l))
}

# The floor of the V: the central difference at whichever step lands closest.
.sv_fd_floor <- function(rep, eq, f, x, d, steps = 10^-(3:7)) {
  fd <- vapply(steps, function(h) rep[[eq]](f, x, h), 0)
  fd[which.min(abs(fd - d))]
}

# Supply profiles and traits that reach the collar bracket's regimes.
.sv_regimes <- list(
  list(ps = c(1, 1.5, 2), tr = leaf_traits()),                    # continuity binds
  list(ps = c(3, 3.5, 4), tr = leaf_traits()),
  list(ps = c(4, 4.5), tr = leaf_traits(), D = 0.5),               # boundary-crit
  list(ps = c(3, 3.5, 4), tr = leaf_traits(), T = 40),             # boundary-soil
  list(ps = c(2, 2.5, 3), tr = leaf_traits(root_P50 = 2)),         # root limit binds, interior
  list(ps = c(1, 1.5, 2), tr = leaf_traits(root_P50 = 1.5)),       # root limit binds, pinned
  list(ps = c(3, 3.2), tr = leaf_traits(root_P50 = 2)),            # root limit binds, pinned
  list(ps = c(3, 3.2), tr = leaf_traits(root_P50 = 1.5)),          # inverted: shut-down
  list(ps = c(6, 6.5), tr = leaf_traits()),                        # wettest past psi_crit
  list(ps = c(5.8, 8, 8), tr = leaf_traits(root_P50 = 4)),         # uptake negative at psi_crit
  list(ps = c(5.8, 8, 8), tr = leaf_traits(root_P50 = 2.5)))       # continuity root at psi_crit

.registry_probes_solve <- list(

  # The package's ci against one Newton step of the printed residual from it,
  # gate off over the temperature grid and the light that reaches the
  # compensation point (PPFD 0), and gate on at the leaf temperature the
  # transpiration gives.
  intercellular_co2 = function(rep, temps) {
    out <- list()
    for (Q in c(0, 100, 900, 1800)) {
      l <- .sv_single(PPFD = Q)
      for (T in temps) {
        set_drivers(l, psi_soil = 1, PPFD = Q, leaf_temp = T, root_network = series_resistance(1e3))
        for (psi in c(1.2, 2, 3.5, 5.5)) {
          ci <- l$psi_stem_to_ci(psi, 1)
          gc <- l$stom_cond_CO2(psi, 1)
          out[[length(out) + 1L]] <- .probe_pairs(
            ci, rep$intercellular_co2(ci, .sv_A(rep, l), gc, l$ca_, l$atm_kpa_, l$gamma_))
        }
      }
    }
    for (Ta in c(15, 30, 45, 50)) for (Q in c(0, 300, 1800)) {
      e <- .sv_eb(PPFD = Q, Tair = Ta)
      for (psi in c(1.3, 2.5, 4.5)) {
        e$set_leaf_states_rates_from_psi_stem(psi, 1)
        out[[length(out) + 1L]] <- .probe_pairs(
          e$ci_, rep$intercellular_co2(e$ci_, .sv_A(rep, e), e$stom_cond_CO2_, e$ca_,
                                       e$atm_kpa_, e$gamma_))
      }
    }
    do.call(rbind, out)
  },

  # Read right after set_drivers, when the block is at the driver temperature
  # (with the gate on as well: the air-temperature baseline).
  assim_max = function(rep, temps) {
    out <- list()
    for (Q in c(0, 50, 900, 1800)) for (ca in c(20, 40, 80)) {
      l <- leaf_model()
      for (T in temps) {
        set_drivers(l, psi_soil = 1, PPFD = Q, ca = ca, leaf_temp = T)
        out[[length(out) + 1L]] <- .probe_pairs(
          l$assim_max_, rep$assim_max(l$ca_, l$vcmax_, l$electron_transport_, l$gamma_,
                                      l$km_, l$curv_fact_colim, l$R_d_, l$atm_kpa_))
      }
    }
    for (Ta in c(15, 30, 45)) {
      e <- .sv_eb(Tair = Ta)
      out[[length(out) + 1L]] <- .probe_pairs(
        e$assim_max_, rep$assim_max(e$ca_, e$vcmax_, e$electron_transport_, e$gamma_,
                                    e$km_, e$curv_fact_colim, e$R_d_, e$atm_kpa_))
    }
    do.call(rbind, out)
  },

  # Both branches (0.8 and 1 are no-flow against an upstream potential of 1),
  # gate off over the grid and at PPFD 0, gate on at three air temperatures.
  leaf_state_chain = function(rep, temps) {
    one <- function(l, psi, up) {
      l$set_leaf_states_rates_from_psi_stem(psi, up)
      pkg <- c(l$transpiration_, l$stom_cond_CO2_, l$ci_, l$assim_colimited_)
      A <- .sv_A(rep, l)
      mine <- rep$leaf_state_chain(psi, up, l$transpiration(psi, up), A, l$ca_, l$atm_kpa_,
                                   l$vpd_leaf_, l$gamma_, pkg[3])
      .probe_pairs(pkg, mine)
    }
    out <- list()
    for (Q in c(0, 900)) {
      l <- .sv_single(PPFD = Q)
      for (T in temps[seq(1, length(temps), by = 2)]) {
        set_drivers(l, psi_soil = 1, PPFD = Q, leaf_temp = T, root_network = series_resistance(1e3))
        for (psi in c(0.8, 1, 1.5, 2.5, 4, 5.5)) out[[length(out) + 1L]] <- one(l, psi, 1)
      }
    }
    for (Ta in c(15, 30, 45)) {
      e <- .sv_eb(Tair = Ta)
      for (psi in c(0.9, 1, 1.5, 3)) out[[length(out) + 1L]] <- one(e, psi, 1)
    }
    do.call(rbind, out)
  },

  # Each reachable no-flow exit of the collar solve, gate off and on.
  no_flow_state = function(rep, temps) {
    state <- function(l) c(l$opt_psi_stem_, l$opt_root_psi_, l$transpiration_, l$stom_cond_CO2_,
                           l$ci_, l$assim_colimited_, l$profit_, l$E_up_)
    shut <- function(l, collar) {
      l$optimise()
      pkg <- state(l)
      C <- l$hydraulic_cost_TF(l$psi_crit)
      .probe_pairs(pkg, rep$no_flow_state("hydraulic-shutdown", collar, l$psi_crit, l$R_d_,
                                          l$gamma_, l$atm_kpa_, C))
    }
    shade <- function(l, ps) {
      l$optimise()
      pkg <- state(l)
      x <- pkg[2]; h <- 1e-6
      Eup <- l$E_up_
      l$E_from_Soil_to_Root_Collar(x + h, ps); up <- l$E_up_
      l$E_from_Soil_to_Root_Collar(x - h, ps); dn <- l$E_up_
      C <- l$hydraulic_cost_TF(x)
      mine <- rep$no_flow_state("shade-death", x, l$psi_crit, l$R_d_, l$gamma_, l$atm_kpa_, C,
                                Eup, (up - dn) / (2 * h))
      # E_up is ~1e-20 kg m-2 s-1 here; compared through the collar it implies.
      .probe_pairs(pkg[-8], mine[-8])
    }
    out <- list()
    for (T in c(15, 25, 40)) {
      m <- .sv_multi(c(6, 6.5), T = T)                                      # wettest past psi_crit
      out[[length(out) + 1L]] <- shut(m, m$psi_crit)
      m <- .sv_multi(c(5.8, 8, 8), T = T, traits = leaf_traits(root_P50 = 4))  # E_col(psi_crit) < 0
      out[[length(out) + 1L]] <- shut(m, m$root_psi_crit)
      m <- .sv_multi(c(3, 3.2), T = T, traits = leaf_traits(root_P50 = 1.5))   # inverted bracket
      out[[length(out) + 1L]] <- shut(m, m$root_psi_crit)
      s <- .sv_single(psi_soil = 6, T = T)
      out[[length(out) + 1L]] <- shut(s, s$psi_crit)
      for (ps in list(1, 2.5)) {
        s <- .sv_single(psi_soil = ps, PPFD = 0, T = T)
        out[[length(out) + 1L]] <- shade(s, ps)
      }
      m <- .sv_multi(c(1, 1.5, 2), PPFD = 0, T = T)
      out[[length(out) + 1L]] <- shade(m, c(1, 1.5, 2))
    }
    e <- .sv_eb(psi_soil = 6, Tair = 35)
    out[[length(out) + 1L]] <- shut(e, e$psi_crit)
    e <- .sv_eb(psi_soil = 1, PPFD = 0, Tair = 35)
    out[[length(out) + 1L]] <- shade(e, 1)
    do.call(rbind, out)
  },

  # Prescribed points on both routes: targets inside and outside the interval,
  # four difference curves on the stem route, two on the collar route.
  profit = function(rep, temps) {
    seat <- function(l, curve, route) {
      if (curve == "CF77") l$CF77_lambda_ <- 1.5e5
      l$set_model(curve, route)
    }
    out <- list()
    for (curve in c("TF24", "CF77", "JS22", "CMax")) for (ps in c(1, 3, 6)) {
      s <- .sv_single(psi_soil = ps)
      seat(s, curve, "stem")
      for (target in c(0.5, 1.5, 3, 5, 7)) {
        s$evaluate_psi_stem_at(target)
        out[[length(out) + 1L]] <- .probe_pairs(
          c(s$opt_psi_stem_, s$profit_),
          rep$profit(target, ps, s$psi_crit, s$assim_colimited_, s$hydraulic_cost_))
      }
    }
    sets <- list(list(single = TRUE, ps = 1), list(single = TRUE, ps = 2.5),
                 list(single = FALSE, ps = c(1, 1.5, 2)), list(single = FALSE, ps = c(3, 3.5, 4)))
    for (curve in c("TF24", "CF77")) for (st in sets) {
      l <- if (st$single) .sv_single(psi_soil = st$ps) else .sv_multi(st$ps)
      seat(l, curve, "collar")
      b <- .sv_bracket_inputs(l, st$ps)
      a <- b$root_zero_E; z <- min(b$root_crit, b$supply_crit)
      for (frac in c(-0.5, 0.1, 0.5, 0.9, 1.5)) {
        target <- a + frac * (z - a)
        l$evaluate_root_collar_psi(target)
        out[[length(out) + 1L]] <- .probe_pairs(
          c(l$opt_root_psi_, l$profit_),
          rep$profit(target, a, z, l$assim_colimited_, l$hydraulic_cost_))
      }
    }
    do.call(rbind, out)
  },

  # Gate off. Feasible rows and the compensation-point rows (PPFD 0, where the
  # flag is 0 but the value is a derivative); no-flow sentinels are excluded.
  foc_stem = function(rep, temps) {
    out <- list()
    for (curve in c("TF24", "CF77", "JS22", "CMax")) for (T in c(15, 25, 35)) for (Q in c(0, 300, 1800)) {
      s <- .sv_single(PPFD = Q, T = T)
      if (curve == "CF77") s$CF77_lambda_ <- 1.5e5
      s$set_model(curve, "stem")
      for (psi in c(1.2, 2, 3, 4.5, 5.5)) {
        v <- s$dprofit_dpsi_stem_checked(psi)
        if (v[2] != 1 && Q != 0) next
        out[[length(out) + 1L]] <- .probe_pairs(
          v[1], .sv_fd_floor(rep, "foc_stem", function(x) s$evaluate_psi_stem_at(x), psi, v[1]))
      }
    }
    do.call(rbind, out)
  },

  # Gate off, interior collars, single and multi-layer supply.
  foc_collar = function(rep, temps) {
    sets <- list(list(single = TRUE, ps = 1), list(single = TRUE, ps = 2.5),
                 list(single = FALSE, ps = c(1, 1.5, 2)), list(single = FALSE, ps = c(2, 2.5, 3)))
    out <- list()
    for (curve in c("TF24", "CF77", "JS22")) for (T in c(15, 25, 35)) for (st in sets) {
      l <- if (st$single) .sv_single(psi_soil = st$ps, T = T) else .sv_multi(st$ps, T = T)
      if (curve == "CF77") l$CF77_lambda_ <- 1.5e5
      l$set_model(curve, "collar")
      b <- .sv_bracket_inputs(l, st$ps)
      a <- b$root_zero_E; z <- min(b$root_crit, b$supply_crit)
      for (frac in c(0.05, 0.5, 0.95)) {
        x <- a + frac * (z - a)
        v <- l$dprofit_droot_collar_psi_checked(x)
        if (v[2] != 1) next
        out[[length(out) + 1L]] <- .probe_pairs(
          v[1], .sv_fd_floor(rep, "foc_collar", function(y) l$evaluate_root_collar_psi(y), x, v[1]))
      }
    }
    do.call(rbind, out)
  },

  # Both roots on a single layer behind three resistances and on three profiles.
  continuity_root = function(rep, temps) {
    cases <- list(list(l = .sv_single(psi_soil = 1, r = 1e2), ps = 1),
                  list(l = .sv_single(psi_soil = 1), ps = 1),
                  list(l = .sv_single(psi_soil = 3, r = 1e4), ps = 3),
                  list(l = .sv_multi(c(1, 1.5, 2)), ps = c(1, 1.5, 2)),
                  list(l = .sv_multi(c(0.2, 0.3)), ps = c(0.2, 0.3)),
                  list(l = .sv_multi(c(3, 3.5, 4)), ps = c(3, 3.5, 4)))
    do.call(rbind, lapply(cases, function(cs) {
      l <- cs$l; ps <- cs$ps
      up <- function(x) { l$E_from_Soil_to_Root_Collar(x, ps); l$E_up_ }
      rc <- l$find_root_psi(min(ps), ps, 1)
      r0 <- l$find_root_psi(min(ps), ps, 0)
      .probe_pairs(c(rc, r0),
                   c(rep$continuity_root(rc, function(x) up(x) - l$transpiration(l$psi_crit, x)),
                     rep$continuity_root(r0, up)))
    }))
  },

  # Collars across the bracket, including its wet end, where the flux is ~0.
  stem_from_collar = function(rep, temps) {
    cases <- list(list(l = .sv_single(psi_soil = 1), ps = 1),
                  list(l = .sv_single(psi_soil = 3, r = 1e4), ps = 3),
                  list(l = .sv_multi(c(1, 1.5, 2)), ps = c(1, 1.5, 2)),
                  list(l = .sv_multi(c(3, 3.5, 4)), ps = c(3, 3.5, 4)))
    do.call(rbind, lapply(cases, function(cs) {
      l <- cs$l; ps <- cs$ps
      b <- .sv_bracket_inputs(l, ps)
      a <- b$root_zero_E; z <- min(b$root_crit, b$supply_crit)
      do.call(rbind, lapply(c(0.001, 0.05, 0.3, 0.6, 0.9, 1), function(frac) {
        x <- a + frac * (z - a)
        s <- l$find_psi_stem_from_psi_root(x, ps)
        Eup <- l$E_up_
        # Relative: the fluxes are ~1e-5 kg m-2 s-1, so the residual is read in MPa.
        .probe_pairs(s, rep$stem_from_collar(s, Eup, function(y) l$transpiration(y, x)))
      }))
    }))
  },

  # The collar evaluate_root_collar_psi returns for a far-wet and a far-dry
  # target, in every regime of the bracket, against the ends the replica builds.
  collar_bracket = function(rep, temps) {
    one <- function(l, ps) {
      b <- .sv_bracket_inputs(l, ps)
      l$evaluate_root_collar_psi(0); wet <- l$opt_root_psi_
      l$evaluate_root_collar_psi(100); dry <- l$opt_root_psi_
      .probe_pairs(c(wet, dry),
                   rep$collar_bracket(b$wettest, l$psi_crit, b$E_at_crit, b$root_crit,
                                      b$root_zero_E, b$supply_crit, 1e-3))
    }
    out <- lapply(.sv_regimes, function(rg) one(.sv_multi(rg$ps, traits = rg$tr), rg$ps))
    out[[length(out) + 1L]] <- one(.sv_single(psi_soil = 1), 1)
    out[[length(out) + 1L]] <- one(.sv_single(psi_soil = 1, r = 1e-2), 1)   # collapsed
    out[[length(out) + 1L]] <- one(.sv_single(psi_soil = 6), 6)
    do.call(rbind, out)
  },

  # The location of each solved answer, read as a kind, against the kind the
  # replica derives from the bracket and the end gradients; interior answers also
  # pair dJ/dpsi there with zero.
  operating_point_kind = function(rep, temps) {
    observed <- function(l, a, z, supply_crit) {
      o <- l$opt_root_psi_
      if (l$transpiration_ == 0 && l$opt_psi_stem_ == l$psi_crit) return(6)
      if (l$transpiration_ == 0 && l$opt_psi_stem_ == o) return(7)
      if (is.na(a) || is.na(z)) return(NA_real_)
      if (o == 0.5 * (a + z)) return(5)
      if (o == supply_crit) return(4)
      if (abs(o - a) <= 2e-4 * (1 + a)) return(2)
      if (abs(o - z) <= 2e-4 * (1 + z)) return(3)
      1
    }
    one <- function(l, ps) {
      b <- .sv_bracket_inputs(l, ps)
      a <- b$root_zero_E; z <- min(b$root_crit, b$supply_crit)
      f_lo <- f_hi <- NA_real_
      if (!is.na(a) && !is.na(z) && z > a) {
        f_lo <- l$dprofit_droot_collar_psi_checked(a + 1e-6 * (z - a))[1]
        f_hi <- l$dprofit_droot_collar_psi_checked(z)[1]
      }
      l$optimise()
      k <- observed(l, a, z, b$supply_crit)
      mine <- rep$operating_point_kind(b$wettest, l$psi_crit, b$E_at_crit,
                                       if (is.na(b$root_crit)) Inf else b$root_crit,
                                       b$root_zero_E, b$supply_crit, l$assim_max_,
                                       f_lo, f_hi, 1e-3)
      if (k == 1) {
        # The distance to the stationary point, J'/J'', in MPa: J' alone is
        # J'' (~10) times the root-find's 1e-12.
        o <- l$opt_root_psi_; h <- 1e-6
        g <- l$dprofit_droot_collar_psi_checked(o)[1]
        gp <- (l$dprofit_droot_collar_psi_checked(o + h)[1] -
               l$dprofit_droot_collar_psi_checked(o - h)[1]) / (2 * h)
        return(.probe_pairs(c(k, g / gp), c(mine, 0)))
      }
      .probe_pairs(k, mine)
    }
    out <- list()
    for (rg in .sv_regimes) {
      for (T in c(25, if (!is.null(rg$T)) rg$T else 35)) {
        out[[length(out) + 1L]] <- one(.sv_multi(rg$ps, T = T, traits = rg$tr,
                                                 atm_vpd = if (is.null(rg$D)) 2 else rg$D), rg$ps)
      }
    }
    for (Q in c(300, 1800)) for (D in c(0.5, 2)) for (ps in c(0.5, 1, 3)) {
      out[[length(out) + 1L]] <- one(.sv_single(psi_soil = ps, PPFD = Q, atm_vpd = D), ps)
    }
    for (T in c(25, 40)) for (D in c(0.5, 2)) {
      out[[length(out) + 1L]] <- one(.sv_multi(c(4, 4.5), PPFD = 1800, T = T, atm_vpd = D), c(4, 4.5))
    }
    out[[length(out) + 1L]] <- one(.sv_single(psi_soil = 1, r = 1e-2), 1)   # determined
    out[[length(out) + 1L]] <- one(.sv_single(psi_soil = 1, PPFD = 0), 1)   # shade death
    do.call(rbind, out)
  }
)
