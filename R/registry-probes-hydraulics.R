# registry_check()'s probes for the vulnerability curves and the water supply.
#
# Same contract as R/registry-probes.R: one probe per equation id in
# inst/model-registry.yml, each driving a Leaf and pairing what the package
# reports with the registry's replica. registry_probes() collects every
# `.registry_probes_*` list in the namespace, so this file needs no other wiring.

# Probes for the hydraulics fragment (sections `vulnerability` and `supply`).
# Sourced with the phylloptim namespace as parent. Each probe reads the package
# through a binding and pairs it with the replica; traits move only through
# leaf_model(leaf_traits()) or set_traits() (hazard 10).

# Stem and root (P50, c) pairs, deliberately different from each other so a
# swapped curve cannot pass (hazard 1).
.hy_stem_sets <- list(c(3.4, 2.680147), c(1.5, 1.2), c(6, 5))
.hy_root_sets <- list(c(3.4, 2.680147), c(2.2, 3.5), c(5, 1.4))
.hy_n <- function() leaf_control()$vulnerability_curve_ncontrol

.hy_traits <- function(st, rt) {
  leaf_traits(stem_P50 = st[1], stem_c = st[2], root_P50 = rt[1], root_c = rt[2])
}

# A single-potential leaf at the given curves and stem conductance.
.hy_single <- function(st = .hy_stem_sets[[1]], rt = .hy_root_sets[[1]], kmax = 1) {
  l <- leaf_model(traits = .hy_traits(st, rt), supply = leaf_supply_singlelayer())
  set_drivers(l, psi_soil = 1, leaf_specific_conductance_max = kmax,
              root_network = series_resistance(1e3))
  l
}

# A multi-layer leaf with a network built directly from resistances.
.hy_network <- function(rH, rV) {
  net <- RootNetwork()
  net$r_R_H_min <- rH
  net$r_R_V_sum <- rV
  net
}
.hy_multi <- function(psi_soil, soil_depth, rH, rV, rt = .hy_root_sets[[1]]) {
  l <- leaf_model(traits = .hy_traits(.hy_stem_sets[[1]], rt))
  set_drivers(l, psi_soil = psi_soil, soil_depth = soil_depth,
              root_network = .hy_network(rH, rV))
  l
}

# Value, or NA for a refusal, so a refused point pairs with a refused replica.
.hy_or_refused <- function(expr) tryCatch(expr, error = function(e) NA_real_)

.registry_probes_hydraulics <- list(

  # --- vulnerability ---------------------------------------------------------
  # Both seating routes: the constructor and set_traits().
  weibull_scale = function(rep, temps) {
    out <- list()
    for (k in seq_along(.hy_stem_sets)) {
      st <- .hy_stem_sets[[k]]; rt <- .hy_root_sets[[k]]
      fresh <- leaf_model(traits = .hy_traits(st, rt))
      moved <- leaf_model()
      set_traits(moved, .hy_traits(st, rt))
      out[[k]] <- .probe_pairs(
        c(fresh$stem_b, fresh$root_b, moved$stem_b, moved$root_b),
        c(rep$weibull_scale(st[1], st[2]), rep$weibull_scale(rt[1], rt[2]),
          rep$weibull_scale(st[1], st[2]), rep$weibull_scale(rt[1], rt[2])))
    }
    do.call(rbind, out)
  },
  # The quantile, and the conductivity the stem curve reports there.
  critical_potential = function(rep, temps) {
    out <- list()
    for (k in seq_along(.hy_stem_sets)) {
      st <- .hy_stem_sets[[k]]; rt <- .hy_root_sets[[k]]
      fresh <- leaf_model(traits = .hy_traits(st, rt))
      moved <- leaf_model()
      set_traits(moved, .hy_traits(st, rt))
      out[[k]] <- .probe_pairs(
        c(fresh$psi_crit, fresh$root_psi_crit, moved$psi_crit, moved$root_psi_crit,
          fresh$proportion_of_conductivity(fresh$psi_crit)),
        c(rep$critical_potential(st[1], st[2]), rep$critical_potential(rt[1], rt[2]),
          rep$critical_potential(st[1], st[2]), rep$critical_potential(rt[1], rt[2]),
          0.05))
    }
    do.call(rbind, out)
  },
  weibull_vulnerability = function(rep, temps) {
    psi <- c(0, 0.1, 0.5, 1, 2, 3.4, 4, 5.87, 6.9, 8, 12, 50)
    do.call(rbind, lapply(.hy_stem_sets, function(st) {
      l <- leaf_model(traits = .hy_traits(st, .hy_root_sets[[1]]))
      .probe_pairs(vapply(psi, l$proportion_of_conductivity, 0),
                   rep$weibull_vulnerability(psi, st[1], st[2]))
    }))
  },
  # The last knot of the root grid, found as the potential past which the
  # clamped root conductivity read stops changing: on the equal-potential branch
  # one layer's flux is -rho g z/(r_H,min/f_r), so it is constant exactly where
  # root_vuln_at clamps. Bisected to adjacent doubles.
  vulnerability_psi_max = function(rep, temps) {
    do.call(rbind, lapply(.hy_root_sets, function(rt) {
      l <- .hy_multi(1, 1, 1e-5, 0, rt = rt)
      flux <- function(T) { l$E_from_Soil_to_Root_Collar(T, T); l$E_up_ }
      flat <- flux(1000)
      lo <- rt[1]; hi <- 1000
      repeat {
        mid <- 0.5 * (lo + hi)
        if (mid <= lo || mid >= hi) break
        if (flux(mid) == flat) hi <- mid else lo <- mid
      }
      .probe_pairs(hi, rep$vulnerability_psi_max(rep$weibull_scale(rt[1], rt[2]), rt[2]))
    }))
  },
  # Where the code reads the closed form: at the knots (the seeds), past P99,
  # and at 1000 MPa (the limit). transpiration(psi, 0) at k_max = 1 is F(psi).
  cumulative_vulnerability_integral = function(rep, temps) {
    n <- .hy_n()
    do.call(rbind, lapply(.hy_stem_sets, function(st) {
      l <- .hy_single(st)
      b <- rep$weibull_scale(st[1], st[2])
      pm <- rep$vulnerability_psi_max(b, st[2])
      i <- c(1, 2, 7, seq(20, n - 20, by = 20), n - 1)
      at <- c(pm * i / n, pm, pm * 1.0001, pm * 1.1, pm * 2, 50, 1000)
      .probe_pairs(vapply(at, function(p) l$transpiration(p, 0), 0),
                   rep$cumulative_vulnerability_integral(at, b, st[2]))
    }))
  },
  # Between the knots, where the interpolant is not the closed form.
  vulnerability_spline = function(rep, temps) {
    n <- .hy_n()
    set.seed(166)
    do.call(rbind, lapply(.hy_stem_sets, function(st) {
      l <- .hy_single(st)
      b <- rep$weibull_scale(st[1], st[2])
      pm <- rep$vulnerability_psi_max(b, st[2])
      at <- c(pm * (seq(0.5, n - 0.5, by = 7)) / n, stats::runif(40, 0, pm))
      .probe_pairs(vapply(at, function(p) l$transpiration(p, 0), 0),
                   rep$vulnerability_spline(b, st[2], n, "integral")(at))
    }))
  },
  # In the grid, past it, refused below zero, and after perturb_stem_P50()
  # (the s != 1 branch), whose spline still describes the old scale.
  stem_vulnerability_integral = function(rep, temps) {
    n <- .hy_n()
    set.seed(167)
    out <- list()
    for (st in .hy_stem_sets) {
      l <- .hy_single(st)
      b <- rep$weibull_scale(st[1], st[2])
      pm <- rep$vulnerability_psi_max(b, st[2])
      at <- c(stats::runif(25, 0, pm), pm * c(0.999, 1, 1.02, 1.5), 30)
      out[[length(out) + 1L]] <- .probe_pairs(
        vapply(at, function(p) l$transpiration(p, 0), 0),
        rep$stem_vulnerability_integral(at, b, st[2], n))
      # Refused below zero, on both sides.
      out[[length(out) + 1L]] <- .probe_pairs(
        is.na(.hy_or_refused(l$transpiration(-0.5, 0))),
        is.na(.hy_or_refused(rep$stem_vulnerability_integral(-0.5, b, st[2], n))))
      for (scale in c(0.9, 1.15)) {
        l$perturb_stem_P50(st[1] * scale)
        b1 <- rep$weibull_scale(st[1] * scale, st[2])
        out[[length(out) + 1L]] <- .probe_pairs(
          c(l$stem_b, vapply(at, function(p) l$transpiration(p, 0), 0)),
          c(b1, rep$stem_vulnerability_integral(at, b1, st[2], n, b0 = b)))
      }
    }
    do.call(rbind, out)
  },
  # The clamped conductivity read, through the equal-potential uptake branch
  # with r_V = 0: E_up = M_w * (-rho g z) / (r_H,min / f_r).
  root_vulnerability = function(rep, temps) {
    n <- .hy_n()
    set.seed(168)
    rH <- 1e-5
    do.call(rbind, lapply(.hy_root_sets, function(rt) {
      l <- .hy_multi(1, 1, rH, 0, rt = rt)
      b <- rep$weibull_scale(rt[1], rt[2])
      pm <- rep$vulnerability_psi_max(b, rt[2])
      at <- c(0, pm * c(1, 50, 133, 399.5) / n, stats::runif(20, 0, pm), pm, pm * 1.05, 40)
      pkg <- vapply(at, function(T) {
        l$E_from_Soil_to_Root_Collar(T, T)
        c(l$E_up_, l$soil_consumption_[1])
      }, c(0, 0))
      Ei <- -(9.8e-3 * 0.5) / (rH / rep$root_vulnerability(at, rt[1], rt[2], n))
      .probe_pairs(as.vector(pkg), as.vector(rbind(Ei * 0.018015, Ei)))
    }))
  },
  # The capped integral, through the general uptake branch with r_V = 0:
  # E_i = (psi_c - psi_s - rho g z) * I / (r_H,min * span). Both ends in the
  # grid, on the linear extension, and past the cap.
  root_vulnerability_integral = function(rep, temps) {
    n <- .hy_n()
    rH <- 1e-4; g <- 9.8e-3 * 0.5
    soils <- c(0.2, 1, 3, 6, 7, 7.2, 9)
    collars <- c(0.5, 2, 5, 6.95, 7.1, 7.6, 12, 40)
    do.call(rbind, lapply(.hy_root_sets, function(rt) {
      l <- .hy_multi(1, 1, rH, 0, rt = rt)
      G <- function(p) rep$root_vulnerability_integral(p, rt[1], rt[2], n)
      gr <- expand.grid(s = soils, c = collars)
      gr <- gr[abs(gr$c - gr$s) > 1e-3 & abs(gr$c - gr$s - g) > 1e-3, ]
      pkg <- mapply(function(s, cc) { l$E_from_Soil_to_Root_Collar(cc, s); l$soil_consumption_[1] },
                    gr$s, gr$c)
      lo <- pmin(gr$s, gr$c); hi <- pmax(gr$s, gr$c)
      mine <- (gr$c - gr$s - g) * (G(hi) - G(lo)) / (rH * (hi - lo))
      .probe_pairs(pkg, mine)
    }))
  },

  # --- supply ----------------------------------------------------------------
  stem_transpiration = function(rep, temps) {
    n <- .hy_n()
    pairs <- expand.grid(s = c(0.3, 1.5, 3, 5.5, 6.8, 7.2, 15), u = c(0, 0.5, 2, 6.5, 8))
    out <- list()
    for (st in .hy_stem_sets[1:2]) {
      b <- rep$weibull_scale(st[1], st[2])
      for (kmax in c(1, 2.5, 3.14e-5)) {
        l <- .hy_single(st, kmax = kmax)
        out[[length(out) + 1L]] <- .probe_pairs(
          mapply(l$transpiration, pairs$s, pairs$u),
          rep$stem_transpiration(pairs$s, pairs$u, kmax, b, st[2], n))
      }
    }
    do.call(rbind, out)
  },
  # At integration_tol = 1e-12, where the quadrature has converged to the closed
  # form; the sweep over the tolerance is in the entry's note.
  stem_transpiration_quadrature = function(rep, temps) {
    pairs <- expand.grid(s = c(0.3, 1.5, 3, 5.5, 6.8, 7.2, 15), u = c(0, 0.5, 2, 6.5, 8))
    do.call(rbind, lapply(.hy_stem_sets, function(st) {
      l <- .hy_single(st, kmax = 1)
      l$initialize_integrator(21L, 1e-12)
      b <- rep$weibull_scale(st[1], st[2])
      .probe_pairs(mapply(l$transpiration_full_integration, pairs$s, pairs$u),
                   rep$stem_transpiration_quadrature(pairs$s, pairs$u, 1, b, st[2]))
    }))
  },
  # Feasible targets pair by value; refused ones (a negative target, or at or
  # past F(inf)) pair as -1 against the replica's own refusal.
  stem_transpiration_inverse = function(rep, temps) {
    n <- .hy_n()
    out <- list()
    for (st in .hy_stem_sets) {
      b <- rep$weibull_scale(st[1], st[2])
      Finf <- (b / st[2]) * gamma(1 / st[2])
      l <- .hy_single(st, kmax = 2)
      g <- expand.grid(w = c(0, 1e-6, 0.01, 0.3, 0.6, 0.9, 0.99, 0.999, 1.01) * Finf,
                       u = c(0, 0.4, 2, 5, 6.85, 7.5, 20))
      g$E <- 2 * g$w
      g <- rbind(g, data.frame(w = NA, u = c(1, 3), E = -2 * Finf))
      pk <- mapply(function(E, u) .hy_or_refused(l$transpiration_to_psi_stem(E, u)), g$E, g$u)
      mi <- mapply(function(E, u) .hy_or_refused(rep$stem_transpiration_inverse(E, u, 2, b, st[2], n)),
                   g$E, g$u)
      out[[length(out) + 1L]] <- .probe_pairs(ifelse(is.na(pk), -1, pk), ifelse(is.na(mi), -1, mi))
    }
    do.call(rbind, out)
  },
  layer_midpoints = function(rep, temps) {
    profiles <- list(0.7, 1:3, c(0.02, 0.3, 1.5), c(0.1, 0.3, 0.6, 1, 1.5, 2))
    do.call(rbind, lapply(profiles, function(d) {
      l <- leaf_model()
      set_drivers(l, psi_soil = rep(1, length(d)), soil_depth = d)
      .probe_pairs(c(l$z_soil_mid_, l$soil_number_of_depths_, l$soil_depth_),
                   c(rep$layer_midpoints(d), length(d), d))
    }))
  },
  # The network the function returns, and the fields the leaf seats from it.
  root_network_from_carbon = function(rep, temps) {
    cases <- list(
      list(C = rep(20 / 3, 3), d = 1:3, bH = 340, bV = 9400),
      list(C = c(2, 3, 0, 1, 0, 0), d = c(0.1, 0.3, 0.6, 1, 1.5, 2), bH = 340, bV = 9400),
      list(C = c(5, 1e-3), d = c(0.02, 1.5), bH = 100, bV = 2e4),
      list(C = c(0.4, 0.4, 0.4, 0), d = c(0.25, 0.5, 1, 3), bH = 50, bV = 500))
    fields <- c("c_r_V", "c_r_H", "r_R_H_min", "r_R_V", "r_R_V_sum")
    do.call(rbind, lapply(cases, function(k) {
      net <- root_network_from_carbon(k$C, k$d, k$bH, k$bV)
      l <- leaf_model()
      set_drivers(l, psi_soil = rep(1, length(k$d)), soil_depth = k$d, root_network = net)
      m <- rep$root_network_from_carbon(k$C, k$d, k$bH, k$bV)
      .probe_pairs(
        c(unlist(lapply(fields, function(f) net[[f]])),
          l$c_r_V_, l$c_r_H_, l$r_R_H_min, l$r_R_V, l$r_R_V_sum, l$max_soil_layer),
        c(unlist(m[fields]), unlist(m[fields]), m$max_soil_layer))
    }))
  },
  # Small resistances so fluxes are O(1); unequal, unrooted and zero-resistance
  # layers; collars wetter and drier than each layer, equal to one (the
  # equal-potential branch), at one's gravity balance, below zero, and past the
  # root grid. Plus the carbon-derived default network.
  multilayer_uptake = function(rep, temps) {
    n <- .hy_n()
    cases <- list(
      list(psi = c(0.5, 2, 4, 6.5), d = c(0.1, 0.4, 1, 2),
           rH = c(0.01, 0.02, 0), rV = c(0.005, 0.01, 0.02), rt = .hy_root_sets[[1]]),
      list(psi = c(1, 1.5, 3), d = c(0.3, 0.6, 3),
           rH = c(0.05, 0.002, 0.03), rV = c(0, 0.04, 0.05), rt = .hy_root_sets[[3]]))
    out <- list()
    for (k in cases) {
      l <- .hy_multi(k$psi, k$d, k$rH, k$rV, rt = k$rt)
      z1 <- (k$d[1]) / 2
      collars <- c(-0.2, 0, 0.3, k$psi[1], k$psi[1] + 9.8e-3 * z1, 1.2, k$psi[2],
                   2.7, 5, 7, 7.5, 12)
      for (cc in collars) {
        l$E_from_Soil_to_Root_Collar(cc, k$psi)
        out[[length(out) + 1L]] <- .probe_pairs(
          c(l$E_up_, l$soil_consumption_),
          rep$multilayer_uptake(cc, k$psi, k$d, k$rH, k$rV, k$rt[1], k$rt[2], n))
      }
    }
    d <- c(0.5, 1, 2)
    net <- root_network_from_carbon(c(1, 1, 0), d)
    l <- leaf_model()
    set_drivers(l, psi_soil = c(0.8, 1.6, 2.4), soil_depth = d, root_network = net)
    for (cc in c(0.5, 1.2, 2, 3)) {
      l$E_from_Soil_to_Root_Collar(cc, c(0.8, 1.6, 2.4))
      out[[length(out) + 1L]] <- .probe_pairs(
        c(l$E_up_, l$soil_consumption_),
        rep$multilayer_uptake(cc, c(0.8, 1.6, 2.4), d, net$r_R_H_min, net$r_R_V_sum,
                              3.4, 2.680147, n))
    }
    do.call(rbind, out)
  },
  single_series_supply = function(rep, temps) {
    g <- expand.grid(h = c(0, 0.05), r = c(800, 0.01), s = c(0.3, 2))
    do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
      l <- leaf_model(supply = leaf_supply_singlelayer(g$h[i]))
      set_drivers(l, psi_soil = g$s[i], root_network = series_resistance(g$r[i]))
      do.call(rbind, lapply(c(0, 0.2, g$s[i], g$s[i] + g$h[i], 1.1, 3.5), function(cc) {
        l$E_from_Soil_to_Root_Collar(cc, g$s[i])
        .probe_pairs(c(l$E_up_, l$soil_consumption_, l$single_psi_soil_,
                       l$single_resistance_, l$single_gravity_head_),
                     c(rep$single_series_supply(cc, g$s[i], g$r[i], g$h[i]),
                       g$s[i], g$r[i], g$h[i]))
      }))
    }))
  }
)
