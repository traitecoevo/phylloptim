# The finite-difference referee: the trait gradient against the model.
#
# Every (curve, route) pair, on both supplies where the route allows it, in the
# parameters its curve owns and in those every curve reads, against a central
# difference of its whole solve. The difference is taken through leaf_solve(),
# which shares nothing with the gradient but the model, so a parameter the
# gradient fails to write (as `LeastCost_beta` once was) is a disagreement of
# order one here, not a plausible zero. gradient_golden.tsv cannot see that: it
# is generated from the gradient and would have recorded the zero as correct.
#
# The step is swept and the best agreement is read (hazard 15): error against
# step is a V, and only its floor is a derivative. An output the curve holds
# invariant (CF77's A in a stem trait: its optimum is where dA/dE meets the
# price, whatever the stem) is noise on both sides, so agreement is measured
# against a floor of 1e-5 of the output per unit of the parameter.
#
# Measured on macOS/arm64: worst 5.2e-4 (CF77's g_min on the leaf route, an
# output at that floor), 99th percentile 1.3e-4, median 2e-9, over 195 rows.

referee_outputs <- c("A", "gc", "psi_leaf", "collar", "profit")

# Interior on the multi-layer collar route for every curve; at their aligned
# defaults CF77 and LeastCost pin there (see helper-gradient.R).
referee_prices <- list(CF77 = list(CF77_lambda = 3e5),
                       TF24_floor = list(TF24_floor_lambda_o = 1.5e5),
                       LeastCost = list(LeastCost_beta = 5))

# g_min away from zero: centred at zero, the difference would step onto a
# negative residual conductance, which the leaf refuses.
referee_traits <- leaf_traits(g_min = 0.002)

referee_setups <- list(
  collar_multi = list(supply = leaf_supply_multilayer(), route = "collar",
                      network = root_network_from_carbon(1 / 0.05,
                                                         soil_depth = 1),
                      depth = 1),
  collar_single = list(supply = leaf_supply_singlelayer(), route = "collar",
                       resistance = 1e4),
  leaf_single = list(supply = leaf_supply_singlelayer(), route = "leaf",
                     resistance = 1e4))

# The drivers both sides take, for a setup and a resistance.
referee_drivers <- function(su, resistance = su$resistance) {
  d <- list(psi_soil = 1.0, PPFD = 900, atm_vpd = 2)
  if (is.null(su$resistance)) {
    d$soil_depth <- su$depth
    d$root_network <- su$network
  } else {
    d$root_network <- series_resistance(resistance)
  }
  d
}

# The five outputs, solved by leaf_solve() with one parameter set to `value`.
referee_solve <- function(curve, su, par = NULL, value = NULL) {
  tr <- referee_traits
  prices <- referee_prices[[curve]]
  kmax <- kmax_default
  resistance <- su$resistance
  if (!is.null(par)) {
    if (par %in% names(tr)) {
      tr[[par]] <- value
    } else if (par == "leaf_specific_conductance_max") {
      kmax <- value
    } else if (par == "resistance") {
      resistance <- value
    } else {
      prices[[phylloptim:::.owned_argument(par)]] <- value
    }
  }
  args <- c(referee_drivers(su, resistance),
            list(leaf_specific_conductance_max = kmax, traits = tr,
                 supply = su$supply, model = curve, route = su$route),
            prices)
  unlist(do.call(leaf_solve, args)[1, referee_outputs])
}

# A parameter's base value, and the scale its step is relative to: the
# gradient's own rule (relative, floored at 1 but for the two drivers).
referee_value <- function(curve, su, par) {
  if (par %in% names(referee_traits)) return(referee_traits[[par]])
  if (par == "leaf_specific_conductance_max") return(kmax_default)
  if (par == "resistance") return(su$resistance)
  given <- referee_prices[[curve]][[phylloptim:::.owned_argument(par)]]
  if (!is.null(given)) given else .level_defaults()[[par]]
}

referee_scale <- function(par, value) {
  if (par %in% c("leaf_specific_conductance_max", "resistance")) {
    abs(value)
  } else {
    max(abs(value), 1)
  }
}

test_that("every curve's gradient on every route agrees with a difference of its solve", {
  params <- cost_curve_parameters()
  rows <- list()
  for (curve in cost_curve_names()) {
    for (sn in names(referee_setups)) {
      su <- referee_setups[[sn]]
      pars <- c("vcmax_25", "stem_P50", "stem_c", "R_d_25", "g_min",
                "leaf_specific_conductance_max",
                params$parameter[params$curve == curve],
                if (identical(su$route, "collar") && !is.null(su$resistance))
                  "resistance")
      x <- leaf_model(referee_traits, supply = su$supply)
      for (f in names(referee_prices[[curve]])) {
        x[[if (f == "CF77_lambda") "CF77_lambda_" else f]] <-
          referee_prices[[curve]][[f]]
      }
      g <- do.call(leaf_gradient,
                   c(referee_drivers(su),
                     list(x = x, traits = referee_traits, pars = pars,
                          model = curve, route = su$route)))
      y0 <- referee_solve(curve, su)
      for (p in pars) {
        v0 <- referee_value(curve, su, p)
        scale <- referee_scale(p, v0)
        floor_y <- 1e-5 * abs(y0) / scale
        gp <- g$gradient[p, referee_outputs]
        best <- rep(Inf, length(referee_outputs))
        steps <- 0L
        for (k in 2:6) {
          h <- scale * 10^-k
          # g_min is refused below zero (leaf_traits() checks it; the edit in
          # referee_solve() does not), so its stencil stays above zero.
          if (p == "g_min" && h >= v0) next
          # A large step can leave the range the solve brackets; the leaf
          # refuses that as infeasible, and the step is skipped.
          fd <- tryCatch(
            (referee_solve(curve, su, p, v0 + h) -
             referee_solve(curve, su, p, v0 - h)) / (2 * h),
            phylloptim_infeasible = function(e) NULL)
          if (is.null(fd)) next
          steps <- steps + 1L
          err <- abs(gp - fd) / pmax(abs(gp), abs(fd), floor_y)
          # A leaf route's collar is not solved for: NaN on both sides.
          err[!is.finite(gp) & !is.finite(fd)] <- 0
          best <- pmin(best, err)
        }
        rows[[length(rows) + 1L]] <- data.frame(
          curve = curve, setup = sn, status = g$status, par = p, steps = steps,
          worst = max(best), output = referee_outputs[which.max(best)])
      }
    }
  }
  d <- do.call(rbind, rows)
  w <- d[which.max(d$worst), ]
  cat(sprintf(paste0("\nreferee: %d pairs, %d rows (%d interior, %d pinned); ",
                     "worst %.3g (%s %s, %s in %s), median %.3g\n"),
              nrow(unique(d[c("curve", "setup")])), nrow(d),
              sum(d$status == "interior"), sum(d$status == "pinned"),
              w$worst, w$curve, w$setup, w$output, w$par, stats::median(d$worst)))
  # Hazard 16: the rows must reach the assertion, and the composite must run on
  # most of them rather than the finite-difference arm.
  expect_identical(nrow(unique(d[c("curve", "setup")])), 27L)
  expect_gte(sum(d$status == "interior"), 170L)
  expect_gte(min(d$steps), 3L)
  bad <- d[d$worst > 5e-3, ]
  expect_identical(nrow(bad), 0L,
                   label = paste("rows off the solve:",
                                 paste(bad$curve, bad$setup, bad$par,
                                       collapse = "; ")))
})
