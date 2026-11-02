# Extracted from test-gradient.R:466

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
            supply = leaf_supply_singlelayer(),
            root_network = series_resistance(1e4))
solve_at <- function(par, value) {
    a <- d
    # `resistance` is a driver now, so perturbing it means rebuilding the
    # root_network driver rather than rebuilding the supply object -- which is the
    # whole of the consistency change, seen from the caller's side.
    if (identical(par, "resistance")) {
      a$root_network <- series_resistance(value)
    } else {
      a[[par]] <- value
    }
    x <- do.call(leaf_solve, a)
    c(A = x$A, gc = x$gc, psi_stem = x$psi_stem, collar = x$collar,
      profit = x$profit)
  }
g <- do.call(leaf_gradient,
               c(d, list(pars = c("leaf_specific_conductance_max",
                                  "resistance"), method = "ift")))
expect_identical(g$status, "interior")
for (p in c("leaf_specific_conductance_max", "resistance")) {
    v <- if (identical(p, "resistance")) 1e4 else 3.14e-5
    h <- v * 1e-5
    ref <- (solve_at(p, v + h) - solve_at(p, v - h)) / (2 * h)
    # Six digits, which is well inside the reference's own precision. Measured
    # agreement is 8 digits for the driver and 9 for the resistance.
    expect_equal(g$gradient[p, ], ref, tolerance = 1e-6)
    # And not accidentally zero, which is how this would fail if the setter did
    # nothing at all -- the failure mode psi_crit genuinely has (see above).
    expect_gt(abs(g$gradient[p, "A"]), 1e-8)
  }
