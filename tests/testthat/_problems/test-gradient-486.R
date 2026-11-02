# Extracted from test-gradient.R:486

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
            supply = leaf_supply_singlelayer(),
            root_network = series_resistance(1e4))
for (p in c("leaf_specific_conductance_max", "resistance")) {
    vals <- vapply(10^-(4:7), function(st) {
      do.call(leaf_gradient, c(d, list(pars = p, step = st,
                                       method = "ift")))$gradient[p, "A"]
    }, numeric(1))
    expect_lt(diff(range(vals)) / abs(mean(vals)), 1e-6)
  }
