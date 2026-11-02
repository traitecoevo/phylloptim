# Extracted from test-gradient.R:503

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
dry <- list(psi_soil = 5.5, PPFD = 900, atm_vpd = 3.0,
              supply = leaf_supply_singlelayer(),
              root_network = series_resistance(1e4))
auto <- do.call(leaf_gradient,
                  c(dry, list(pars = c("leaf_specific_conductance_max",
                                       "resistance"))))
expect_true(auto$status %in% c("pinned", "no-gradient"))
expect_identical(auto$method, "fd")
expect_true(all(is.finite(auto$gradient)))
