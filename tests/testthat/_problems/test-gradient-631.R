# Extracted from test-gradient.R:631

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- grid_drivers(0.5, ppfd = 1500, vpd = 2.0, layers = 1L)
g <- do.call(leaf_gradient, c(d, list(pars = "vcmax_25")))
