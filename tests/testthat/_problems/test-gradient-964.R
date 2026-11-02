# Extracted from test-gradient.R:964

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- grid_drivers(2.0)
a <- do.call(leaf_gradient, c(d, list(pars = "vcmax_25")))
