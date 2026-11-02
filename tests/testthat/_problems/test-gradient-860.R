# Extracted from test-gradient.R:860

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- grid_drivers(2.0)
pars <- c("vcmax_25", "stem_P50", "TF24_cost_scale")
a <- do.call(leaf_gradient, c(d, list(pars = pars)))
