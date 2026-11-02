# Extracted from test-gradient.R:54

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
l <- leaf_model(supply = leaf_supply_multilayer())
do.call(set_drivers, c(list(l), grid_drivers(2.0)))
