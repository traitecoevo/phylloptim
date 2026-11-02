# Extracted from test-gradient.R:32

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
traits <- leaf_traits(vcmax_25 = 110, stem_P50 = 4.2, root_P50 = 3.5,
                        TF24_cost_scale = 8.0)
d <- grid_drivers(2.0)
fresh <- leaf_model(traits, supply = leaf_supply_multilayer())
do.call(set_drivers, c(list(fresh), d))
