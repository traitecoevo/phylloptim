# Extracted from test-gradient.R:1278

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
tr <- leaf_traits()
net <- series_resistance(1e4)
l <- leaf_model(tr, leaf_control(), leaf_supply_singlelayer())
collar <- leaf_gradient(psi_soil = 1.5, PPFD = 1500, pars = "TF24_cost_scale",
                          model = "TF24", supply = leaf_supply_multilayer())
