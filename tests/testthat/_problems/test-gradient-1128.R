# Extracted from test-gradient.R:1128

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_true("TF24_floor_lambda_o" %in% gradient_par_names())
expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900, model = "TF24",
                             pars = "TF24_floor_lambda_o", supply = leaf_supply_multilayer()),
               "not a parameter of the TF24 model")
