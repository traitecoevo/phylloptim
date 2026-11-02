# Extracted from test-gradient.R:749

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
tr <- leaf_traits(vcmax_25 = 105, stem_P50 = 4.2)
args <- list(psi_soil = 2.0, PPFD = 900, traits = tr, pars = "vcmax_25")
fresh <- do.call(leaf_gradient, args)
wrong <- leaf_model(traits = leaf_traits(vcmax_25 = 60, stem_P50 = 4.6), supply = leaf_supply_multilayer())
reused <- do.call(leaf_gradient, c(args, list(x = wrong)))
