# Extracted from test-gradient.R:722

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
tr <- leaf_traits(vcmax_25 = 105, stem_P50 = 4.2)
args <- list(psi_soil = 2.0, PPFD = 900, atm_vpd = 1.5, traits = tr,
               pars = c("vcmax_25", "stem_P50", "TF24_cost_scale"))
fresh <- do.call(leaf_gradient, args)
l <- leaf_model(traits = tr, supply = leaf_supply_multilayer())
reused <- do.call(leaf_gradient, c(args, list(x = l)))
