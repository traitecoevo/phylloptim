# Extracted from test-gradient.R:764

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
tr <- leaf_traits()
args <- list(psi_soil = 1.0, PPFD = 900, traits = tr,
               pars = c("vcmax_25", "stem_P50"))
fresh <- do.call(leaf_gradient, args)
seated <- leaf_model(traits = tr, supply = leaf_supply_multilayer())
seated$set_model("CF77", "collar")
reused <- do.call(leaf_gradient, c(args, list(x = seated)))
