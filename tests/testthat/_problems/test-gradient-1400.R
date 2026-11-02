# Extracted from test-gradient.R:1400

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
sp <- leaf_supply_singlelayer()
net <- series_resistance(1e4)
outs <- c("A", "gc", "psi_stem", "collar", "profit")
g0 <- leaf_gradient(1.5, pars = "vcmax_25", supply = leaf_supply_multilayer(), model = "TF24", route = "collar")
expect_identical(c(g0$model, g0$route), c("TF24_floor", "collar"))
