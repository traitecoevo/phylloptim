# Extracted from test-gradient.R:274

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
g <- grid_gradient(16.0, pars = c("vcmax_25", "stem_P50", "root_P50", "R_d_25"))
expect_identical(g$status, "no-gradient")
expect_identical(g$method, "fd")
expect_lt(abs(g$value[["A"]]), 1e-12)
expect_lt(abs(g$gradient["vcmax_25", "A"]), 1e-9)
expect_lt(abs(g$gradient["R_d_25", "A"]), 1e-9)
expect_lt(abs(g$gradient["R_d_25", "profit"]), 1e-9)
l <- leaf_model(supply = leaf_supply_multilayer())
do.call(set_drivers, c(list(l), grid_drivers(16.0)))
