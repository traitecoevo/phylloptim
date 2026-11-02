# Extracted from test-gradient.R:1179

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
nms <- cost_curve_names()
expect_identical(nms, c("TF24", "CF77", "JS22", "CMax", "SOX", "JW26",
                          "ProfitMax", "TF24_floor", "LeastCost"))
expect_false(exists("cost_curve_has_derivative",
                      where = asNamespace("phylloptim"), inherits = FALSE))
expect_identical(phylloptim:::.gradient_models(), c("collar", nms))
