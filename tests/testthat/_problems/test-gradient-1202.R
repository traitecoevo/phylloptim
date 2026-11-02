# Extracted from test-gradient.R:1202

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "phylloptim", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
tf24_collar <- function(f) {
  function(...) {
    a <- list(...)
    if (is.null(a$model)) {
      a$model <- "TF24"
      if (is.null(a$route)) a$route <- "collar"
    }
    if (is.null(a$supply) && is.null(a$x)) a$supply <- leaf_supply_multilayer()
    do.call(f, a)
  }
}
grid_leaf_gradient <- tf24_collar(leaf_gradient)
grid_leaf_solve <- tf24_collar(leaf_solve)
grid_set_drivers <- function(l, ...) {
  a <- list(...)
  a$supply <- NULL
  do.call(grid_set_drivers, c(list(l), a))
}

# test -------------------------------------------------------------------------
nms <- cost_curve_names()
expect_identical(nms, c("TF24", "CF77", "JS22", "CMax", "SOX", "JW26",
                          "ProfitMax", "TF24_floor", "LeastCost"))
expect_false(exists("cost_curve_has_derivative",
                      where = asNamespace("phylloptim"), inherits = FALSE))
expect_identical(phylloptim:::.gradient_models(), c("collar", nms))
