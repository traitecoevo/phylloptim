# Extracted from test-gradient.R:297

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
g <- grid_gradient(16.0, pars = c("vcmax_25", "stem_P50", "root_P50", "R_d_25"))
expect_identical(g$status, "no-gradient")
expect_identical(g$method, "fd")
expect_lt(abs(g$value[["A"]]), 1e-12)
expect_lt(abs(g$gradient["vcmax_25", "A"]), 1e-9)
expect_lt(abs(g$gradient["R_d_25", "A"]), 1e-9)
expect_lt(abs(g$gradient["R_d_25", "profit"]), 1e-9)
l <- leaf_model(supply = leaf_supply_multilayer())
do.call(grid_set_drivers, c(list(l), grid_drivers(16.0)))
