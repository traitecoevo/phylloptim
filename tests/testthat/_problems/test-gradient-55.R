# Extracted from test-gradient.R:55

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
traits <- leaf_traits(vcmax_25 = 110, stem_P50 = 4.2, root_P50 = 3.5,
                        TF24_cost_scale = 8.0)
d <- grid_drivers(2.0)
fresh <- leaf_model(traits, supply = leaf_supply_multilayer())
do.call(grid_set_drivers, c(list(fresh), d))
