# Extracted from test-gradient.R:924

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
d <- grid_drivers(2.0)
a <- do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25")))
psi_off <- a$psi + 0.2
g <- do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25", psi = psi_off)))
expect_identical(g$status, "prescribed")
expect_gt(g$stationarity, 1e-3)
l <- leaf_model(leaf_traits(), leaf_control(), leaf_supply_multilayer())
v <- leaf_traits()$vcmax_25
h <- max(abs(v), 1) * 1e-6
at <- function(x, field) {
    set_traits(l, leaf_traits(vcmax_25 = x))
    do.call(grid_set_drivers, c(list(l), d))
    l$evaluate_root_collar_psi(psi_off)
    l[[field]]
  }
