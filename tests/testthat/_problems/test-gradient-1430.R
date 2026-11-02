# Extracted from test-gradient.R:1430

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
sp <- leaf_supply_singlelayer()
net <- series_resistance(1e4)
outs <- c("A", "gc", "psi_stem", "collar", "profit")
g0 <- leaf_gradient(1.5, pars = "vcmax_25", supply = leaf_supply_multilayer(), model = "TF24", route = "collar")
expect_identical(c(g0$model, g0$route), c("TF24_floor", "collar"))
expect_identical(leaf_gradient(1.5, supply = sp, root_network = net,
                                 pars = "vcmax_25", model = "SOX")$route,
                   "collar")
expect_identical(leaf_solve(1.5, supply = sp, root_network = net,
                              model = "SOX"),
                   leaf_solve(1.5, supply = sp, root_network = net,
                              model = "SOX", route = "collar"))
