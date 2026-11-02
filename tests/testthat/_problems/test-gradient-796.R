# Extracted from test-gradient.R:796

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
tr <- leaf_traits()
args <- list(psi_soil = 1.0, PPFD = 900, traits = tr,
               pars = c("vcmax_25", "stem_P50"))
fresh <- do.call(grid_leaf_gradient, args)
seated <- leaf_model(traits = tr, supply = leaf_supply_multilayer())
seated$set_model("CF77", "collar")
reused <- do.call(grid_leaf_gradient, c(args, list(x = seated)))
expect_identical(fresh$status, "interior")
expect_identical(reused$status, fresh$status)
expect_identical(reused$gradient, fresh$gradient)
expect_identical(reused$value, fresh$value)
tf <- leaf_model(traits = tr, supply = leaf_supply_multilayer())
set_drivers(tf, psi_soil = 1.0, PPFD = 900)
tf$optimise()
expect_identical(operating_point(seated), operating_point(tf))
