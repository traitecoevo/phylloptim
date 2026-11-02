# Extracted from test-gradient.R:1118

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
expect_true("CF77_lambda_" %in% gradient_par_names())
i <- match("CF77_lambda_", gradient_par_names())
expect_identical(gradient_par_names()[i + 0:2],
                   c("CF77_lambda_", "TF24_floor_lambda_o", "LeastCost_beta"))
expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900,
                             pars = "CF77_lambda_", supply = leaf_supply_multilayer(), model = "TF24", route = "collar"),
               "not a parameter of the TF24_floor model")
