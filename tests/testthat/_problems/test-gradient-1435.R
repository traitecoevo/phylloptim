# Extracted from test-gradient.R:1435

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
expect_identical(leaf_gradient(1.5, supply = sp, root_network = net,
                                 pars = "vcmax_25", model = "SOX")$route,
                   "collar")
expect_identical(leaf_solve(1.5, supply = sp, root_network = net,
                              model = "SOX"),
                   leaf_solve(1.5, supply = sp, root_network = net,
                              model = "SOX", route = "collar"))
tf24 <- leaf_solve(1.5, model = "TF24", route = "collar", supply = leaf_supply_multilayer())
expect_identical(leaf_solve(1.5, supply = leaf_supply_multilayer(), model = "TF24", route = "collar"), tf24)
expect_identical(leaf_solve(1.5, model = "TF24", route = "collar", supply = leaf_supply_multilayer()), tf24)
expect_identical(
    g0$gradient,
    leaf_gradient(1.5, pars = "vcmax_25", model = "TF24",
                  route = "collar", supply = leaf_supply_multilayer())$gradient)
for (m in c("JS22", "SOX", "ProfitMax")) {
    g <- leaf_gradient(1.5, pars = "vcmax_25", model = m, route = "collar", supply = leaf_supply_multilayer())
    r <- leaf_solve(1.5, model = m, route = "collar", supply = leaf_supply_multilayer())
    expect_identical(unname(g$value[outs]), unlist(r[1, outs], use.names = FALSE),
                     label = m)
    expect_identical(g$status, "interior", label = m)
    expect_true(is.finite(g$gradient["vcmax_25", "collar"]), label = m)
  }
b <- leaf_batch(psi_soil = c(1.0, 1.5), supply = leaf_supply_multilayer())
expect_identical(leaf_gradient_batch(b, pars = "vcmax_25", model = "SOX",
                                       route = "collar")$status,
                   c("interior", "interior"))
expect_error(leaf_gradient(1.5, pars = "vcmax_25", model = "TF24", route = "collar",
                             route = "stem", supply = leaf_supply_multilayer()),
               "model = \"TF24\"")
