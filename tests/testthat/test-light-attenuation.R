# light_depth_ > 0: a leaf whose layers share c_i, light falling as exp(-t).

attenuated_leaf <- function(depth, theta_layer = 0.9, ppfd = 800, ci_fixed = NULL) {
  l <- leaf_model(leaf_traits(curv_fact_elec_trans = theta_layer, R_d_25 = 0))
  l$light_depth_ <- depth
  set_drivers(l, psi_soil = 0.5, PPFD = ppfd)
  l
}

# The same leaf by brute force: layers of a non-rectangular hyperbola, summed.
layer_sum <- function(A_sat, A_L, tau, theta, n = 20000) {
  t <- (seq_len(n) - 0.5) / n * tau
  top <- A_L / (1 - exp(-tau))
  light <- top * exp(-t); s <- A_sat / tau
  sum((light + s - sqrt((light + s)^2 - 4 * theta * light * s)) / (2 * theta)) * tau / n
}

test_that("attenuation is off by default and leaves the solve unchanged", {
  l0 <- leaf_model(); set_drivers(l0, psi_soil = 0.5, PPFD = 800); l0$optimise()
  expect_equal(l0$light_depth_, 0)
  l1 <- leaf_model(); l1$light_depth_ <- 0; set_drivers(l1, psi_soil = 0.5, PPFD = 800); l1$optimise()
  expect_identical(operating_point(l0), operating_point(l1))
})

test_that("the attenuated kernel is the integral over layers", {
  for (depth in c(0.5, 2, 5)) for (ppfd in c(50, 400, 2000)) {
    l <- attenuated_leaf(depth, theta_layer = 0.9, ppfd = ppfd)
    ci <- 25
    gs <- l$gamma_ * l$atm_kpa_ * 1e-3
    wj <- (ci - gs) / (4 * (ci + 2 * gs)); wc <- (ci - gs) / (ci + l$km_)
    a <- l$vcmax_ * wc; m <- l$jmax_ * wj; th <- l$curv_fact_colim
    A_sat <- (a + m - sqrt((a + m)^2 - 4 * th * a * m)) / (2 * th)
    A_L <- l$phi_J * l$alpha_PAR_ * ppfd * wj
    expect_equal(l$assim_colimited(ci), layer_sum(A_sat, A_L, depth, 0.9), tolerance = 1e-6,
                 info = paste("depth", depth, "PPFD", ppfd))
  }
})

test_that("respiration is subtracted once, after the integral over layers", {
  l <- leaf_model(leaf_traits(curv_fact_elec_trans = 0.9))
  l$light_depth_ <- 2
  set_drivers(l, psi_soil = 0.5, PPFD = 600)
  expect_gt(l$R_d_, 0)
  ci <- 25
  gs <- l$gamma_ * l$atm_kpa_ * 1e-3
  wj <- (ci - gs) / (4 * (ci + 2 * gs)); wc <- (ci - gs) / (ci + l$km_)
  a <- l$vcmax_ * wc; m <- l$jmax_ * wj; th <- l$curv_fact_colim
  A_sat <- (a + m - sqrt((a + m)^2 - 4 * th * a * m)) / (2 * th)
  A_L <- l$phi_J * l$alpha_PAR_ * 600 * wj
  expect_equal(l$assim_colimited(ci), layer_sum(A_sat, A_L, 2, 0.9) - l$R_d_, tolerance = 1e-6)
})

test_that("depth concentrates light in the upper layers", {
  ci <- 25
  # in weak light every layer is light-limited and the leaf uses all it absorbs
  thin <- attenuated_leaf(0.02, ppfd = 20)
  gs <- thin$gamma_ * thin$atm_kpa_ * 1e-3
  wj <- (ci - gs) / (4 * (ci + 2 * gs))
  expect_equal(thin$assim_colimited(ci), thin$phi_J * thin$alpha_PAR_ * 20 * wj, tolerance = 0.02)
  # at fixed capacity and absorbed light, a deeper stack puts more of the light
  # on its upper layers, which saturate: the rate falls with depth
  r <- vapply(c(0.5, 2, 6), function(d) attenuated_leaf(d, ppfd = 1500)$assim_colimited(ci), 0)
  expect_true(all(diff(r) < 0))
})

test_that("the solve runs with attenuation on, and the closed form refuses it", {
  l <- attenuated_leaf(3, ppfd = 1200)
  l$optimise()
  op <- operating_point(l)
  expect_true(is.finite(op$A) && op$A > 0)
  expect_true(op$ci > 0 && op$ci < 40)
  c <- leaf_model(leaf_traits(curv_fact_elec_trans = 0.9))
  c$light_depth_ <- 3
  set_drivers(c, psi_soil = 0.5, PPFD = 1200)
  c$set_model("TF24", "leaf", "closed")
  expect_error(c$optimise(), "light_depth_")
})

test_that("bad settings are refused by name", {
  expect_error(attenuated_leaf(-1), "light_depth_")
  expect_error(attenuated_leaf(NaN), "light_depth_")
  expect_error(attenuated_leaf(Inf), "light_depth_")
  expect_error(attenuated_leaf(2, theta_layer = 1), "curv_fact_elec_trans")
  # set after the drivers, a bad depth is refused when the leaf is solved
  l <- leaf_model(); set_drivers(l, psi_soil = 0.5, PPFD = 800)
  l$light_depth_ <- -1
  expect_error(l$optimise(), "light_depth_")
  # the single layer's electron-limited rate is not the attenuated leaf's
  expect_error(attenuated_leaf(2)$assim_electron_limited(25), "light_depth_")
})

test_that("the trait gradient matches a difference of the solve with attenuation on", {
  # every trait the attenuated kernel reads: the capacities, the light (phi_J),
  # each layer's curvature and respiration
  pars <- c("vcmax_25", "jmax_25", "phi_J", "curv_fact_elec_trans", "curv_fact_colim",
            "R_d_25")
  base <- leaf_traits(curv_fact_elec_trans = 0.9, curv_fact_colim = 0.95)
  solve_A <- function(tr) {
    l <- leaf_model(tr)
    l$light_depth_ <- 3
    set_drivers(l, psi_soil = 0.5, PPFD = 900)
    l$optimise()
    operating_point(l)$A
  }
  l <- leaf_model(base); l$light_depth_ <- 3
  g <- leaf_gradient(psi_soil = 0.5, PPFD = 900, x = l, traits = base, pars = pars)
  for (p in pars) {
    h <- 1e-3 * base[[p]]
    up <- base; up[[p]] <- base[[p]] + h
    dn <- base; dn[[p]] <- base[[p]] - h
    fd <- (solve_A(up) - solve_A(dn)) / (2 * h)
    expect_true(fd != 0, info = p)
    expect_equal(g$gradient[p, "A"], fd, tolerance = 1e-4, info = p)
  }
})
