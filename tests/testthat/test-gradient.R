# Trait gradients: set_traits() and leaf_gradient() (issue #4, PLAN 11d stage 1).
#
# Two things are being pinned here, and they are not equally hard.
#
#   1. That the implicit-function composite is right where its premise holds.
#      This is the easy half: the reference values below were arbitrated against
#      a least-squares slope of the solved output over a +-2% trait span at
#      n = 41, and the composite agrees to about four digits.
#   2. That it is NOT used where its premise fails. This is the half that took
#      the work. At a wet-pinned optimum the bare composite returns O(1) where
#      the truth is ~1e-08 -- plausible-looking and wrong by seven orders of
#      magnitude -- so the test that matters is the one asserting the guard
#      fires. `method = "ift"` exists so that failure can be provoked on demand
#      rather than described in a comment.

# The golden grid's drivers, so a row here is the same operating point the C++
# suite and test-golden.R use. Taken from tests/cpp/test_golden.cpp.
grid_drivers <- function(psi_soil, ppfd = 900, vpd = 2.0, layers = 1L) {
  theta <- 0.000157
  area_leaf <- 0.05
  list(psi_soil = psi_soil + 0.25 * (seq_len(layers) - 1),
       PPFD = ppfd,
       soil_depth = 1.0 * seq_len(layers),
       root_carbon_per_leaf_area = rep(1 / layers / area_leaf, layers),
       leaf_specific_conductance_max = 1.0 * theta / 5.0,
       atm_vpd = vpd, ca = 40, leaf_temp = 25, atm_o2_kpa = 21,
       atm_kpa = 101.3)
}

grid_gradient <- function(psi_soil, ppfd = 900, vpd = 2.0, layers = 1L, ...) {
  do.call(leaf_gradient, c(grid_drivers(psi_soil, ppfd, vpd, layers), list(...)))
}

test_that("set_traits() on a used leaf equals a leaf built with those traits", {
  # The R-side statement of the C++ suite's test_set_traits_matches_a_fresh_leaf,
  # made at the boundary a calibration loop actually crosses. Bit-exact, because
  # the two routes share no code: any derived state set_traits fails to refresh
  # -- the two vulnerability splines, or vcmax_/jmax_/R_d_ behind
  # set_physiology's temperature cache -- shows up as a difference.
  traits <- leaf_traits(vcmax_25 = 110, stem_b = 4.2, root_b = 3.5,
                        cost_scale_TF24 = 8.0)
  d <- grid_drivers(2.0)

  fresh <- leaf_model(traits)
  do.call(set_drivers, c(list(fresh), d))
  fresh$find_root_collar_psi()

  # Solved at the DEFAULTS first, so every cache is warm and pointing at the old
  # traits. Solving first is what makes this a test rather than a coincidence.
  reused <- leaf_model()
  do.call(set_drivers, c(list(reused), d))
  reused$find_root_collar_psi()
  set_traits(reused, traits)
  do.call(set_drivers, c(list(reused), d))
  reused$find_root_collar_psi()

  expect_identical(operating_point(reused), operating_point(fresh))
  expect_identical(reused$vcmax_, fresh$vcmax_)
  expect_identical(reused$R_d_, fresh$R_d_)
})

test_that("set_traits() requires the drivers to be set again", {
  # Not a nicety: the derived photosynthetic parameters really are unknown until
  # set_physiology runs, so the operating point must not be readable as if it
  # were still valid. setup_clean_leaf puts the NA sentinels back.
  l <- leaf_model()
  do.call(set_drivers, c(list(l), grid_drivers(2.0)))
  l$find_root_collar_psi()
  expect_true(is.finite(l$assim_colimited_))

  set_traits(l, leaf_traits(vcmax_25 = 110))
  expect_false(is.finite(l$assim_colimited_))
})

test_that("set_traits() enforces the one representation for psi (#25)", {
  # A bare field write would bypass this, which is one of the four reasons the
  # traits are not bound as settable fields.
  l <- leaf_model()
  for (p in c("psi_crit", "stem_b", "root_b", "root_psi_crit")) {
    tr <- leaf_traits()
    tr[[p]] <- -tr[[p]]
    expect_error(set_traits(l, structure(tr, class = class(leaf_traits()))),
                 "positive magnitude", label = p)
  }
})

test_that("the composite reproduces the arbitrated reference gradients", {
  # psi_soil = 2, PPFD = 900, VPD = 2, one layer, default traits. H = -8.9561.
  # Ratios against a least-squares slope over +-2% at n = 41 were 0.9979-1.0000
  # when these were established; the tolerance below is that agreement, not the
  # composite's own precision, which is finer.
  ref <- rbind(
    vcmax_25        = c(collar =  1.7459e-03, A =  1.7209e-02),
    jmax_25         = c(collar =  1.0862e-04, A =  9.1320e-04),
    cost_scale_TF24 = c(collar = -7.6704e-02, A = -3.9520e-01),
    beta2           = c(collar = -4.2610e-02, A = -2.1954e-01),
    stem_b          = c(collar =  5.5247e-01, A =  2.8465e+00),
    stem_c          = c(collar = -1.6390e-01, A = -8.4448e-01),
    root_b          = c(collar =  2.6656e-04, A =  2.3895e-02),
    beta_R_H        = c(collar =  8.1973e-06, A = -2.8802e-04))

  g <- grid_gradient(2.0, pars = rownames(ref))
  expect_identical(g$status, "interior")
  expect_identical(g$method, "ift")
  expect_equal(g$H, -8.9561, tolerance = 1e-4)

  expect_equal(g$gradient[rownames(ref), "collar"], ref[, "collar"],
               tolerance = 5e-3)
  expect_equal(g$gradient[rownames(ref), "A"], ref[, "A"], tolerance = 5e-3)
})

test_that("the composite and the finite difference agree at interior points", {
  # The cross-check that needs no external reference: two routes to the same
  # derivative, sharing only the model. They are checked across soil moisture
  # because the indirect term's share of the answer varies with it.
  pars <- c("vcmax_25", "jmax_25", "stem_b", "cost_scale_TF24", "root_b")
  for (psi_soil in c(0.5, 1.0, 2.0, 3.0)) {
    ift <- grid_gradient(psi_soil, pars = pars, method = "ift")
    fd <- grid_gradient(psi_soil, pars = pars, method = "fd")
    expect_identical(ift$status, "interior", label = paste("psi_soil", psi_soil))
    # 1e-3 relative, which is the finite difference's precision rather than the
    # composite's: differencing the solved argmax carries the solve's own ~1e-09
    # floor divided by the step. Do not read agreement here as either route
    # being good to more than three digits.
    expect_equal(ift$gradient, fd$gradient, tolerance = 1e-3,
                 label = paste("psi_soil", psi_soil))
  }
})

test_that("dcollar/dtheta is dpsi*/dtheta, computed two ways", {
  # The `collar` column is psi* itself, so the composite reports -M/H there while
  # the fallback differences the solved argmax. They are different computations of
  # the same quantity, which is why the column is worth reporting at all.
  pars <- c("vcmax_25", "stem_b", "beta2")
  ift <- grid_gradient(2.0, pars = pars, method = "ift")
  fd <- grid_gradient(2.0, pars = pars, method = "fd")
  expect_equal(ift$gradient[, "collar"], fd$gradient[, "collar"],
               tolerance = 1e-3)
})

test_that("a pinned optimum takes the fallback, and the composite would be wrong", {
  # The hazard, pinned to specific operating points. These three are rows of the
  # golden grid; the wet-pinned pair are where the bare composite is worst.
  #
  # ⚠️ The truth at a wet-pinned point is ~1e-07 or smaller -- the leaf is stuck
  # at the wettest feasible collar and barely responds -- so the ASSERTION IS
  # ABSOLUTE, not relative. A ratio test here would be meaningless in both
  # directions: the fallback's relative error is large and does not matter, and
  # the composite's O(1) answer is catastrophic precisely because the truth is
  # not O(1).
  wet <- list(list(psi_soil = 4, vpd = 2.0, layers = 5L),
              list(psi_soil = 3, vpd = 4.0, layers = 3L))
  for (w in wet) {
    lab <- sprintf("psi_soil=%g vpd=%g layers=%d", w$psi_soil, w$vpd, w$layers)
    auto <- grid_gradient(w$psi_soil, vpd = w$vpd, layers = w$layers,
                          pars = c("stem_b", "vcmax_25"))
    expect_identical(auto$status, "pinned", label = lab)
    expect_identical(auto$method, "fd", label = lab)
    # The guard fired with room to spare: measured stationarity at these points
    # is ~6e-06 to ~9e-06 against a 1e-08 threshold.
    expect_gt(auto$stationarity, 1e-7)
    expect_lt(abs(auto$gradient["stem_b", "A"]), 1e-3)

    # Forcing the composite here does not produce the wrong number -- it fails.
    # That was not the design and is worth stating as a measurement: psi* sits
    # 1e-06 of a bracket width from its bound at a pinned point, so the step in
    # psi cannot be centred without clamping, and that is checked. The grid test
    # below shows it holds at all 42 pinned rows, which means the composite's
    # O(1) answer is not reachable through this function.
    expect_error(grid_gradient(w$psi_soil, vpd = w$vpd, layers = w$layers,
                               pars = "stem_b", method = "ift"),
                 "narrower than one step", label = lab)
  }

  # Pinned DRY is the milder case and is worth separating: the composite is only
  # a factor of 1.6 out for vcmax_25 there rather than seven orders, so a test
  # that only looked at the wet rows would let a threshold drift past it.
  dry <- grid_gradient(4, vpd = 0.5, layers = 3L,
                       pars = c("stem_b", "vcmax_25"))
  expect_identical(dry$status, "pinned")
  expect_equal(dry$gradient["vcmax_25", "A"], 0.014357, tolerance = 1e-3)
  expect_equal(dry$gradient["stem_b", "A"], 4.6991, tolerance = 1e-3)
})

test_that("at a pinned optimum the BOUND's own trait carries the gradient", {
  # The sharpest statement of why the fallback is not merely a fallback, and it
  # is a stronger claim than "the composite is inaccurate there".
  #
  # `psi_crit` does not appear in the profit function at all -- it only sets the
  # dry end of the feasible collar interval. So at an INTERIOR optimum it has no
  # effect and its gradient is exactly zero, which is right. At a DRY-PINNED
  # optimum it IS the binding constraint, so moving it moves the answer directly.
  #
  # ⚠️ The composite cannot see that, and not by a small margin: psi* enters
  # `dprofit` nowhere, so M = 0, dpsi*/dtheta = 0, the direct term is 0, and the
  # composite returns EXACTLY ZERO where the truth is 1.26. A silent zero is
  # arguably a worse failure than the wet-pinned 3.5e+07 -- it tells an optimiser
  # the parameter does nothing, and nothing about it looks wrong.
  interior <- grid_gradient(2.0, pars = c("psi_crit", "root_psi_crit"))
  expect_identical(interior$status, "interior")
  expect_equal(interior$gradient["psi_crit", "A"], 0)
  expect_equal(interior$gradient["root_psi_crit", "A"], 0)

  pinned <- grid_gradient(4, vpd = 0.5, layers = 3L, pars = "psi_crit")
  expect_identical(pinned$status, "pinned")
  expect_identical(pinned$method, "fd")
  # Arbitrated against a least-squares slope of the solved A over +-2% of
  # psi_crit at n = 41, which gives 1.2605 (R^2 = 0.998; the fit is imperfect
  # because the response has a kink in it, which is the whole point).
  expect_equal(pinned$gradient["psi_crit", "A"], 1.2605, tolerance = 5e-3)
  expect_gt(pinned$gradient["psi_crit", "collar"], 0)
})

test_that("a shut-down operating point reports no gradient and still differences", {
  # psi_soil = 6 is drier than psi_crit, so there is no optimisation to
  # differentiate: dprofit returns a sentinel zero rather than a derivative, and
  # H comes back zero with it. The composite has nothing to stand on and says so
  # -- but the gradient itself is not zero, because R_d still depends on
  # vcmax_25, so the fallback still has work to do.
  g <- grid_gradient(6.0, pars = c("vcmax_25", "stem_b"))
  expect_identical(g$status, "no-gradient")
  expect_identical(g$method, "fd")
  expect_equal(g$gradient["vcmax_25", "A"], -0.015, tolerance = 1e-4)
  expect_equal(g$gradient["stem_b", "A"], 0)

  # Forcing the composite here is an error rather than a wrong number: unlike a
  # pinned point, there is no curvature to divide by at all.
  expect_error(grid_gradient(6.0, pars = "vcmax_25", method = "ift"),
               "usable curvature")
})

test_that("the active-set classification splits the golden grid as measured", {
  # The whole 288-point grid, which is what makes `stationarity_tol` a measured
  # threshold rather than a guess. Two claims:
  #
  #   * the counts: 198 interior, 42 pinned, 48 with no gradient at all;
  #   * the GAP. The worst interior point is ~5e-11 and the mildest pinned one
  #     ~6e-06, so the default 1e-08 sits in an empty band four orders wide on
  #     each side. That is the number worth defending, because a threshold in a
  #     populated region would be a judgement call instead.
  psi_soils <- c(0.5, 1.0, 2.0, 3.0, 4.0, 6.0)
  vpds <- c(0.5, 1.0, 2.0, 4.0)
  layer_counts <- c(1L, 3L, 5L)
  ppfds <- c(100, 500, 900, 1500)

  rows <- expand.grid(layers = layer_counts, vpd = vpds, ppfd = ppfds,
                      psi_soil = psi_soils)
  res <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    # One trait only: the classification is a property of the operating point,
    # not of which trait is being differentiated.
    g <- grid_gradient(r$psi_soil, ppfd = r$ppfd, vpd = r$vpd,
                       layers = r$layers, pars = "vcmax_25")
    # Does forcing the composite here produce a number or refuse? Recorded per
    # row so the claim below is about the whole grid rather than a sample.
    forced <- tryCatch({
      grid_gradient(r$psi_soil, ppfd = r$ppfd, vpd = r$vpd, layers = r$layers,
                    pars = "vcmax_25", method = "ift")
      TRUE
    }, error = function(e) FALSE)
    data.frame(psi_soil = r$psi_soil, status = g$status,
               stationarity = g$stationarity, forced_ift_ran = forced)
  }))

  expect_identical(nrow(res), 288L)
  expect_identical(sum(res$status == "interior"), 198L)
  expect_identical(sum(res$status == "pinned"), 42L)
  expect_identical(sum(res$status == "no-gradient"), 48L)

  # Pinned points are all at the dry end, which is where a drought calibration
  # wanders and is the reason this is not a curiosity.
  expect_setequal(unique(res$psi_soil[res$status == "pinned"]), c(3, 4))

  worst_interior <- max(res$stationarity[res$status == "interior"])
  mildest_pinned <- min(res$stationarity[res$status == "pinned"])
  expect_lt(worst_interior, 1e-9)
  expect_gt(mildest_pinned, 1e-7)

  # The second guard, over the whole grid: the composite runs on every interior
  # row and refuses on every other one. So a caller who forces `ift` cannot
  # obtain the silently-wrong answer -- they get an error. This is a measured
  # property of where these points sit, not a designed one, and it is asserted
  # here so that it is noticed if it stops being true.
  expect_true(all(res$forced_ift_ran[res$status == "interior"]))
  expect_false(any(res$forced_ift_ran[res$status != "interior"]))
})

test_that("leaf_gradient() works on the single-potential supply path", {
  g <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                     supply = leaf_supply_single(resistance = 1e3),
                     pars = c("vcmax_25", "stem_b"))
  expect_identical(g$status, "interior")
  fd <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                      supply = leaf_supply_single(resistance = 1e3),
                      pars = c("vcmax_25", "stem_b"), method = "fd")
  expect_equal(g$gradient, fd$gradient, tolerance = 1e-3)
})

test_that("leaf_gradient() leaves the traits it was given", {
  # The perturbation loop mutates one Leaf, so it has to put it back. The Leaf is
  # internal, so what this really checks is that two calls agree -- a loop that
  # leaked its last perturbation would drift.
  first <- grid_gradient(2.0, pars = "vcmax_25")
  second <- grid_gradient(2.0, pars = "vcmax_25")
  expect_identical(first$gradient, second$gradient)
  expect_identical(first$value, second$value)
})

test_that("GSS_tol_abs does not reach the production path", {
  # Not a gradient test, but it belongs with the ones above: it is the claim that
  # makes them possible. Since PLAN 11a the collar solve solves dprofit == 0 at
  # 1e-12 instead of searching profit to GSS_tol_abs, so this control has no
  # effect on the answer -- and an argmax determined only to 1e-3 could not be
  # differentiated at all.
  #
  # Asserted because it was documented WRONG in two places until this branch: the
  # vignette demonstrated 1e-3 against 1e-5 and said they "differ in the fourth
  # decimal", and leaf_control()'s roxygen said this sets how well the operating
  # point is determined. Both were true before 11a. Bit-exact, because "no
  # effect" is exactly what is being claimed.
  loose <- leaf_solve(psi_soil = 2.0, PPFD = 900)
  tight <- leaf_solve(psi_soil = 2.0, PPFD = 900,
                      control = leaf_control(GSS_tol_abs = 1e-5))
  expect_identical(loose, tight)
})

test_that("the documented examples classify the way their comments claim", {
  # Written after the first version of ?leaf_gradient offered
  # `psi_soil = c(4, 4.25, 4.5), atm_vpd = 2` as "a pinned optimum" and it came
  # back interior. An example whose comment contradicts its output is worse than
  # no example, and nothing else here would have caught it: examples are run by
  # R CMD check but their *classifications* are not asserted anywhere. So assert
  # them, using exactly the calls in the roxygen block.
  expect_identical(leaf_gradient(psi_soil = 2.0, PPFD = 900)$method, "ift")

  dry <- leaf_gradient(psi_soil = 4.5, PPFD = 900, atm_vpd = 3.0)
  expect_identical(dry$method, "fd")
  expect_identical(dry$status, "pinned")
})

test_that("the two non-trait parameters agree with a resolved reference", {
  # #44. Half of leaf-calibration's free parameters are not traits: `K_total` and
  # `f_plant` reach the leaf as leaf_specific_conductance_max (a driver) and the
  # single-potential resistance. Nothing in the derivation cares -- the implicit
  # function theorem is applied to dprofit/dpsi = 0, and any parameter profit
  # depends on goes through it identically -- so the thing to check is that the
  # plumbing applies the perturbation through the right setter.
  #
  # Arbitrated the way 11c arbitrated the traits: against a central difference of
  # the WHOLE SOLVE, read where the step sweep plateaus rather than at one step.
  d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
            supply = leaf_supply_single(resistance = 1e4))

  solve_at <- function(par, value) {
    a <- d
    if (identical(par, "resistance")) {
      a$supply <- leaf_supply_single(resistance = value)
    } else {
      a[[par]] <- value
    }
    x <- do.call(leaf_solve, a)
    c(A = x$A, gc = x$gc, psi_stem = x$psi_stem, collar = x$collar)
  }

  g <- do.call(leaf_gradient,
               c(d, list(pars = c("leaf_specific_conductance_max",
                                  "resistance"), method = "ift")))
  expect_identical(g$status, "interior")

  for (p in c("leaf_specific_conductance_max", "resistance")) {
    v <- if (identical(p, "resistance")) 1e4 else 3.14e-5
    h <- v * 1e-5
    ref <- (solve_at(p, v + h) - solve_at(p, v - h)) / (2 * h)
    # Six digits, which is well inside the reference's own precision. Measured
    # agreement is 8 digits for the driver and 9 for the resistance.
    expect_equal(g$gradient[p, ], ref, tolerance = 1e-6)
    # And not accidentally zero, which is how this would fail if the setter did
    # nothing at all -- the failure mode psi_crit genuinely has (see above).
    expect_gt(abs(g$gradient[p, "A"]), 1e-8)
  }
})

test_that("the non-trait gradients are stable across decades of step", {
  # The mixed partial is a difference of an exact quantity, so it should be
  # step-insensitive over a wide band; a parameter whose step rule was wrong
  # would show up here as drift rather than as a wrong answer.
  #
  # ⚠️ This is the test that would have caught the step-rule trap:
  # leaf_specific_conductance_max defaults to 3.14e-05, so the traits' rule of
  # flooring the step at 1 would perturb it by 3% and measure a secant.
  d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
            supply = leaf_supply_single(resistance = 1e4))
  for (p in c("leaf_specific_conductance_max", "resistance")) {
    vals <- vapply(10^-(4:7), function(st) {
      do.call(leaf_gradient, c(d, list(pars = p, step = st,
                                       method = "ift")))$gradient[p, "A"]
    }, numeric(1))
    expect_lt(diff(range(vals)) / abs(mean(vals)), 1e-6)
  }
})

test_that("the active-set guard covers the non-trait parameters too", {
  # Both move the FEASIBLE COLLAR INTERVAL and not only the profit function --
  # resistance directly, the conductance through the supply curve -- so a
  # perturbation can push a bound past psi* more readily than a trait
  # perturbation can. The guard must not be trait-specific.
  dry <- list(psi_soil = 5.5, PPFD = 900, atm_vpd = 3.0,
              supply = leaf_supply_single(resistance = 1e4))
  auto <- do.call(leaf_gradient,
                  c(dry, list(pars = c("leaf_specific_conductance_max",
                                       "resistance"))))
  expect_true(auto$status %in% c("pinned", "no-gradient"))
  expect_identical(auto$method, "fd")
  expect_true(all(is.finite(auto$gradient)))

  # ...and forcing the composite there fails loudly rather than returning a
  # confident number, which is the whole reason `auto` exists.
  expect_error(do.call(leaf_gradient,
                       c(dry, list(pars = "resistance", method = "ift"))))
})

test_that("leaf_gradient() rejects bad arguments", {
  expect_error(grid_gradient(2.0, pars = "not_a_trait"),
               "cannot differentiate")
  # `resistance` is a real parameter -- on the other supply path -- so the
  # message says which path it belongs to rather than calling it unknown.
  expect_error(grid_gradient(2.0, pars = "resistance"),
               "single-potential path only")
  expect_error(grid_gradient(2.0, step = -1), "positive number")
  expect_error(grid_gradient(2.0, traits = list(vcmax_25 = 96)),
               "leaf_traits")
  expect_error(grid_gradient(2.0, method = "magic"))
})

test_that("the gradient is reported for every output the fit needs", {
  # leaf-calibration fits three responses -- A, gs and psi_leaf -- so all three
  # are differentiated, not just A. `collar` comes along because it is psi*.
  g <- grid_gradient(2.0, pars = "vcmax_25")
  expect_identical(colnames(g$gradient), c("A", "gc", "psi_stem", "collar"))
  expect_identical(names(g$value), c("A", "gc", "psi_stem", "collar"))
  expect_true(all(is.finite(g$gradient)))
  # Raising vcmax_25 raises assimilation and opens the stomata, and the leaf pays
  # for it with a more negative water potential (a larger positive magnitude).
  expect_gt(g$gradient["vcmax_25", "A"], 0)
  expect_gt(g$gradient["vcmax_25", "gc"], 0)
  expect_gt(g$gradient["vcmax_25", "psi_stem"], 0)
})
