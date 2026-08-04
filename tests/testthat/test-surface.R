# The stage 2 R surface: leaf_traits / leaf_control / leaf_model / set_drivers /
# operating_point / leaf_solve.
#
# The load-bearing tests here are the ones that check the convenient path and
# the raw path agree. Everything in R/leaf-model.R is a convenience over the
# generated glue, and the failure mode of a convenience layer is not that it
# errors -- it is that it quietly reorders an argument and returns a plausible
# number. test-golden.R pins the raw path against the C++ baseline; these tests
# pin the friendly path against the raw one, so the chain reaches the golden
# file from both.

test_that("leaf_traits() and leaf_control() partition the C++ constructor", {
  # Between them they must cover the 19 constructor arguments exactly once --
  # no argument silently dropped (which would leave it at whatever the
  # positional call happens to supply) and none duplicated.
  ctor_args <- setdiff(names(formals(Leaf)), "")
  covered <- c(names(leaf_traits()), names(leaf_control()))

  expect_setequal(setdiff(ctor_args, covered), character(0))
  expect_setequal(setdiff(covered, ctor_args),
                  c("integration_rule", "integration_tol"))
  expect_length(intersect(names(leaf_traits()), names(leaf_control())), 0)

  # And the split is the one the issue asked for: tolerances on the control
  # side, physiology on the trait side.
  expect_true(all(c("GSS_tol_abs", "ci_abs_tol", "ci_niter",
                    "vulnerability_curve_ncontrol") %in% names(leaf_control())))
  expect_false(any(grepl("tol|niter", names(leaf_traits()))))
})

test_that("leaf_model() and the raw Leaf() constructor agree", {
  # The reason leaf_model() exists is that mapping 15 traits and 4 tolerances
  # onto 19 positional slots is exactly the kind of thing that goes wrong once
  # and is never noticed. So check it against a hand-written positional call
  # with the same values, on a full solve rather than on the arguments.
  raw <- Leaf(96, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245, 5.870283,
              1.5, 157.44, 0.30, 0.7, 0.99, 1e-3, 100, 1e-3, 1000, 7.5,
              3.4e2, 9.4e3)
  raw$initialize_integrator(21, 1e-8)
  friendly <- leaf_model()

  for (l in list(raw, friendly)) {
    l$set_physiology(20, 900, 2.0, 1.0, golden_kmax, 2.0, 40.0, 25, 21, 101.3)
    l$find_root_collar_psi()
  }
  expect_identical(operating_point(raw), operating_point(friendly))
})

test_that("a non-default trait reaches the model through leaf_model()", {
  # The previous test would pass even if leaf_model() ignored `traits` entirely
  # and always used the defaults, because the defaults are what it compares. So
  # move one and check it lands in the right slot -- vcmax_25, which raises
  # assimilation, against stem_b, which moves the vulnerability curve.
  base <- leaf_solve(psi_soil = 2.0, PPFD = 900)
  hi_vcmax <- leaf_solve(psi_soil = 2.0, PPFD = 900,
                         traits = leaf_traits(vcmax_25 = 150))
  expect_gt(hi_vcmax$A, base$A)

  brittle <- leaf_model(leaf_traits(stem_b = 2.0))
  expect_lt(brittle$proportion_of_conductivity(2.0),
            leaf_model()$proportion_of_conductivity(2.0))
})

test_that("a control setting reaches the model and is not treated as a trait", {
  # vulnerability_curve_ncontrol changes the spline resolution, so it moves the
  # answer slightly without being a trait. That it moves the answer at all is
  # what shows it reached the constructor rather than being dropped.
  coarse <- leaf_solve(psi_soil = 2.0, PPFD = 900,
                       control = leaf_control(vulnerability_curve_ncontrol = 20))
  fine <- leaf_solve(psi_soil = 2.0, PPFD = 900)
  expect_false(identical(coarse$A, fine$A))
  # Still the same model, though -- a resolution change is not a trait change.
  expect_equal(coarse$A, fine$A, tolerance = 1e-2)
})

test_that("set_drivers() agrees with a positional set_physiology()", {
  raw <- leaf_model()
  raw$set_physiology(rep(20 / 3, 3), 900, c(1, 1.25, 1.5), c(1, 2, 3),
                     golden_kmax, 2.0, 40.0, 25, 21, 101.3)
  raw$find_root_collar_psi()

  friendly <- leaf_model()
  set_drivers(friendly, psi_soil = c(1, 1.25, 1.5), PPFD = 900,
              leaf_specific_conductance_max = golden_kmax)
  friendly$find_root_collar_psi()

  # The defaults set_drivers() supplies -- 1 m layers, root carbon split evenly
  # -- must be exactly what was passed by hand above.
  expect_identical(operating_point(raw), operating_point(friendly))
})

test_that("set_drivers() keeps the signed-psi_soil rejection", {
  # #25's input-boundary assertion is the only thing between a pre-#25 script and
  # a plausible wrong number. A friendly wrapper is precisely where someone would
  # be tempted to "helpfully" abs() it.
  l <- leaf_model()
  expect_error(set_drivers(l, psi_soil = -2.0), "positive magnitudes in MPa")
  expect_error(set_drivers(l, psi_soil = c(1, -2)), "positive magnitudes in MPa")
  expect_error(leaf_solve(psi_soil = -2.0), "positive magnitudes in MPa")
})

test_that("set_drivers() rejects mismatched layer counts", {
  l <- leaf_model()
  expect_error(set_drivers(l, psi_soil = c(1, 2), soil_depth = 1),
               "one entry per soil layer")
  expect_error(
    set_drivers(l, psi_soil = c(1, 2), root_carbon_per_leaf_area = c(1, 2, 3)),
    "one entry per soil layer")
  expect_error(set_drivers(l, psi_soil = numeric(0)), "at least one layer")
})

test_that("leaf_traits() and leaf_control() reject non-scalars", {
  expect_error(leaf_traits(vcmax_25 = c(90, 100)), "single finite number")
  expect_error(leaf_traits(vcmax_25 = NA_real_), "single finite number")
  expect_error(leaf_control(ci_niter = NULL), "single finite number")
  expect_error(leaf_model(traits = list()), "must come from leaf_traits")
  expect_error(leaf_model(control = list()), "must come from leaf_control")
})

test_that("leaf_solve() reproduces the stateful path exactly", {
  stateful <- leaf_model()
  set_drivers(stateful, psi_soil = 2.0, PPFD = 900,
              leaf_specific_conductance_max = golden_kmax)
  stateful$find_root_collar_psi()

  one_call <- leaf_solve(psi_soil = 2.0, PPFD = 900,
                         leaf_specific_conductance_max = golden_kmax)
  expect_identical(one_call[, names(operating_point(stateful))],
                   operating_point(stateful))
})

test_that("leaf_solve() vectorises and recycles", {
  d <- leaf_solve(psi_soil = c(0.5, 1, 2), PPFD = 900)
  expect_s3_class(d, "data.frame")
  expect_identical(nrow(d), 3L)
  expect_identical(d$psi_soil, c(0.5, 1, 2))
  expect_identical(unique(d$PPFD), 900)

  # Two vectors of the same length pair up elementwise; a scalar recycles.
  pairs <- leaf_solve(psi_soil = c(1, 2), PPFD = c(400, 1200))
  expect_identical(nrow(pairs), 2L)
  expect_identical(pairs$PPFD, c(400, 1200))

  expect_error(leaf_solve(psi_soil = c(1, 2, 3), PPFD = c(400, 1200)),
               "cannot be recycled")
})

test_that("leaf_solve() takes a list for multi-layer profiles", {
  d <- leaf_solve(psi_soil = list(c(1, 1.5, 2), c(3, 3.5, 4)), PPFD = 900)
  expect_identical(nrow(d), 2L)
  expect_identical(d$layers, c(3L, 3L))
  # psi_soil is reported as the wettest layer, so the column stays numeric.
  expect_identical(d$psi_soil, c(1, 3))
  # Drier profile, less assimilation.
  expect_lt(d$A[[2]], d$A[[1]])
})

test_that("reuse = TRUE and reuse = FALSE give identical answers", {
  # This is hazard 8 checked at the R level, and it is the reason leaf_solve()
  # is allowed to reuse one object for speed. Leaf is a value type that plant
  # reuses for every individual in a patch, and three code paths have historically
  # left an output holding the PREVIOUS solve's value. The drivers below cross
  # from a wet, productive point to a shut-down one and back, which is the
  # transition that exposed those bugs.
  drivers <- list(psi_soil = c(0.5, 6.0, 1.0, 6.0, 0.5),
                  PPFD = c(1500, 100, 900, 1500, 100),
                  atm_vpd = c(0.5, 4.0, 2.0, 4.0, 0.5))
  shared <- do.call(leaf_solve, c(drivers, list(reuse = TRUE)))
  fresh <- do.call(leaf_solve, c(drivers, list(reuse = FALSE)))
  expect_identical(shared, fresh)

  # And the sweep really does reach the shut-down branch, so the check above is
  # not vacuous -- reading a bit-identical result as evidence without checking
  # the branch was reached is a mistake this project has made before.
  expect_true(any(shared$E == 0))
  expect_true(any(shared$E > 0))
})

test_that("operating_point() reports lambda and g1_eff from the solved state", {
  d <- leaf_solve(psi_soil = 1.0, PPFD = 900)
  expect_true(is.finite(d$lambda))
  expect_true(is.finite(d$g1_eff))
  expect_gt(d$lambda, 0)

  # g1_eff is computed from the solve rather than supplied, so it must vary with
  # conditions. A flat curve here would mean it had become disconnected from the
  # solved state -- reported, but not actually reporting anything.
  sweep <- leaf_solve(psi_soil = seq(0.5, 4, length.out = 8), PPFD = 900)
  live <- sweep[sweep$gc > 0 & is.finite(sweep$g1_eff), ]
  expect_gt(nrow(live), 2)
  expect_gt(diff(range(live$g1_eff)) / mean(live$g1_eff), 0.1)
})

test_that("atm_kpa is not decorative", {
  # 10c: the ppm-to-Pa conversion is derived from atm_kpa. This was hard-coded at
  # 101.3 until #15, and plant's own driver default is 100.5 -- which turned out
  # to move TF24 offspring production by +2.4%, against +0.10% for the entire
  # rest of the swap. So the R surface must actually pass it through, and a test
  # that only ever uses 101.3 would not notice if it stopped.
  a <- leaf_solve(psi_soil = 2.0, PPFD = 900, atm_kpa = 101.3)
  b <- leaf_solve(psi_soil = 2.0, PPFD = 900, atm_kpa = 100.5)
  expect_false(identical(a$A, b$A))
})
