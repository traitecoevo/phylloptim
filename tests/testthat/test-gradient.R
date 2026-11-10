# Trait gradients: set_traits() and leaf_gradient() (#4).
#
# The gradient is one pass of the derivative surface at the solved point, which
# places the decision variable from the kind the solve recorded: by the implicit
# function theorem at an interior optimum, on the bound that holds it at a
# pinned one. Two things are pinned here:
#
#   1. That the rows are the model's: against arbitrated references, and
#      against a central difference of the whole solve (solve_difference(), in
#      helper-gradient.R), at interior AND pinned points.
#   2. That the status says which of those the point is. At a pinned point the
#      implicit function theorem's answer would be O(1) where the truth is
#      ~1e-08, so the status is what tells a caller the rows came off a bound.

# `grid_drivers()` and `grid_gradient()` moved to helper-gradient.R when
# test-gradient-batch.R needed the same operating points: each test file gets its
# own environment, so a second copy would have been free to drift and the two
# files would have pinned different points while appearing to pin the same ones.

# The production path these tests were written on: TF24 on the collar route
# over the multi-layer supply. A call that names its own model, route or supply
# keeps it.
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
# set_drivers() takes the drivers, not the supply the leaf was built with.
grid_set_drivers <- function(l, ...) {
  a <- list(...)
  a$supply <- NULL
  do.call(set_drivers, c(list(l), a))
}

test_that("set_traits() on a used leaf equals a leaf built with those traits", {
  # The R-side statement of the C++ suite's test_set_traits_matches_a_fresh_leaf,
  # made at the boundary a calibration loop actually crosses. Bit-exact, because
  # the two routes share no code: any derived state set_traits fails to refresh
  # -- the two vulnerability splines, or vcmax_/jmax_/R_d_ behind
  # set_physiology's temperature cache -- shows up as a difference.
  traits <- leaf_traits(vcmax_25 = 110, stem_P50 = 4.2, root_P50 = 3.5,
                        TF24_cost_scale = 8.0)
  d <- grid_drivers(2.0)

  fresh <- leaf_model(traits, supply = leaf_supply_multilayer())
  do.call(grid_set_drivers, c(list(fresh), d))
  fresh$find_root_collar_psi()

  # Solved at the DEFAULTS first, so every cache is warm and pointing at the old
  # traits. Solving first is what makes this a test rather than a coincidence.
  reused <- leaf_model(supply = leaf_supply_multilayer())
  do.call(grid_set_drivers, c(list(reused), d))
  reused$find_root_collar_psi()
  set_traits(reused, traits)
  do.call(grid_set_drivers, c(list(reused), d))
  reused$find_root_collar_psi()

  expect_identical(operating_point(reused), operating_point(fresh))
  expect_identical(reused$vcmax_, fresh$vcmax_)
  expect_identical(reused$R_d_, fresh$R_d_)
})

test_that("set_traits() requires the drivers to be set again", {
  # Not a nicety: the derived photosynthetic parameters really are unknown until
  # set_physiology runs, so the operating point must not be readable as if it
  # were still valid. setup_clean_leaf puts the NA sentinels back.
  l <- leaf_model(supply = leaf_supply_multilayer())
  do.call(grid_set_drivers, c(list(l), grid_drivers(2.0)))
  l$find_root_collar_psi()
  expect_true(is.finite(l$assim_colimited_))

  set_traits(l, leaf_traits(vcmax_25 = 110))
  expect_false(is.finite(l$assim_colimited_))
})

test_that("set_traits() enforces the one representation for psi (#25)", {
  # A bare field write would bypass this, which is one of the four reasons the
  # traits are not bound as settable fields.
  l <- leaf_model(supply = leaf_supply_multilayer())
  for (p in c("stem_P50", "root_P50")) {
    tr <- leaf_traits()
    tr[[p]] <- -tr[[p]]
    expect_error(set_traits(l, structure(tr, class = class(leaf_traits()))),
                 "positive magnitude", label = p)
  }
})

test_that("the gradient reproduces the arbitrated reference gradients", {
  # psi_soil = 2, PPFD = 900, VPD = 2, one layer, default traits.
  # Ratios against a least-squares slope over +-2% at n = 41 were 0.9979-1.0000
  # when these were established; the tolerance below is that agreement, not the
  # gradient's own precision, which is finer.
  # ⚠️ THE CURVE ROWS ARE NOT RE-ARBITRATED NUMBERS -- they are the arbitrated
  # b-gradients carried through the reparameterisation by its own chain rule, so
  # this table asserts the chain rule as well as the gradient. With
  # b = P50 / (ln 2)^(1/c),
  #
  #     d/dP50 |_c = (ln 2)^(-1/c) . d/db
  #     d/dc |_P50 = d/dc |_b  +  (db/dc) . d/db,   db/dc = b.ln(ln 2)/c^2
  #
  # and the derived psi_crit contributes NOTHING at this row because it is an
  # INTERIOR optimum, where dY/dpsi_crit is exactly zero. That is why the identity
  # is clean here and is not at a pinned optimum -- see the pinned test below,
  # which asserts the other half.
  tr <- leaf_traits()
  stem_b <- tr$stem_P50 / log(2)^(1 / tr$stem_c)
  d_dP50 <- log(2)^(-1 / tr$stem_c)
  d_dc   <- stem_b * log(log(2)) / tr$stem_c^2
  b_stem <- c(collar = 5.5247e-01, A = 2.8465e+00)   # arbitrated d/dstem_b
  b_root <- c(collar = 2.6656e-04, A = 2.3895e-02)   # arbitrated d/droot_b
  ref <- rbind(
    vcmax_25        = c(collar =  1.7459e-03, A =  1.7209e-02),
    jmax_25         = c(collar =  1.0862e-04, A =  9.1320e-04),
    TF24_cost_scale = c(collar = -7.6704e-02, A = -3.9520e-01),
    TF24_beta2           = c(collar = -4.2610e-02, A = -2.1954e-01),
    stem_P50        = d_dP50 * b_stem,
    stem_c          = c(collar = -1.6390e-01, A = -8.4448e-01) + d_dc * b_stem,
    root_P50        = (log(2)^(-1 / tr$root_c)) * b_root)
  # beta_R_H used to be an eighth row here (collar 8.1973e-06, A -2.8802e-04).
  # #33 removed it from the trait vector along with the root-architecture model,
  # so there is no longer a `pars` name for it. The recorded values are kept in
  # this comment rather than deleted, because they are what a hand-rolled
  # two-network difference should reproduce if anyone needs to check that route.

  g <- grid_gradient(2.0, pars = rownames(ref))
  expect_identical(g$status, "interior")
  # ⚠️ THE MODEL'S OWN NUMBER, NOT AN ARBITRATED ONE, and it had drifted from the
  # model by more than this tolerance: the literal read -8.9578 and the comment
  # above it -8.9561, where the model gives -8.9554 at 100 knots. It is pinned at
  # the CONVERGED value, which holds at every resolution measured -- -8.955414,
  # -8.955057, -8.955110, -8.955133, -8.955141 at 100, 200, 400, 800 and 1600
  # knots, a spread of 4e-05 relative and so well inside this tolerance. The
  # seven rows below are the arbitrated ones and did not move.
  expect_equal(g$H, -8.9551, tolerance = 1e-4)

  expect_equal(g$gradient[rownames(ref), "collar"], ref[, "collar"],
               tolerance = 5e-3)
  expect_equal(g$gradient[rownames(ref), "A"], ref[, "A"], tolerance = 5e-3)
})

test_that("the gradient agrees with a difference of the solve at interior points", {
  # The cross-check that needs no external reference: the surface and the solve
  # share only the model. Checked across soil moisture because the indirect
  # term's share of the answer varies with it. The `collar` column is psi*
  # itself, so it is the implicit function theorem's -M/H against the moved
  # argmax: the same quantity computed two ways.
  #
  # 1e-6 is the difference's precision, not the gradient's: differencing the
  # solved argmax carries the solve's ~1e-09 floor divided by the step.
  # Measured agreement is 1e-10 to 6e-09 (mean relative).
  pars <- c("vcmax_25", "jmax_25", "stem_P50", "TF24_cost_scale", "root_P50",
            "TF24_beta2")
  for (psi_soil in c(0.5, 1.0, 2.0, 3.0)) {
    g <- grid_gradient(psi_soil, pars = pars)
    expect_identical(g$status, "interior", label = paste("psi_soil", psi_soil))
    expect_equal(g$gradient, solve_difference(grid_solve_args(psi_soil), pars),
                 tolerance = 1e-6, label = paste("psi_soil", psi_soil))
    expect_equal(g$gradient[, "collar"], -g$M / g$H, tolerance = 1e-12)
  }
})

test_that("a pinned optimum's rows are the bound's, and agree with the solve", {
  # The hazard, pinned to specific operating points. These three are rows of the
  # golden grid; the boundary-soil pair are where the implicit function
  # theorem's answer would be worst, O(1) against a truth of ~1e-07.
  #
  # ⚠️ The truth at a boundary-soil point is ~1e-07 or smaller -- the leaf is stuck
  # at the wettest feasible collar and barely responds -- so the ASSERTION IS
  # ABSOLUTE, not relative: a ratio against a near-zero row measures the
  # difference's noise. Measured agreement is 5e-10 or better on every row.
  wet <- list(list(psi_soil = 4, vpd = 2.0, layers = 5L),
              list(psi_soil = 3, vpd = 3.0, layers = 5L))
  pars <- c("stem_P50", "vcmax_25")
  for (w in wet) {
    lab <- sprintf("psi_soil=%g vpd=%g layers=%d", w$psi_soil, w$vpd, w$layers)
    g <- grid_gradient(w$psi_soil, vpd = w$vpd, layers = w$layers, pars = pars)
    expect_identical(g$status, "pinned", label = lab)
    # No curvature is taken off an interior point, so there is none to report.
    expect_true(is.na(g$H), label = lab)
    ref <- solve_difference(grid_solve_args(w$psi_soil, vpd = w$vpd,
                                            layers = w$layers), pars)
    expect_lt(max(abs(g$gradient - ref)), 1e-8, label = lab)
  }

  # Pinned DRY is the milder case and is worth separating: the implicit
  # function theorem would be only a factor of 1.6 out for vcmax_25 there
  # rather than seven orders, so a test that only looked at the wet rows would
  # let a misclassification past it.
  # The dry pin (#183): a dry soil in humid air, where the leaf transpires all
  # its stem can carry and the stem is held at the end of its curve,
  # b x_end^(1/c). Its water is the stem's capacity, which no photosynthetic
  # trait moves, so the conductance and collar rows in vcmax_25 and R_d_25 are
  # exactly zero while A's are not; and the stem moves with stem_P50 exactly as
  # the end of the curve does, (40 / ln 2)^(1/c) -- the term through the bound,
  # which the surface carries by placing the stem on the end.
  dry_pars <- c("stem_P50", "vcmax_25", "R_d_25")
  dry <- grid_gradient(6, ppfd = 500, vpd = 0.5, layers = 1L, pars = dry_pars)
  expect_identical(dry$status, "pinned")
  for (p in c("vcmax_25", "R_d_25")) {
    expect_identical(dry$gradient[p, "gc"], 0, label = p)
    expect_identical(dry$gradient[p, "collar"], 0, label = p)
    expect_true(dry$gradient[p, "A"] != 0, label = p)
  }
  tr <- leaf_traits()
  expect_equal(dry$gradient["stem_P50", "psi_leaf"],
               (40 / log(2))^(1 / tr$stem_c), tolerance = 1e-6)
  ref <- solve_difference(grid_solve_args(6, ppfd = 500, vpd = 0.5), dry_pars)
  expect_lt(max(abs(dry$gradient - ref)), 1e-8)
})

test_that("at a pinned optimum the moving BOUND carries part of the gradient", {
  # The sharpest statement of why the fallback is not merely a fallback, and it
  # is a stronger claim than "the implicit function theorem is inaccurate there".
  #
  # `psi_crit` sets the dry end of the feasible collar interval and appears in the
  # profit function nowhere else. It is no longer a trait -- it is the
  # 5%-conductivity quantile of the stem curve -- so moving `stem_P50` moves it,
  # and the effect that used to be a column of its own is now a TERM inside the
  # stem gradient. Where the optimum is interior that term is exactly zero; where
  # it is pinned at the dry bound the term is what the answer is made of.
  #
  # ⚠️ The implicit function theorem cannot see it, and not by a small margin:
  # psi* enters `dprofit` nowhere, so M = 0, dpsi*/dtheta = 0 and the direct
  # term through the bound is 0. Placed by the first-order condition here, the
  # pinned gradient would lose the whole bound term, and nothing about the
  # number would look wrong. The surface places a dry pin on the end itself.
  tr <- leaf_traits()
  d_dP50 <- log(2)^(-1 / tr$stem_c)
  bound  <- (log(1 / 0.05) / log(2))^(1 / tr$stem_c)

  # Interior: the bound term drops out, so the gradient is the b-gradient scaled.
  interior <- grid_gradient(2.0, pars = "stem_P50")
  expect_identical(interior$status, "interior")
  expect_equal(interior$gradient["stem_P50", "A"], d_dP50 * 2.8465,
               tolerance = 5e-3)

  # Pinned: the bound term is present and is 40% of the total.
  # At a dry pin the stem is held at the end of its curve, so its row in
  # stem_P50 is the end's own, (40 / ln 2)^(1/c): carried by the bound, not by
  # the curve (#183).
  pinned <- grid_gradient(6, ppfd = 500, vpd = 0.5, layers = 1L, pars = "stem_P50")
  expect_identical(pinned$status, "pinned")
  expect_equal(pinned$gradient["stem_P50", "psi_leaf"],
               (40 / log(2))^(1 / tr$stem_c), tolerance = 1e-6)
  expect_false(isTRUE(all.equal(pinned$gradient["stem_P50", "psi_leaf"],
                                interior$gradient["stem_P50", "psi_leaf"],
                                tolerance = 1e-3)))
  expect_gt(pinned$gradient["stem_P50", "collar"], 0)
})

test_that("a shut-down operating point is pinned at its wet bound", {
  # Past the end of the stem curve the leaf shuts (#183): there is no interval
  # to optimise over, and it sits at its wet bound with its stomata shut. With
  # one layer and g_min = 0 that bound is the zero-uptake collar, the soil
  # potential itself, so the surface places it there and its rows are the
  # bound's. No curvature is taken: there is no interior to take it in.
  g <- grid_gradient(16.0, pars = c("vcmax_25", "stem_P50", "root_P50", "R_d_25"))
  expect_identical(g$status, "pinned")
  expect_true(is.na(g$H))

  # A shut leaf with no conductance nets the root of A(ci) = 0 in light: the
  # g_min -> 0 limit of a leaf whose cuticle still passes CO2, and the limit the
  # transpiring branch reaches as E -> 0. That is zero to the ci solve's
  # tolerance whatever the photosynthetic traits are, so their rows are rounding.
  expect_lt(abs(g$value[["A"]]), 1e-12)
  expect_lt(abs(g$gradient["vcmax_25", "A"]), 1e-9)
  expect_lt(abs(g$gradient["R_d_25", "A"]), 1e-9)
  expect_lt(abs(g$gradient["R_d_25", "profit"]), 1e-9)

  l <- leaf_model(supply = leaf_supply_multilayer())
  do.call(grid_set_drivers, c(list(l), grid_drivers(16.0)))
  l$find_root_collar_psi()
  expect_identical(l$transpiration_, 0)
  expect_identical(g$value[["profit"]],
                   l$assim_colimited_ - l$hydraulic_cost_)

  # The shut stem sits AT the zero-uptake collar -- the soil potential plus the
  # layer's gravity head -- rather than being held at a critical potential (#183
  # decision 1). So neither curve's P50 moves it, and the stem-curve P50 reaches
  # profit only through the cost of a stem at that potential, which falls as the
  # stem gets harder to embolise.
  expect_identical(l$opt_psi_leaf_, l$opt_root_psi_)
  expect_lt(abs(l$opt_psi_leaf_ - 16.0), 0.01)
  expect_identical(g$gradient["stem_P50", "psi_leaf"], 0)
  expect_identical(g$gradient["root_P50", "psi_leaf"], 0)
  expect_identical(g$gradient["root_P50", "profit"], 0)
  # And past the end of the curve the hydraulic cost has saturated, so not even
  # profit sees the stem trait, to rounding: a re-solve agrees.
  expect_lt(abs(g$gradient["stem_P50", "profit"]), 1e-12)
  P50s <- leaf_traits()$stem_P50 + c(-0.05, 0, 0.05)
  profits <- vapply(P50s, function(P50) {
    ll <- leaf_model(leaf_traits(stem_P50 = P50), supply = leaf_supply_multilayer())
    do.call(grid_set_drivers, c(list(ll), grid_drivers(16.0)))
    ll$find_root_collar_psi()
    ll$profit_
  }, numeric(1))
  expect_lt(max(abs(diff(profits))), 1e-12)
})

test_that("the active-set classification splits the golden grid as measured", {
  # The whole 288-point grid. The status is the kind the solve recorded, read
  # rather than tested for: interior where the first-order condition was
  # root-found inside the feasible interval, pinned where a bound holds the
  # collar. 237 interior and 51 pinned, and no point without a gradient.
  #
  # psi_soil = 6 is past psi_crit, which bounds nothing under TF24 (#183), so
  # none of its 48 rows is shut: some are interior optima, and the rest are
  # pinned, at the wet bound where the leaf declines to open or at the dry end
  # where it transpires all its stem can carry.
  psi_soils <- c(0.5, 1.0, 2.0, 3.0, 4.0, 6.0)
  vpds <- c(0.5, 1.0, 2.0, 3.0)  # 3, not 4: below saturation at 25 C
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
    data.frame(psi_soil = r$psi_soil, status = g$status, H = g$H,
               finite = all(is.finite(g$gradient)))
  }))

  expect_identical(nrow(res), 288L)
  expect_identical(sum(res$status == "interior"), 237L)
  expect_identical(sum(res$status == "pinned"), 51L)
  expect_identical(sum(res$status == "no-gradient"), 0L)

  # Open pinned points are all at the dry end, which is where a drought
  # calibration wanders and is the reason this is not a curiosity.
  expect_setequal(unique(res$psi_soil[res$status == "pinned"]), c(3, 4, 6))

  # Every row has a gradient; the curvature is reported where the implicit
  # function theorem was used and nowhere else.
  expect_true(all(res$finite))
  expect_true(all(res$H[res$status == "interior"] < 0))
  expect_true(all(is.na(res$H[res$status == "pinned"])))
})

test_that("leaf_gradient() works on the single-potential supply path", {
  a <- list(psi_soil = 1.5, PPFD = 900, supply = leaf_supply_singlelayer(),
            root_network = series_resistance(1e3), model = "TF24",
            route = "collar")
  pars <- c("vcmax_25", "stem_P50")
  g <- do.call(leaf_gradient, c(a, list(pars = pars)))
  expect_identical(g$status, "interior")
  expect_equal(g$gradient, solve_difference(a, pars), tolerance = 1e-6)
})

test_that("a network rooted above the bottom of the soil profile is differentiated", {
  # A network is sized to its deepest rooted layer, so it can be shorter than
  # the soil profile it is given. The pass builds the soil and the network at
  # their own lengths; indexing one by the other read past the network's end.
  net <- root_network_from_carbon(c(20, 0, 0), soil_depth = c(0.5, 1, 1.5))
  expect_length(net$r_R_H_min, 1L)
  a <- list(psi_soil = c(1, 2, 3), soil_depth = c(0.5, 1, 1.5), PPFD = 900,
            root_network = net, supply = leaf_supply_multilayer(),
            model = "TF24", route = "collar")
  pars <- c("vcmax_25", "root_P50")
  g <- do.call(leaf_gradient, c(a, list(pars = pars)))
  expect_identical(g$status, "interior")
  b <- a
  b$psi_soil <- list(a$psi_soil)
  b$soil_depth <- list(a$soil_depth)
  expect_equal(g$gradient, solve_difference(b, pars), tolerance = 1e-6)
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
  # makes them possible. The collar solve solves dprofit == 0 at
  # 1e-12 instead of searching profit to GSS_tol_abs, so this control has no
  # effect on the answer -- and an argmax determined only to 1e-3 could not be
  # differentiated at all.
  #
  # Asserted because it was documented WRONG in two places until this branch: the
  # vignette demonstrated 1e-3 against 1e-5 and said they "differ in the fourth
  # decimal", and leaf_control()'s roxygen said this sets how well the operating
  # point is determined. Both were true before 11a. Bit-exact, because "no
  # effect" is exactly what is being claimed.
  loose <- leaf_solve(psi_soil = 2.0, PPFD = 900, supply = leaf_supply_multilayer(), model = "TF24", route = "collar")
  tight <- leaf_solve(psi_soil = 2.0, PPFD = 900,
                      control = leaf_control(GSS_tol_abs = 1e-5), supply = leaf_supply_multilayer(), model = "TF24", route = "collar")
  expect_identical(loose, tight)
})

test_that("the documented examples classify the way their comments claim", {
  # Written after the first version of ?leaf_gradient offered
  # `psi_soil = c(4, 4.25, 4.5), atm_vpd = 2` as "a pinned optimum" and it came
  # back interior. An example whose comment contradicts its output is worse than
  # no example, and nothing else here would have caught it: examples are run by
  # R CMD check but their *classifications* are not asserted anywhere. So assert
  # them, using exactly the calls in the roxygen block.
  expect_identical(leaf_gradient(psi_soil = 2.0, PPFD = 900)$status, "interior")

  dry <- leaf_gradient(psi_soil = 4.5, PPFD = 900, atm_vpd = 3.0,
                       supply = leaf_supply_multilayer(), route = "collar")
  expect_identical(dry$status, "pinned")
})

test_that("the two non-trait parameters agree with a resolved reference", {
  # #44. Half of a companion study's free parameters are not traits: `K_total` and
  # `f_plant` reach the leaf as leaf_specific_conductance_max (a driver) and the
  # single-potential resistance. Nothing in the derivation cares -- the implicit
  # function theorem is applied to dprofit/dpsi = 0, and any parameter profit
  # depends on goes through it identically -- so the thing to check is that the
  # plumbing applies the perturbation through the right setter.
  #
  # Arbitrated the way 11c arbitrated the traits: against a central difference of
  # the WHOLE SOLVE, read where the step sweep plateaus rather than at one step.
  d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0,
            supply = leaf_supply_singlelayer(),
            root_network = series_resistance(1e4))

  solve_at <- function(par, value) {
    a <- d
    # `resistance` is a driver now, so perturbing it means rebuilding the
    # root_network driver rather than rebuilding the supply object -- which is the
    # whole of the consistency change, seen from the caller's side.
    if (identical(par, "resistance")) {
      a$root_network <- series_resistance(value)
    } else {
      a[[par]] <- value
    }
    x <- do.call(grid_leaf_solve, a)
    c(A = x$A, gc = x$gc, psi_leaf = x$psi_leaf, collar = x$collar,
      profit = x$profit)
  }

  g <- do.call(grid_leaf_gradient,
               c(d, list(pars = c("leaf_specific_conductance_max",
                                  "resistance"))))
  expect_identical(g$status, "interior")

  for (p in c("leaf_specific_conductance_max", "resistance")) {
    v <- if (identical(p, "resistance")) 1e4 else kmax_default
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

test_that("a pinned point's non-trait rows are the bound's too", {
  # Both move the FEASIBLE COLLAR INTERVAL and not only the profit function --
  # resistance directly, the conductance through the supply curve -- so their
  # rows at a pinned point run through the bound, as a trait's do.
  dry <- list(psi_soil = 5.5, PPFD = 900, atm_vpd = 3.0,
              supply = leaf_supply_singlelayer(),
              root_network = series_resistance(1e4))
  g <- do.call(grid_leaf_gradient,
               c(dry, list(pars = c("leaf_specific_conductance_max",
                                    "resistance"))))
  expect_identical(g$status, "pinned")
  expect_true(all(is.finite(g$gradient)))
  # Against the solve, absolutely: gc is near zero here, and the difference
  # reads it with ~2e-9 of noise.
  ref <- solve_difference(c(dry, list(model = "TF24", route = "collar")),
                          "leaf_specific_conductance_max")
  expect_lt(max(abs(g$gradient["leaf_specific_conductance_max", ] - ref[1, ])),
            1e-5)
  expect_equal(g$gradient["leaf_specific_conductance_max", "psi_leaf"],
               ref[1, "psi_leaf"], tolerance = 1e-4)
})

test_that("a gradient leaves the leaf's stem curve as it found it", {
  # A gradient in a stem trait reads the curve at the seated traits and
  # perturbs nothing, so every later solve is on the curve the traits name.
  # BIT-identity is the right test: anything a call left behind would show.
  d <- list(psi_soil = 1.5, PPFD = 900, atm_vpd = 2.0)
  before <- do.call(grid_leaf_solve, d)

  invisible(do.call(grid_leaf_gradient, c(d, list(pars = "stem_P50"))))
  invisible(do.call(grid_leaf_gradient, c(d, list(pars = c("stem_P50", "vcmax_25")))))
  # ...and one that names the stem trait last.
  invisible(do.call(grid_leaf_gradient, c(d, list(pars = c("vcmax_25", "stem_P50")))))

  expect_identical(do.call(grid_leaf_solve, d), before)
})

test_that("leaf_gradient() rejects bad arguments", {
  expect_error(grid_gradient(2.0, pars = "not_a_trait"),
               "cannot differentiate")
  # `resistance` is a real parameter -- on the other supply path -- so the
  # message says which path it belongs to rather than calling it unknown.
  expect_error(grid_gradient(2.0, pars = "resistance"),
               "single-potential path only")
  expect_error(grid_gradient(2.0, traits = list(vcmax_25 = 96)),
               "leaf_traits")
  # The differencing controls are gone, not ignored: there is no step and no
  # method to choose, so naming one is an error.
  for (a in c("step", "method", "stationarity_tol", "fast_stem_curve")) {
    expect_error(do.call(leaf_gradient, stats::setNames(list(2.0, 1e-6),
                                                         c("psi_soil", a))),
                 "unused argument", label = a)
  }
})

test_that("the gradient is reported for every output the fit needs", {
  # a companion study fits three responses -- A, gs and psi_leaf -- so all three
  # are differentiated, not just A. `collar` comes along because it is psi*, and
  # `profit` because it is what a DEMOGRAPHIC consumer bills: plant's carbon is
  # `leaf.profit_`, not `assim_colimited_`, so without it the four columns above
  # -- the calibration set -- were disjoint from plant's (#87).
  g <- grid_gradient(2.0, pars = "vcmax_25")
  expect_identical(colnames(g$gradient),
                   c("A", "gc", "psi_leaf", "collar", "profit"))
  expect_identical(names(g$value),
                   c("A", "gc", "psi_leaf", "collar", "profit"))
  expect_true(all(is.finite(g$gradient)))
  # Raising vcmax_25 raises assimilation and opens the stomata, and the leaf pays
  # for it with a more negative water potential (a larger positive magnitude).
  expect_gt(g$gradient["vcmax_25", "A"], 0)
  expect_gt(g$gradient["vcmax_25", "gc"], 0)
  expect_gt(g$gradient["vcmax_25", "psi_leaf"], 0)
  # More photosynthetic capacity is worth having: profit is A minus the hydraulic
  # cost, and the extra carbon exceeds the extra cost. It is BELOW dA/dvcmax_25,
  # which is the statement that the cost is not zero -- the cheapest available
  # check that this column is the objective and not a copy of A.
  expect_gt(g$gradient["vcmax_25", "profit"], 0)
  expect_lt(g$gradient["vcmax_25", "profit"], g$gradient["vcmax_25", "A"])
})

test_that("profit's gradient is the direct term alone at an interior optimum", {
  # The envelope theorem, which is the ONE place this package uses it. At an
  # interior optimum dprofit/dpsi == 0, so the indirect term
  # (dprofit/dpsi)(dpsi*/dtheta) vanishes IDENTICALLY, and dprofit/dtheta is the
  # direct partial at fixed psi. The surface encodes that by reading profit at
  # the collar HELD at psi*, so the column is the direct term and
  # dY_dpsi["profit"] is the analytic zero; the check is that the reported
  # column equals a direct central difference with the collar held at psi*.
  #
  # ⚠️ THE OPERATING POINT IS CHOSEN, NOT ARBITRARY, and choosing it is what makes
  # this a test rather than a formality. The assignment only matters where the term
  # it removes is big enough to see, and that term is NOISE rather than an h^2
  # truncation: `profit` is the maximum, so it is flat, and a central difference
  # of it divides the solve's ~1e-9 floor by a ~1e-6 step. So the point has to be
  # one where that noise is LARGE; at the suite's usual `grid_drivers(2.0)` the same
  # term is ~1e-10 and this test would pass either way, which is the version of it
  # written first.
  #
  # ⚠️ THE POINT DOES NOT SURVIVE A RESULTS CHANGE, so RE-DERIVE IT rather than
  # adjusting the bounds whenever this fails on the `fd` assertion: sweep the golden
  # grid's interior rows and take the largest |fd|.
  #
  #   for each (psi_soil, ppfd, vpd, layers) in the golden grid:
  #     g <- leaf_gradient(..., pars = "vcmax_25"); keep if g$status == "interior"
  #     solve, then fd <- central difference of profit_ in the collar at psi*
  #   pick argmax |fd|
  #
  # A failure saying `abs(fd) <= 1e-5` means the point went quiet, NOT that the
  # feature broke. A failure on `abs(exact)` is the real signal.
  #
  # ⚠️ Those figures are macOS/arm64's, and the second half of this test says why
  # that matters. Read it before adding an assertion on a magnitude here.
  d <- grid_drivers(0.5, ppfd = 1500, vpd = 2.0, layers = 1L)
  g <- do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25")))
  expect_identical(g$status, "interior")

  l <- leaf_model(leaf_traits(), leaf_control(), leaf_supply_multilayer())
  do.call(grid_set_drivers, c(list(l), d))
  l$find_root_collar_psi()
  psi_star <- l$opt_root_psi_

  # The two instruments, both read on a leaf still at base traits.
  hp <- max(abs(psi_star), 1) * 1e-6
  exact <- l$dprofit_droot_collar_psi(psi_star)
  l$evaluate_root_collar_psi(psi_star + hp)
  hi <- l$profit_
  l$evaluate_root_collar_psi(psi_star - hp)
  lo <- l$profit_
  fd <- (hi - lo) / (2 * hp)

  # The answer itself, against a direct difference at fixed collar. THIS IS THE
  # PART THAT HOLDS EVERYWHERE -- it is the claim the feature makes, and it does
  # not depend on how big the term that was dropped happens to be.
  v <- leaf_traits()$vcmax_25
  h <- max(abs(v), 1) * 1e-6
  at <- function(x) {
    set_traits(l, leaf_traits(vcmax_25 = x))
    do.call(grid_set_drivers, c(list(l), d))
    l$evaluate_root_collar_psi(psi_star)
    l$profit_
  }
  direct <- (at(v + h) - at(v - h)) / (2 * h)
  expect_equal(g$gradient["vcmax_25", "profit"], direct, tolerance = 1e-7)

  # ⚠️ THE MAGNITUDES BELOW ARE PLATFORM-SPECIFIC, AND THE FIRST VERSION OF THIS
  # TEST ASSERTED THEM EVERYWHERE. It passed on macOS/arm64 and failed on Linux
  # CI, where the point then in use gave `exact` = 6.6e-11 and `fd` = 8.0e-10
  # rather than macOS's 2.4e-15 and 2.1e-04. At this point macOS gives 1.2e-14
  # and 3.3e-04.
  #
  # That is not a different answer, it is a different NOISE FLOOR: which side of
  # the collar solver's tolerance the root-find lands on is set by libm's exp/pow,
  # and those are not reproducible between Apple's arm64 libm and glibc on
  # x86-64. The package already has a name for that -- it is the same reason the
  # golden files compare bit-exactly on one platform and with a tolerance
  # elsewhere -- so the same predicate gates it here.
  #
  # ⚠️ So this test has TEETH only on the platform the measurement was made on.
  # Off it, the check above still holds and the one below is skipped; a Linux-only
  # run would not catch the assignment being removed. Said plainly rather than left
  # for someone to discover from a green CI.
  skip_if_not(golden_bit_exact_platform(),
              "the noise floor these two numbers measure is macOS/arm64's")
  # The exact instrument says stationary to solver precision; the differenced one
  # says nothing of the kind. Asserting BOTH is the point -- it is the difference
  # between "the dropped term is small" and "the dropped term is unmeasurable by
  # the route that would have supplied it".
  #
  # ⚠️ THE BOUND ON `exact` IS THE SOLVER FLOOR (~1e-9), not a landing point. Pinning
  # it tighter pins which side of its own tolerance the collar root-find lands on,
  # which moves for reasons unrelated to this feature.
  #
  # What carries the test is the GAP, and it is asserted as a RATIO rather than as
  # a magnitude on `fd`.
  #
  # ⚠️ `fd`'s absolute size is not a property of the feature -- it is where the
  # collar root-find happened to land inside its own tolerance, divided by the
  # step. The (P50, c) reparameterisation moved every parameter by ~5e-08 and
  # `fd` at this point went from 2.7e-04 to 3.6e-10, four orders, with nothing
  # about the model changed. A bound on `fd` therefore pins a coincidence; the
  # ratio is the claim. Measured across eight operating points it runs 1e+04 to
  # 1e+09, so 1e+03 is inside every one of them.
  expect_lt(abs(exact), 1e-9)
  expect_gt(abs(fd) / abs(exact), 1e3)
})

# ---------------------------------------------------------------------------
# Reusing a Leaf across gradients (#52)
# ---------------------------------------------------------------------------

test_that("leaf_gradient(x =) matches building a leaf per call", {
  # The whole claim: reuse is a pure cost saving. Bit-identical, not merely close --
  # the setter applies the same traits and drivers either way, so anything else means
  # the reused leaf carried state into the answer.
  # ⚠️ stem_b is raised, not lowered. #38: the vulnerability spline is built over
  # [0, b*log(100)^(1/c)] and never consults psi_crit, so shrinking stem_b below
  # ~3.4 at the default psi_crit = 5.87 puts the solve outside its own domain and it
  # dies naming neither. stem_b = 3.1 does exactly that; 4.2 gives a domain of 7.4.
  tr <- leaf_traits(vcmax_25 = 105, stem_P50 = 4.2)
  args <- list(psi_soil = 2.0, PPFD = 900, atm_vpd = 1.5, traits = tr,
               pars = c("vcmax_25", "stem_P50", "TF24_cost_scale"))
  fresh <- do.call(grid_leaf_gradient, args)

  l <- leaf_model(traits = tr, supply = leaf_supply_multilayer())
  reused <- do.call(grid_leaf_gradient, c(args, list(x = l)))
  expect_identical(reused$gradient, fresh$gradient)
  expect_identical(reused$value, fresh$value)
  expect_identical(reused$status, fresh$status)
  expect_identical(reused$H, fresh$H)

  # And again on the SAME object, which is the point of it: a second call must not
  # inherit the first call's last perturbation.
  again <- do.call(grid_leaf_gradient, c(args, list(x = l)))
  expect_identical(again$gradient, fresh$gradient)
})

test_that("a reused leaf is differentiated at `traits`, not at its own state", {
  # `x` is a vessel. A leaf built at OTHER traits must give the same answer as a
  # fresh one, because the setter applies `traits` before anything is read. If it
  # did not, `value` and `psi_star` -- and so the whole gradient -- would describe
  # the wrong point, plausibly and with no symptom. This is the failure mode that
  # made seating the leaf with reset() rather than set_drivers() load-bearing.
  tr <- leaf_traits(vcmax_25 = 105, stem_P50 = 4.2)
  args <- list(psi_soil = 2.0, PPFD = 900, traits = tr, pars = "vcmax_25")
  fresh <- do.call(grid_leaf_gradient, args)

  # Deliberately different, and deliberately still inside the #38 domain -- the
  # vessel's traits are overwritten before any solve, but leaf_model() does build
  # its splines, and an invalid vessel would confuse the diagnosis if this failed.
  wrong <- leaf_model(traits = leaf_traits(vcmax_25 = 60, stem_P50 = 4.6), supply = leaf_supply_multilayer())
  reused <- do.call(grid_leaf_gradient, c(args, list(x = wrong)))
  expect_identical(reused$gradient, fresh$gradient)
  expect_identical(reused$value, fresh$value)
})

test_that("a reused leaf is differentiated on the route's model, not the one it carries", {
  # Each route seats its own model before its closures read the seated curve. A
  # leaf arriving seated on another curve would otherwise be solved on one curve
  # and tested against the other's first-order condition, and read as pinned.
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
  # And it is handed back on the route's model, solved there: TF24's collar.
  tf <- leaf_model(traits = tr, supply = leaf_supply_multilayer())
  set_drivers(tf, psi_soil = 1.0, PPFD = 900)
  tf$set_model("TF24", "collar")
  tf$optimise()
  expect_identical(operating_point(seated), operating_point(tf))
})

test_that("a leaf-route gradient leaves the caller's leaf solved on that route", {
  # The base point is restored on the route's own model: a collar solve here would
  # hand back a leaf on the leaf route holding a different operating point.
  tr <- leaf_traits()
  sup <- leaf_supply_singlelayer()
  l <- leaf_model(traits = tr, supply = sup)
  set_drivers(l, psi_soil = 1.0, PPFD = 900)
  l$set_model("TF24", "leaf")
  l$optimise()
  before <- operating_point(l)
  leaf_gradient(psi_soil = 1.0, PPFD = 900, x = l, traits = tr, model = "TF24",
                route = "leaf", pars = "vcmax_25")
  expect_identical(operating_point(l), before)
})

test_that("leaf_gradient(x =) leaves the caller's leaf solved at the base point", {
  # hazard 8: an output no path rewrites reads as the previous caller's. The setter
  # leaves the leaf at the last perturbation, so this function restores it.
  tr <- leaf_traits()
  l <- leaf_model(traits = tr, supply = leaf_supply_multilayer())
  set_drivers(l, psi_soil = 2.0, PPFD = 900)
  l$find_root_collar_psi()
  before <- operating_point(l)

  leaf_gradient(psi_soil = 2.0, PPFD = 900, x = l, traits = tr, pars = "vcmax_25", model = "TF24", route = "collar")
  expect_identical(operating_point(l), before)

  # Restored on an ERROR path too -- there are several stop()s downstream of the
  # first perturbation, and a half-perturbed leaf handed back would be worse than
  # the error.
  expect_error(
    leaf_gradient(psi_soil = 2.0, PPFD = 900, x = l, traits = tr,
                  pars = "not_a_parameter", model = "TF24", route = "collar"),
    "cannot differentiate")
  expect_identical(operating_point(l), before)
})

test_that("leaf_gradient() refuses the combinations that would disagree", {
  tr <- leaf_traits()
  l <- leaf_model(traits = tr, supply = leaf_supply_multilayer())
  expect_error(leaf_gradient(psi_soil = 2, x = l, model = "TF24", route = "collar"), "`traits` must be given with")
  expect_error(leaf_gradient(psi_soil = 2, x = l, traits = tr,
                             control = leaf_control(), model = "TF24", route = "collar"), "`control` comes from")
  expect_error(leaf_gradient(psi_soil = 2, x = l, traits = tr,
                             supply = leaf_supply_multilayer(), model = "TF24", route = "collar"),
               "`supply` comes from")
  expect_error(leaf_gradient(psi_soil = 2, x = leaf_traits(), traits = tr, model = "TF24", route = "collar"),
               "must be a Leaf")

  # The supply path is read off the object, so `resistance` is differentiable when
  # `x` is on the single-potential path and not when it is not -- without the caller
  # naming the path twice.
  s <- leaf_model(supply = leaf_supply_singlelayer())
  g <- leaf_gradient(psi_soil = 1.5, x = s, traits = tr,
                     root_network = series_resistance(1e3), pars = "resistance", model = "TF24", route = "collar")
  expect_true(is.finite(g$gradient["resistance", "A"]))
  expect_error(leaf_gradient(psi_soil = 2, x = l, traits = tr, pars = "resistance", model = "TF24", route = "collar"),
               "single-potential path only")
})

test_that("the setter's positional trait call cannot drift in arity", {
  # Traits are applied POSITIONALLY (hazard 13): set_traits() and C++'s
  # `gradient::apply()` pass them in the C++ setter's argument order, and the
  # gradient's slots take their order from the C++ enumeration. A trait added to
  # one list and not the others would silently pass the wrong value for every
  # argument after it, so the three lists are asserted to agree.
  expect_length(leaf_traits(), 16L)
  expect_length(formals(leaf_model(supply = leaf_supply_multilayer())$set_traits), 16L)
  expect_identical(names(leaf_traits()), names(formals(leaf_model(supply = leaf_supply_multilayer())$set_traits)))
  expect_identical(names(leaf_traits()),
                   gradient_par_names()[seq_along(leaf_traits())])
})

# --- a collar potential the caller supplies (#88) -----------------------------

test_that("a prescribed psi at psi* reproduces the solving path", {
  # ⚠️ THE LOAD-BEARING TEST OF THE WHOLE FEATURE. Give the prescribed path back
  # the collar the solve found, and the response of that collar the solve
  # derived, and it must give the solving path's rows: both are the same
  # surface, one with the collar placed by the first-order condition and one
  # with it carried by the caller's dpsi_dtheta. Equal to rounding, not to the
  # bit: the solving path reads profit at the held collar, where the prescribed
  # one reads it at the live collar and adds dY/dpsi x dpsi/dtheta. Measured
  # worst 1.4e-14 relative.
  d <- grid_drivers(2.0)
  pars <- c("vcmax_25", "stem_P50", "TF24_cost_scale")
  a <- do.call(grid_leaf_gradient, c(d, list(pars = pars)))
  expect_identical(a$status, "interior")

  b <- do.call(grid_leaf_gradient,
               c(d, list(pars = pars, psi = a$psi, dpsi_dtheta = -a$M / a$H)))
  expect_identical(b$status, "prescribed")
  expect_equal(b$gradient, a$gradient, tolerance = 1e-12)
  expect_equal(b$M, a$M, tolerance = 1e-12)
  # The curvature too, the other coefficient of a tracking caller's ODE.
  expect_equal(b$H, a$H, tolerance = 1e-9)
  # profit's dY/dpsi is the analytic zero on the solving path and the
  # marginal at the live collar here: the solve's residual.
  outs <- c("A", "gc", "psi_leaf", "collar")
  expect_equal(b$dY_dpsi[outs], a$dY_dpsi[outs], tolerance = 1e-12)
  expect_identical(a$dY_dpsi[["profit"]], 0)
  expect_lt(abs(b$dY_dpsi[["profit"]]), 1e-10)
  expect_identical(b$psi, a$psi)

  # `collar` IS psi, so its column must be exactly what was supplied -- the one
  # column whose value the caller controls outright.
  expect_identical(unname(b$gradient[, "collar"]), unname(-a$M / a$H))

  # And the equivalence is not an accident of the default: dpsi_dtheta = 0 gives
  # the direct terms, which differ from the solving path by the indirect term.
  z <- do.call(grid_leaf_gradient, c(d, list(pars = pars, psi = a$psi)))
  expect_identical(unname(z$gradient[, "collar"]), rep(0, length(pars)))
  expect_true(all(abs(z$gradient[, "A"] - a$gradient[, "A"]) > 1e-12))
})

test_that("a prescribed psi off the optimum matches a difference at fixed collar", {
  # The claim the feature makes: at a collar the caller chose, dY/dtheta is the
  # direct partial. Checked against a central difference computed the long way,
  # at a psi deliberately away from psi* so the indirect term the solving path
  # would add is large and its absence is visible.
  d <- grid_drivers(2.0)
  a <- do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25")))
  psi_off <- a$psi + 0.2
  g <- do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25", psi = psi_off)))
  expect_identical(g$status, "prescribed")

  l <- leaf_model(leaf_traits(), leaf_control(), leaf_supply_multilayer())
  v <- leaf_traits()$vcmax_25
  h <- max(abs(v), 1) * 1e-6
  at <- function(x, field) {
    set_traits(l, leaf_traits(vcmax_25 = x))
    do.call(grid_set_drivers, c(list(l), d))
    l$evaluate_root_collar_psi(psi_off)
    l[[field]]
  }
  fields <- c(A = "assim_colimited_", profit = "profit_")
  for (column in names(fields)) {
    ref <- (at(v + h, fields[[column]]) - at(v - h, fields[[column]])) / (2 * h)
    expect_equal(g$gradient["vcmax_25", column], ref, tolerance = 1e-7,
                 label = column)
  }

  # profit's dY/dpsi is the exact dprofit/dpsi here, far from zero because the
  # point is not stationary: the regime this exists for. The surface's
  # marginal against the leaf's own, to rounding.
  set_traits(l, leaf_traits())
  do.call(grid_set_drivers, c(list(l), d))
  expect_equal(g$dY_dpsi[["profit"]], l$dprofit_droot_collar_psi(psi_off),
               tolerance = 1e-12)
  expect_gt(abs(g$dY_dpsi[["profit"]]), 1)
  # And H is the curvature at the collar given, not at the optimum.
  hp <- 1e-5
  curv <- (l$dprofit_droot_collar_psi(psi_off + hp) -
             l$dprofit_droot_collar_psi(psi_off - hp)) / (2 * hp)
  expect_equal(g$H, curv, tolerance = 1e-6)
})

test_that("a psi outside the feasible interval is refused, with no gradient", {
  # ⚠️ NOT AN ERROR, AND NOT A PROJECTION. The leaf refuses a collar outside its
  # feasible interval rather than moving it onto the end, where the outputs
  # would follow the BOUND instead of the caller's dpsi_dtheta. A fit reaches
  # these points routinely, so the row is reported, not thrown.
  d <- grid_drivers(2.0)
  pars <- c("vcmax_25", "stem_P50")

  g <- do.call(grid_leaf_gradient, c(d, list(pars = pars, psi = 99)))
  expect_identical(g$status, "refused")
  expect_identical(g$psi, 99)
  expect_true(all(is.na(g$gradient)))
  expect_true(all(is.na(g$M)))
  expect_true(all(is.na(g$dY_dpsi)))
  # There is no operating point to describe.
  expect_true(all(is.na(g$value)))

  # Just inside the end the collar is feasible and has a value, but the leaf
  # potential it implies lies in the stem curve's tail, where inverting the
  # flux for it divides by a conductance near e^-40: every row would be
  # rounding (measured, 1e17 in psi_leaf). So that point has a value and no
  # gradient. The end, bisected: the last psi with a value.
  a <- do.call(grid_leaf_gradient, c(d, list(pars = pars)))
  at <- function(p) do.call(grid_leaf_gradient, c(d, list(pars = pars, psi = p)))
  inside <- a$psi
  outside <- 99
  for (k in 1:60) {
    m <- (inside + outside) / 2
    if (all(is.finite(at(m)$value))) inside <- m else outside <- m
  }
  edge <- at(inside)
  expect_identical(edge$status, "no-gradient")
  expect_true(all(is.finite(edge$value)))
  expect_true(all(is.na(edge$gradient)))
  # ...while a psi comfortably inside it is not.
  expect_identical(at(a$psi)$status, "prescribed")
})

test_that("an integer psi is the same psi as its double", {
  # An integer `psi` -- `psi = 3L`, or the entirely ordinary `for (p in 2:5)`,
  # since `2:5` yields integers -- is coerced by the checker, so C++ receives
  # the double.
  d <- grid_drivers(2.0)
  a <- do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25")))
  whole <- 3                        # inside the feasible interval at this point
  dbl <- do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25", psi = whole)))
  int <- do.call(grid_leaf_gradient,
                 c(d, list(pars = "vcmax_25", psi = as.integer(whole))))
  expect_identical(dbl$status, "prescribed")
  expect_identical(int$status, "prescribed")
  expect_identical(int$gradient, dbl$gradient)
  expect_identical(int$psi, dbl$psi)

  # And the batch entry point agrees on the same input, which it did not before.
  b <- leaf_batch(psi_soil = 2.0, PPFD = 900, supply = leaf_supply_multilayer())
  expect_identical(
    leaf_gradient_batch(b, pars = "vcmax_25",
                        psi = as.integer(whole), model = "TF24", route = "collar")$status, "prescribed")
})

test_that("a prescribed collar on a shut leaf is no-gradient", {
  # The collar the shut-down state seats is feasible, and the solving path has
  # rows there, from the bound it pinned on. A caller who prescribes that same
  # collar asks a different question: how the leaf responds with the collar
  # held. The leaf potential then lies in the stem curve's tail, where nothing
  # flows and the inversion for it is rounding, so the point has a value and no
  # gradient, and is not reported as pinned: it was imposed, not solved.
  d <- grid_drivers(16.0)          # past the end of the stem curve: shut down
  solved <- do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25")))
  expect_identical(solved$status, "pinned")

  g <- do.call(grid_leaf_gradient,
               c(d, list(pars = "vcmax_25", psi = solved$psi,
                         dpsi_dtheta = 1)))
  # Feasible -- it is the collar the shut-down state itself seated -- and still
  # without a gradient, which is the whole point of the test.
  expect_identical(g$psi, solved$psi)
  expect_identical(g$status, "no-gradient")
  expect_true(all(is.finite(g$value)))
  expect_true(all(is.na(g$gradient)))
  expect_true(all(is.na(g$dY_dpsi)))

  b <- leaf_batch(psi_soil = 16.0, PPFD = 900, supply = leaf_supply_multilayer())
  expect_identical(leaf_gradient_batch(b, pars = "vcmax_25",
                                       psi = solved$psi, model = "TF24", route = "collar")$status, "no-gradient")
})

test_that("psi and dpsi_dtheta are validated together", {
  d <- grid_drivers(2.0)
  gr <- function(...) do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25", ...)))
  expect_error(gr(psi = -1), "positive")
  expect_error(gr(psi = 0), "positive")
  expect_error(gr(psi = c(1, 2)), "single")
  expect_error(gr(psi = NA_real_), "finite")
  # dpsi_dtheta without psi is meaningless: there the collar's response is
  # derived, not supplied.
  expect_error(gr(dpsi_dtheta = 1), "needs `psi`")
  expect_error(gr(psi = 2, dpsi_dtheta = c(1, 2)), "length 1 or one value")
  expect_error(gr(psi = 2, dpsi_dtheta = c(nope = 1)), "names must be exactly")

  # ⚠️ AN INFINITE `dpsi_dtheta` IS REFUSED, AND THE ENVELOPE IS WHY. A row is
  # `direct + dY_dpsi * dpsi_dtheta`, so the columns come back +-Inf -- except
  # where dY_dpsi is zero, as profit's is at a stationary psi, where Inf x 0 is
  # NaN or, assigned from the direct term, finite and plausible beside them.
  # Reading `profit` alone is plant's own case (#87), so the one surviving
  # column was the one most likely to be believed.
  #
  # `anyNA` already caught NA and NaN; `Inf` it did not, and `is.finite` is the
  # single check that covers all three.
  expect_error(gr(psi = 2, dpsi_dtheta = Inf), "finite")
  expect_error(gr(psi = 2, dpsi_dtheta = -Inf), "finite")
  expect_error(gr(psi = 2, dpsi_dtheta = NaN), "finite")
  # At psi* specifically, which is where the envelope fires and where the
  # asymmetry was, rather than only at an arbitrary prescribed collar.
  solved <- do.call(grid_leaf_gradient, c(d, list(pars = "vcmax_25")))
  expect_error(gr(psi = solved$psi, dpsi_dtheta = Inf), "finite")

  # The batch entry point validates the same way, on BOTH shapes it accepts --
  # the matrix form never reaches `.gradient_dpsi_dtheta()`, so one check there
  # would have left it open.
  b <- leaf_batch(psi_soil = c(2.0, 2.5), PPFD = 900, supply = leaf_supply_multilayer())
  expect_error(leaf_gradient_batch(b, pars = "vcmax_25", psi = c(3, 3),
                                   dpsi_dtheta = Inf, model = "TF24", route = "collar"), "finite")
  expect_error(leaf_gradient_batch(b, pars = "vcmax_25", psi = c(3, 3),
                                   dpsi_dtheta = matrix(c(1, Inf), nrow = 2), model = "TF24", route = "collar"),
               "finite")

  # Named is matched by name, not by position -- the same discipline `pars` has,
  # for the same reason.
  pars <- c("vcmax_25", "stem_P50")
  by_name <- do.call(grid_leaf_gradient,
                     c(d, list(pars = pars, psi = 3.0,
                               dpsi_dtheta = c(stem_P50 = 2, vcmax_25 = 1))))
  by_pos <- do.call(grid_leaf_gradient,
                    c(d, list(pars = pars, psi = 3.0, dpsi_dtheta = c(1, 2))))
  expect_identical(by_name$gradient, by_pos$gradient)
})

# --- per-model parameter availability (plan item 4) --------------------------

test_that("CF77_lambda_ is in the enumeration and refused by other models", {
  # It is Cowan & Farquhar's only parameter. Every other curve's lambda is
  # EMERGENT -- derived from that curve's own parameters rather than set -- so the
  # slot exists once in the enumeration and is available for one model.
  expect_true("CF77_lambda_" %in% gradient_par_names())
  # ⚠️ NOT A FIXED POSITION. `TF24_floor_lambda_o` was appended after it,
  # `LeastCost_beta` after that, and the thermal-damage layer's inputs after all
  # three, which is exactly the "appending is safe" the enumeration is designed
  # for -- so what is asserted is the group and its order, adjacent, not where it
  # sits.
  #
  # ⚠️ AND THE THIRD IS NOT A PRICE. `LeastCost_beta` is a dimensionless ratio of
  # two unit costs; least-cost's price of water is emergent from it. What puts it
  # in this group is only that one curve reads it and no other has the slot.
  i <- match("CF77_lambda_", gradient_par_names())
  expect_identical(gradient_par_names()[i + 0:2],
                   c("CF77_lambda_", "TF24_floor_lambda_o", "LeastCost_beta"))

  # ⚠️ REFUSED, not returned as a zero. A zero means "in this objective but
  # inactive at this operating point", which is the psi_crit case and is
  # information. Structurally absent from the objective is a different statement,
  # and an optimiser handed a zero would read the parameter as inert and converge.
  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900,
                             pars = "CF77_lambda_"),
               "not a parameter of the TF24_floor model")
  # And the message names the axis that ruled it out, so a caller is not left
  # looking at the supply path.
  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900,
                             pars = "CF77_lambda_"),
               "emergent")

  # The default must not offer what it would then refuse.
  expect_false("CF77_lambda_" %in%
                 phylloptim:::.gradient_available_pars(TRUE, "TF24_floor"))
  expect_true("CF77_lambda_" %in%
                phylloptim:::.gradient_available_pars(TRUE, "CF77"))

  # The two axes are independent: the supply path still rules out `resistance`.
  multi <- phylloptim:::.gradient_available_pars(FALSE, "CF77")
  expect_false("resistance" %in% multi)
  expect_true("CF77_lambda_" %in% multi)

  # So does the route: the leaf route pins the stem's upstream end at psi_soil
  # and never reads the soil-to-collar resistance, so it offers no row for it,
  # and naming it says which route does.
  expect_false("resistance" %in%
                 phylloptim:::.gradient_available_pars(TRUE, "TF24_floor"))
  expect_true("resistance" %in%
                phylloptim:::.gradient_available_pars(TRUE, "TF24_floor",
                                                      route = "collar"))
  expect_false("resistance" %in% rownames(
    leaf_gradient(psi_soil = 2.0, PPFD = 900)$gradient))
  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900, pars = "resistance"),
               "reaches only the collar route")
  b <- leaf_batch(psi_soil = 2.0, PPFD = 900)
  expect_error(leaf_gradient_batch(b, pars = "resistance"),
               "reaches only the collar route")

  # An unknown model is refused rather than silently offering everything.
  expect_error(phylloptim:::.gradient_available_pars(TRUE, "Medlyn"),
               "no gradient route")
})

test_that("TF24_floor_lambda_o is the second model-owned price, and refused elsewhere", {
  # ⚠️ THE CLASS, NOT THE INSTANCE. There are two prescribed prices in the
  # enumeration now -- Cowan-Farquhar's and the two-term curve's -- so the thing
  # under test is that `.owned_fields()` maps each to its owner rather than
  # that one name is special. On every other curve the price is EMERGENT, derived
  # from that curve's own parameters rather than set.
  expect_true("TF24_floor_lambda_o" %in% gradient_par_names())

  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900, model = "TF24",
                             pars = "TF24_floor_lambda_o"),
               "not a parameter of the TF24 model")
  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900, model = "TF24",
                             pars = "TF24_floor_lambda_o"),
               "emergent")
  # And the refusal names the OWNER, so a caller is told what to pass rather than
  # only what is wrong.
  expect_error(leaf_gradient(psi_soil = 2.0, PPFD = 900, model = "TF24",
                             pars = "TF24_floor_lambda_o"),
               "TF24_floor")

  # The two prices exclude each other, which is what makes this a table rather
  # than a special case.
  cf <- phylloptim:::.gradient_available_pars(TRUE, "CF77")
  tt <- phylloptim:::.gradient_available_pars(TRUE, "TF24_floor")
  expect_true("CF77_lambda_" %in% cf)
  expect_false("TF24_floor_lambda_o" %in% cf)
  expect_true("TF24_floor_lambda_o" %in% tt)
  expect_false("CF77_lambda_" %in% tt)

  # ⚠️ AND THE HYDRAULIC HALF IS NOT IN THIS CLASS AT ALL. `TF24_floor` reads
  # `TF24_cost_scale`, which is a TRAIT and available on every model -- returning
  # an exactly zero gradient off the curves that do not read it, which is the
  # `psi_crit` case and is information rather than an error. The distinction
  # between "zero" and "refused" is the whole design here, so both are checked.
  expect_true("TF24_cost_scale" %in%
                phylloptim:::.gradient_available_pars(TRUE, "TF24_floor"))
  g <- leaf_gradient(psi_soil = 2.0, PPFD = 900, pars = "JS22_gamma", supply = leaf_supply_multilayer(), model = "TF24", route = "collar")
  expect_identical(unname(g$gradient["JS22_gamma", "A"]), 0)
})

test_that("the curve registry is read from C++, not restated in R", {
  # R selects a curve by POSITION in this vector, so it is compared rather than
  # trusted -- the same discipline gradient_par_names() gets.
  nms <- cost_curve_names()
  expect_identical(nms, c("TF24", "CF77", "JS22", "CMax", "SOX", "JW26",
                          "ProfitMax", "TF24_floor", "LeastCost"))

  # ⚠️ THERE IS NO `cost_curve_has_derivative()` ANY MORE, and its absence is the
  # statement. Every curve is `h(A) - C(psi)` and every derivative is
  # `h'(A)*dA/dpsi - dC/dpsi`, so the predicate returned TRUE for all seven and
  # nothing read it -- a constant dressed as a query. Product objectives are the
  # `log` link and ProfitMax the `A/|A|max` one; neither is a special case.
  expect_false(exists("cost_curve_has_derivative",
                      where = asNamespace("phylloptim"), inherits = FALSE))

  # ⚠️ AND THERE IS NO SECOND LIST OF MODEL NAMES EITHER. `.gradient_link_verified`
  # / `.GRADIENT_VERIFIED_LINKS` recorded which routes had been checked against a
  # finite difference of their own solve; once every route was, the predicate was
  # a constant and the `stop()` behind it unreachable -- and that dead branch still
  # named `$optimise_psi_stem_<curve>()`, which no longer exists. Deleted. What
  # replaces it is derived, so a curve added in C++ arrives here on its own.
  expect_identical(phylloptim:::.gradient_models(), nms)
  expect_false(exists(".gradient_link_verified",
                      where = asNamespace("phylloptim"), inherits = FALSE))

  # Exactly one curve owns `CF77_lambda_`; every other model refuses it, and
  # nothing else is refused on the single-potential path.
  for (m in phylloptim:::.gradient_models()) {
    got <- phylloptim:::.gradient_available_pars(TRUE, m)
    expect_identical("CF77_lambda_" %in% got, identical(m, "CF77"),
                     info = m)
  }

  # The C++ derivative itself answers for a product curve rather than refusing --
  # that is the point of the link -- and it is reached by NAME, through the model
  # the leaf is seated on, with no curve index anywhere.
  l <- leaf_model(supply = leaf_supply_singlelayer())
  set_drivers(l, psi_soil = 1.5, PPFD = 900)
  l$set_model("SOX", "leaf"); l$optimise()
  expect_identical(l$model_curve(), "SOX")
  expect_identical(l$model_route(), "leaf")
  d <- l$dprofit_dpsi_leaf_checked(l$opt_psi_leaf_)
  expect_length(d, 2L)
  expect_true(is.finite(d[[1]]))

  # ⚠️ AND A STEM METHOD REFUSES ON THE COLLAR ROUTE rather than answering about a
  # model the leaf is not seated on. This is what the integer argument used to
  # make unrepresentable-by-accident: an index said which curve, and nothing said
  # which route.
  l$set_model("SOX", "collar")
  expect_error(l$dprofit_dpsi_leaf_checked(l$opt_psi_leaf_), "LEAF-route")
  expect_error(l$evaluate_psi_leaf_at(2.0), "LEAF-route")
})

# --- differentiating a leaf-route model's optimum (plan item 5) --------------------

test_that("every additive curve's leaf-route gradient matches differencing its solve", {
  # ⚠️ The gradient is the exact one and the difference is the noisy one -- a
  # difference of a solve divides that solve's ~1e-09 floor by the step -- so
  # this bounds the difference's error rather than the gradient's. What it
  # establishes is that the surface differentiates the model the optimiser
  # solved, on the route where the leaf potential is the decision variable.
  tr <- leaf_traits()
  net <- series_resistance(1e4)
  mk <- function() {
    l <- leaf_model(tr, leaf_control(), leaf_supply_singlelayer())
    l$CF77_lambda_ <- 1.5e5   # in band: TF24's lambda runs 9e4-3e5
    l$TF24_floor_lambda_o <- 1.5e5   # the same band, and the same reason
    l
  }
  # ⚠️ `TF24_floor` APPEARS TWICE ON PURPOSE. It is the one curve with a parameter in
  # each half of the decomposition -- `TF24_cost_scale` in `Theta~(psi)` and
  # `TF24_floor_lambda_o` in the price -- and the two reach the first-order
  # condition by different routes: the hydraulic term is a function of the
  # potential alone, while `lambda_o*kmax*f(psi)` carries the supply. A gradient
  # that had them the right way round on one and not the other would pass on
  # either alone. Note the first is a parameter the curve SHARES with TF24, which
  # is exactly the point of the shared parameterisation.
  cases <- list(c("TF24", "TF24_cost_scale"), c("CF77", "CF77_lambda_"),
                c("JS22", "JS22_gamma"), c("CMax", "CMax_a"),
                c("TF24_floor", "TF24_cost_scale"),
                c("TF24_floor", "TF24_floor_lambda_o"))
  n_checked <- 0L
  for (case in cases) {
    m <- case[[1]]
    p <- case[[2]]
    g <- leaf_gradient(psi_soil = 1.5, PPFD = 1500, root_network = net,
                       x = mk(), traits = tr, pars = p, model = m,
                       route = "leaf")
    expect_identical(g$model, m)
    # The leaf solve records its kind, as the collar solve does, so an interior
    # optimum is classified interior rather than read off a residual.
    expect_identical(g$status, "interior")
    a <- list(psi_soil = 1.5, PPFD = 1500, root_network = net,
              supply = leaf_supply_singlelayer(), traits = tr, model = m,
              route = "leaf")
    if (m == "CF77") a$CF77_lambda <- 1.5e5
    if (m == "TF24_floor") a$TF24_floor_lambda_o <- 1.5e5
    ref <- solve_difference(a, p)
    for (out in c("A", "gc", "psi_leaf", "profit")) {
      expect_equal(g$gradient[p, out], ref[p, out], tolerance = 1e-5,
                   label = paste(m, p, out))
    }
    n_checked <- n_checked + 1L
  }
  expect_identical(n_checked, 6L)
})

test_that("the decision variable follows the route", {
  tr <- leaf_traits()
  net <- series_resistance(1e4)
  l <- leaf_model(tr, leaf_control(), leaf_supply_singlelayer())

  collar <- leaf_gradient(psi_soil = 1.5, PPFD = 1500, pars = "TF24_cost_scale",
                          model = "TF24", route = "collar",
                          supply = leaf_supply_multilayer())
  stem <- leaf_gradient(psi_soil = 1.5, PPFD = 1500, root_network = net, x = l,
                        traits = tr, pars = "TF24_cost_scale", model = "TF24",
                        route = "leaf")

  # On the collar route `collar` IS psi*, so its derivative is assigned.
  expect_identical(collar$route, "collar")
  expect_true(is.finite(collar$gradient[1, "collar"]))
  # On a leaf route the collar is never solved for, so NOT finite is the truthful
  # answer -- reading a number there would be reading stale or absent state.
  expect_false(is.finite(stem$gradient[1, "collar"]))
  expect_true(is.finite(stem$gradient[1, "psi_leaf"]))

  # ⚠️ And they are DIFFERENT MODELS, not one restricted to the other: same
  # drivers, same cost curve, different optimum.
  expect_false(isTRUE(all.equal(collar$psi, stem$psi)))
})

test_that("a product objective's gradient survives an h-sweep", {
  # ⚠️ THE SWEEP IS THE TEST, not a single step (hazard 15): a finite difference
  # through a solve has a noise floor set by the solver, and only a gradient
  # that is flat across six decades of step shows it is not reading that floor.
  tr <- leaf_traits()
  net <- series_resistance(1e4)
  mk <- function(v, p50) {
    l <- leaf_model(leaf_traits(vcmax_25 = v, stem_P50 = p50), leaf_control(),
                    leaf_supply_singlelayer())
    set_drivers(l, psi_soil = 1.5, PPFD = 1500, root_network = net)
    l
  }
  n_checked <- 0L
  for (m in c("SOX", "JW26")) {
    for (par in c("vcmax_25", "stem_P50")) {
      base <- if (par == "vcmax_25") 96 else 3.4
      solve_at <- function(x) {
        l <- if (par == "vcmax_25") mk(x, 3.4) else mk(96, x)
        l$set_model(if (m == "TF") "TF24" else m, "leaf")
        l$optimise()
        l$assim_colimited_
      }
      h <- base * 1e-2               # the floor of the V, not the toe
      fd <- (solve_at(base + h) - solve_at(base - h)) / (2 * h)
      g <- leaf_gradient(psi_soil = 1.5, PPFD = 1500, root_network = net,
                         x = mk(96, 3.4), traits = tr, pars = par, model = m,
                         route = "leaf")
      expect_equal(g$gradient[par, "A"], fd, tolerance = 1e-2,
                   label = paste(m, par))
      n_checked <- n_checked + 1L
    }
  }
  expect_identical(n_checked, 4L)
})

test_that("the remaining route restrictions are the topological ones", {
  tr <- leaf_traits()
  # A leaf route needs the single-potential path, and says so rather than
  # returning a number from the wrong topology.
  expect_error(leaf_gradient(psi_soil = 1.5, PPFD = 1500,
                             pars = "TF24_cost_scale", model = "TF24",
                             route = "leaf", supply = leaf_supply_multilayer()),
               "single-potential supply path")
  expect_error(leaf_gradient(psi_soil = 1.5, PPFD = 1500,
                             pars = "TF24_cost_scale", model = "Medlyn", supply = leaf_supply_multilayer()),
               "unknown model")
})

test_that("ProfitMax's normaliser follows the traits, so the gradient is a total", {
  net <- series_resistance(1e4)
  mkv <- function(v) {
    l <- leaf_model(leaf_traits(vcmax_25 = v), leaf_control(),
                    leaf_supply_singlelayer())
    set_drivers(l, psi_soil = 1.5, PPFD = 1500, root_network = net)
    l
  }
  g <- leaf_gradient(psi_soil = 1.5, PPFD = 1500, root_network = net, x = mkv(96),
                     traits = leaf_traits(), pars = "vcmax_25", model = "ProfitMax",
                     route = "leaf")
  expect_identical(g$status, "interior")
  # ⚠️ IT USED TO BE NaN, THEN A PARTIAL, AND IS NOW A TOTAL -- three states, and
  # the middle one is worth knowing about. `|A|max` is state a perturbation
  # destroys (`set_traits` + `set_physiology` clears it: 16.757 -> NaN), so
  # something must restore it. Pinning it at the base point fixed the NaN but froze
  # the normaliser, which made the gradient a PARTIAL: 0.0766 against a true
  # 0.0579, the omitted term 32% of the answer. That was the right call only while
  # `|A|max` was the argmax of a 500-point scan and so piecewise constant in the
  # traits. It is a root-find now, and the surface carries it on the pack
  # (normalisers_at), so the gradient follows the traits through it.
  expect_true(is.finite(g$gradient[1, "A"]))

  # The gradient and a difference of the solve agree, which is the statement
  # that the split is gone.
  ref <- solve_difference(list(psi_soil = 1.5, PPFD = 1500, root_network = net,
                               supply = leaf_supply_singlelayer(),
                               model = "ProfitMax", route = "leaf"),
                          "vcmax_25")
  expect_equal(g$gradient[1, "A"], ref[1, "A"], tolerance = 1e-5)

  # And a step-swept difference of the solve, which is what a fit sees.
  h <- 96 * 1e-2
  solve_A <- function(v) {
    l <- mkv(v)
    l$set_model("ProfitMax", "leaf"); l$optimise()
    l$assim_colimited_
  }
  swept <- (solve_A(96 + h) - solve_A(96 - h)) / (2 * h)
  expect_equal(g$gradient[1, "A"], swept, tolerance = 1e-3)
})

test_that("`route` names the potential a curve is maximised over (decision D3)", {
  sp <- leaf_supply_singlelayer()
  multi <- leaf_supply_multilayer()
  net <- series_resistance(1e4)
  outs <- c("A", "gc", "psi_leaf", "collar", "profit")
  # One set of defaults: TF24_floor on the leaf route over the single potential.
  g0 <- leaf_gradient(1.5, pars = "vcmax_25")
  expect_identical(c(g0$model, g0$route), c("TF24_floor", "leaf"))
  expect_false(is.finite(g0$gradient["vcmax_25", "collar"]))
  expect_identical(leaf_gradient(1.5, root_network = net, pars = "vcmax_25",
                                 model = "SOX")$route, "leaf")
  expect_identical(leaf_solve(1.5, root_network = net, model = "SOX"),
                   leaf_solve(1.5, supply = sp, root_network = net,
                              model = "SOX", route = "leaf"))
  # At its default price, 0, TF24_floor is TF24, to the bit.
  expect_identical(leaf_solve(1.5), leaf_solve(1.5, model = "TF24"))
  expect_identical(
    g0$gradient,
    leaf_gradient(1.5, pars = "vcmax_25", model = "TF24")$gradient)
  # A curve on the collar route, on the multi-layer supply: the gradient is
  # taken at the point leaf_solve() solves, and the collar is its decision.
  for (m in c("JS22", "SOX", "ProfitMax")) {
    g <- leaf_gradient(1.5, pars = "vcmax_25", model = m, route = "collar",
                       supply = multi)
    r <- leaf_solve(1.5, model = m, route = "collar", supply = multi)
    expect_identical(unname(g$value[outs]), unlist(r[1, outs], use.names = FALSE),
                     label = m)
    expect_identical(g$status, "interior", label = m)
    expect_true(is.finite(g$gradient["vcmax_25", "collar"]), label = m)
  }
  b <- leaf_batch(psi_soil = c(1.0, 1.5), supply = multi)
  expect_identical(leaf_gradient_batch(b, pars = "vcmax_25", model = "SOX",
                                       route = "collar")$status,
                   c("interior", "interior"))
  # Refused by name: "collar" is a route, not a model; a leaf route, the
  # default, needs the single potential; and a route is one of two.
  expect_error(leaf_gradient(1.5, pars = "vcmax_25", model = "collar"),
               "unknown model \"collar\"")
  expect_error(leaf_solve(1.5, model = "SOX", supply = multi),
               "route = \"collar\"")
  expect_error(leaf_solve(1.5, model = "SOX", route = "root"),
               "`route` must be")
  expect_error(leaf_gradient_batch(b, pars = "vcmax_25", model = "SOX"),
               "single-potential supply path")
})

test_that("leaf_gradient() checks what C++ would otherwise read silently", {
  # An empty request is a 0 x 5 gradient, as the batch's is 1 x 0 x 5.
  g <- leaf_gradient(2, pars = character(0), supply = leaf_supply_multilayer(), model = "TF24", route = "collar")
  expect_identical(dim(g$gradient), c(0L, 5L))
  expect_identical(colnames(g$gradient), c("A", "gc", "psi_leaf", "collar",
                                           "profit"))
  # One observation: a vector driver is refused by name.
  expect_error(leaf_gradient(2, PPFD = c(900, 800), pars = "vcmax_25", supply = leaf_supply_multilayer(), model = "TF24", route = "collar"), "PPFD")
  expect_error(leaf_gradient(2, pars = "CF77_lambda_", supply = leaf_supply_multilayer(), model = "TF24", route = "collar"), "route = \"collar\"")
})
