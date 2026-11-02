# The batched trait gradient (#4), which leaf_gradient() is with one
# observation.
#
# WHAT IS BEING PINNED HERE, in order of how much it matters.
#
#   1. That leaf_gradient() and leaf_gradient_batch() hand C++ the same problem.
#      They assemble it separately -- theta from `.gradient_theta()` against
#      `.gradient_theta_matrix()`, one resolved observation against the batch's
#      recycling, a caller's leaf against the batch's -- and a slip in either
#      shows as a difference here, over both supply paths, all three methods and
#      the pinned and shut-down rows. Identical rather than close, since both
#      reach one composite. Whether that composite is RIGHT is
#      test-gradient-referee.R's question, against a difference of the solve.
#
#   2. That the values are what they were. Move the step rule, the stationarity
#      threshold or a solver tolerance and every gradient moves, by less than
#      the referee can see. So `gradient_golden.tsv` records them, over every
#      curve on both routes. ⚠️ It is bit-exact only on macOS/arm64, and
#      compared with the golden file's own per-class tolerance elsewhere, for the
#      reason helper-golden.R gives at length.
#
#   3. That a row the solve cannot handle costs THAT ROW. A proposal during a fit
#      will reach such points; taking out the dataset for one of them would take
#      out the draw.
#
#   4. That the parameter enumeration cannot drift. R passes integer POSITIONS
#      into it, so a reordering silently differentiates the wrong parameter.

# `grid_drivers()` and `batch_drivers()` are in helper-gradient.R, shared with
# test-gradient.R so both files pin the same operating points.

# One row's gradient through the batch, in `leaf_gradient()`'s shape, so the two
# can be compared field by field without every test restating the indexing.
one_row <- function(g) {
  list(gradient = g$gradient[1, , , drop = FALSE],
       value = g$value[1, ],
       status = g$status[[1]],
       method = g$method[[1]],
       H = g$H[[1]],
       stationarity = g$stationarity[[1]])
}

# TF24 on the collar route, the grid's model, unless `args` names one.
expect_entry_points_agree <- function(args, pars, method, label) {
  if (is.null(args$model)) args$model <- "TF24"
  if (is.null(args$route)) args$route <- "collar"
  batch_args <- args[setdiff(names(args), c("model", "route"))]
  ref <- tryCatch(
    do.call(leaf_gradient, c(args, list(pars = pars, method = method))),
    error = function(e) conditionMessage(e))
  if (!is.null(batch_args$psi_soil) && !is.list(batch_args$psi_soil)) {
    batch_args$psi_soil <- list(batch_args$psi_soil)
    if (!is.null(batch_args$soil_depth)) {
      batch_args$soil_depth <- list(batch_args$soil_depth)
    }
  }
  b <- do.call(leaf_batch, batch_args)
  got <- one_row(leaf_gradient_batch(b, pars = pars, method = method,
                                     model = args$model, route = args$route))

  if (is.character(ref)) {
    # leaf_gradient() stops; the batch reports the same refusal per row rather
    # than propagating it, and must not quietly return a number instead.
    expect_identical(got$status, "error", label = label)
    return(invisible(NULL))
  }
  expect_identical(as.vector(got$gradient), as.vector(ref$gradient),
                   label = paste(label, "gradient"))
  expect_identical(unname(got$value), unname(ref$value),
                   label = paste(label, "value"))
  expect_identical(got$status, ref$status, label = paste(label, "status"))
  expect_identical(got$method, ref$method, label = paste(label, "method"))
  expect_identical(got$H, ref$H, label = paste(label, "H"))
  expect_identical(got$stationarity, ref$stationarity,
                   label = paste(label, "stationarity"))
  invisible(ref$status)
}

test_that("the parameter enumeration is the same order in R and in C++", {
  # ⚠️ THE ORDER IS AN INTERFACE AND NOTHING IN THE TYPES ENFORCES IT. R sends
  # integer positions into C++'s `phylloptim::gradient::par_names`, so appending
  # is safe and reordering would differentiate the wrong parameter and report
  # plausible numbers for it. Compared rather than trusted, in both directions.
  r_side <- c(names(leaf_traits()), "leaf_specific_conductance_max",
              "resistance", "CF77_lambda_", "TF24_floor_lambda_o",
              "LeastCost_beta",
              "A_", "I_", "D_", "topt_offset_", "dTopt_max_", "K_A_", "k_i_",
              "k_i_baseline_", "k_rec_", "k_rec_Ea_",
              "c_acclim_maint_", "c_repair_maint_", "c_protect_maint_",
              "c_repair_flux_")
  expect_identical(gradient_par_names(), r_side)
  # And the count, which is what a positional trait call would silently break:
  # sixteen traits, the five that are not traits (two drivers, three model-owned
  # slots), then the fourteen inputs of the thermal-damage layer that reach the
  # solve.
  expect_length(gradient_par_names(), 35L)
  expect_identical(gradient_par_names()[1:16], names(leaf_traits()))
})

test_that("leaf_gradient() is the batch's one observation across the grid", {
  # Eighteen operating points spanning wet to shut-down and one to five soil
  # layers, times all three methods -- so it covers the interior rows the
  # composite runs on, the pinned rows it must decline, the shut-down rows where
  # there is no curvature at all, and the forced-method cases where
  # leaf_gradient() raises an error and the batch reports one per row instead.
  #
  # `stem_c` is in `pars` deliberately: it takes the rebuild route while
  # `stem_P50` beside it moves b and leaves the table.
  pars <- c("vcmax_25", "stem_P50", "TF24_cost_scale", "stem_c")
  rows <- expand.grid(layers = c(1L, 3L, 5L), psi_soil = c(0.5, 2.0, 3.0, 4.0,
                                                          6.0, 16.0),
                      KEEP.OUT.ATTRS = FALSE)
  seen <- character(0)
  for (method in c("auto", "ift", "fd")) {
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      lab <- sprintf("psi_soil=%g layers=%d method=%s", r$psi_soil, r$layers,
                     method)
      # `grid_drivers()`, not `batch_drivers()`: the helper below hands the same
      # arguments to BOTH entry points and puts the profile in a list for the
      # batch itself, since `leaf_gradient()` takes a bare vector.
      st <- expect_entry_points_agree(
        grid_drivers(r$psi_soil, layers = r$layers), pars, method, lab)
      if (method == "auto") {
        seen <- c(seen, st)
      }
    }
  }
  # ⚠️ The grid has to REACH all three classes, or the loop above is only testing
  # the happy path while looking thorough -- which is the way a comparison test
  # fails silently. Asserted by WHERE the classes fall rather than by bare counts,
  # so the claim survives adding a row to the grid: the rows past the end of the
  # stem curve are shut (psi_crit bounds nothing under TF24, #183), and the
  # pinned ones are at the dry end a drought calibration visits.
  expect_length(seen, 18L)
  expect_true(all(seen[rows$psi_soil == 16] == "no-gradient"))
  expect_true(all(rows$psi_soil[seen == "pinned"] %in% c(3, 4, 6)))
  expect_gt(sum(seen == "pinned" & rows$psi_soil < 6), 0L)
  expect_identical(sum(seen == "interior"), 12L)
})

test_that("leaf_gradient() is the batch's one observation on the single-potential path", {
  # The other supply path, and with it the two parameters that are not traits.
  # Both of them move the FEASIBLE COLLAR INTERVAL and not only the profit
  # function -- `resistance` directly, the conductance through the supply curve --
  # so a perturbation can push a bound past psi* more readily than a trait
  # perturbation can, and the refusal detector has to be reached through this
  # route too.
  pars <- c("vcmax_25", "stem_P50", "leaf_specific_conductance_max", "resistance")
  for (psi_soil in c(0.5, 1.5, 3.5, 5.5, 6.5)) {
    for (resistance in c(1e3, 1e4)) {
      for (method in c("auto", "ift", "fd")) {
        args <- list(psi_soil = psi_soil, PPFD = 900, atm_vpd = 3.0,
                     supply = leaf_supply_singlelayer(),
                     root_network = series_resistance(resistance))
        expect_entry_points_agree(
          args, pars, method,
          sprintf("psi_soil=%g r=%g method=%s", psi_soil, resistance, method))
      }
    }
  }
})

test_that("leaf_gradient() is the batch's one observation with the stem curve rebuilt", {
  # `fast_stem_curve = FALSE` rebuilds the vulnerability spline instead of
  # rescaling it, so it is a different code path through `apply()` and
  # not merely a slower one; both entry points must pass it through.
  for (psi_soil in c(1.5, 2.0, 4.0)) {
    for (method in c("auto", "fd")) {
      args <- batch_drivers(psi_soil)
      ref <- do.call(leaf_gradient, c(grid_drivers(psi_soil),
                                      list(pars = "stem_P50", method = method,
                                           fast_stem_curve = FALSE,
                                           model = "TF24", route = "collar")))
      b <- do.call(leaf_batch, args)
      got <- leaf_gradient_batch(b, pars = "stem_P50", method = method,
                                 fast_stem_curve = FALSE,
                                 model = "TF24", route = "collar")
      expect_identical(as.vector(got$gradient[1, , ]),
                       as.vector(ref$gradient),
                       label = sprintf("psi_soil=%g method=%s", psi_soil,
                                       method))
    }
  }
})

test_that("the batch's finite-difference arm takes a step the caller names", {
  # leaf_gradient() differences with a named step and with the route's otherwise;
  # the batch must too, or the two disagree only when a step is named, which the
  # default-step comparisons above cannot see.
  p <- c("vcmax_25", "stem_P50")
  ref <- leaf_gradient(psi_soil = 1.0, PPFD = 900, pars = p, method = "fd",
                       step = 1e-3,
                       supply = leaf_supply_multilayer(), model = "TF24",
                       route = "collar")
  b <- leaf_batch(psi_soil = 1.0, PPFD = 900, supply = leaf_supply_multilayer())
  got <- one_row(leaf_gradient_batch(b, pars = p, method = "fd", step = 1e-3,
                                     model = "TF24", route = "collar"))
  expect_identical(as.vector(got$gradient), as.vector(ref$gradient))
  # The named step reached the arm: it moves the answer off the default's.
  default <- one_row(leaf_gradient_batch(b, pars = p, method = "fd",
                                         model = "TF24", route = "collar"))
  expect_false(identical(as.vector(default$gradient), as.vector(got$gradient)))

  # Under the default "auto" a pinned row takes the same arm, so the named step
  # reaches it there too.
  pinned <- list(psi_soil = 4.5, PPFD = 900, atm_vpd = 3,
                 supply = leaf_supply_multilayer())
  ref_auto <- do.call(leaf_gradient, c(pinned, list(pars = p, step = 1e-3,
                                                    model = "TF24",
                                                    route = "collar")))
  expect_identical(ref_auto$status, "pinned")
  got_auto <- one_row(leaf_gradient_batch(do.call(leaf_batch, pinned), pars = p,
                                          step = 1e-3,
                                          model = "TF24", route = "collar"))
  expect_identical(as.vector(got_auto$gradient), as.vector(ref_auto$gradient))
})

test_that("the recorded gradients have not moved", {
  # The guard the finite-difference referee cannot be: it would pass a change
  # that moves every number by less than the solver's floor. See this file's
  # header, and tools/gradient_golden.R for how to regenerate -- deliberately.
  #
  # ⚠️ Hex floats, read exactly. Do not convert them to decimals: R's decimal
  # parser is not correctly rounded and about 18% of full-precision values come
  # back one ULP out (#13).
  gold <- utils::read.delim(test_path("gradient_golden.tsv"),
                           stringsAsFactors = FALSE)
  # The one list tools/gradient_golden.R wrote the file from.
  cases <- gradient_golden_cases()

  expect_setequal(unique(gold$case), names(cases))
  # ⚠️ BIT-EXACT ONLY ON THE PLATFORM THAT GENERATED IT, macOS/arm64, exactly as
  # tests/cpp/golden/ is -- and the first version of this test asserted equality
  # everywhere, passed here and failed on Linux CI, which is the same mistake
  # test-golden.R's first version made. libm's exp/pow are not bit-reproducible
  # between glibc on x86-64 and Apple's libm on arm64.
  #
  # THE TOLERANCE IS `gradient_golden_tolerance()`, not the solved-output one: a
  # gradient here is a finite difference, so it carries the solve's floor divided by
  # the step, and the smallest-magnitude parameter sets it for the whole file. See
  # that function for the arithmetic.
  #
  # ⚠️ Read the SUMMARY LINE below for a magnitude, never the FAIL lines. The
  # figures in this package's guide were wrong twice because a truncated failure
  # list was read as if it were the distribution.
  worst <- 0
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    rows <- gold[gold$case == nm, ]
    expect_identical(rows$par, cs$pars, label = nm)
    g <- gradient_golden_run(cs)
    expect_identical(g$status[[1]], rows$status[[1]], label = nm)
    expect_identical(g$method[[1]], rows$method[[1]], label = nm)
    for (i in seq_len(nrow(rows))) {
      for (out in c("A", "gc", "psi_stem", "collar", "profit")) {
        got <- g$gradient[1, rows$par[[i]], out]
        expect_golden(got, rows[[out]][[i]], out,
                      paste(nm, rows$par[[i]]),
                      tolerance = gradient_golden_tolerance())
        ref <- as.numeric(rows[[out]][[i]])
        # A stem route's `collar` is NaN: it is not solved for.
        if (is.finite(ref) && ref != 0) {
          worst <- max(worst, abs(got - ref) / abs(ref))
        }
      }
    }
  }
  cat(sprintf("\ngradient_golden.tsv: worst relative difference %.3g (%s)\n",
              worst,
              if (golden_bit_exact_platform()) "bit-exact platform" else
                "compared with tolerance"))
  # And the claim the choice of rows exists to make, stated rather than left
  # implicit in the hex. Where the optimum is pinned at the dry end, psi_stem IS
  # the end of the stem curve, b x_end^(1/c) with x_end = 40 (#183) -- and that
  # moves with stem_P50 by a factor that depends only on stem_c:
  #
  #     dpsi_end/dstem_P50 = (40 / log(2))^(1/stem_c)
  #
  # That is the term a derived end adds to the stem gradient, and the one a
  # hand-written chain rule would be most likely to drop. It is exact at the
  # pinned row, and it is NOT the number at an interior optimum -- which is what
  # makes the comparison a test of the routes rather than of the algebra.
  ratio <- (40 / log(2))^(1 / leaf_traits()$stem_c)
  row <- gold[gold$case == "pinned-dry-1layer" & gold$par == "stem_P50", ]
  expect_equal(as.numeric(row$psi_stem), ratio, tolerance = 1e-6,
               label = "dpsi_stem/dstem_P50 at pinned-dry-1layer")
  interior <- gold[gold$case == "interior-1layer" & gold$par == "stem_P50", ]
  expect_false(isTRUE(all.equal(as.numeric(interior$psi_stem), ratio,
                                tolerance = 1e-3)))

  # The shut-down row sits at its wet bound, the zero-uptake collar, which no
  # stem trait moves (#182): its stem and collar rows in stem_P50 are exactly
  # zero and it moves no gas. Past the end of the stem curve the hydraulic cost
  # has saturated too, so nothing reaches profit either (#183).
  shut <- gold[gold$case == "shutdown-1layer" & gold$par == "stem_P50", ]
  expect_identical(as.numeric(shut$A), 0)
  expect_identical(as.numeric(shut$psi_stem), 0)
  expect_identical(as.numeric(shut$collar), 0)
  expect_identical(as.numeric(shut$profit), 0)
})

test_that("a forced method = \"ift\" at a pinned row is refused by name (#147)", {
  # The batch's own refusal, not leaf_gradient()'s: the two are separate
  # implementations of the same rule, and test-gradient.R checks only the R one.
  # The rows are the pinned ones test-gradient.R names -- two wet, one dry -- and
  # each is checked to be pinned under "auto" first, so the refusal is not read
  # off a row that never reached the branch.
  pinned <- list(list(psi_soil = 4, vpd = 2.0, layers = 5L, ppfd = 900),
                 list(psi_soil = 3, vpd = 3.0, layers = 5L, ppfd = 900),
                 list(psi_soil = 6, vpd = 0.5, layers = 1L, ppfd = 500))
  pars <- c("stem_P50", "vcmax_25")
  for (p in pinned) {
    lab <- sprintf("psi_soil=%g vpd=%g layers=%d", p$psi_soil, p$vpd, p$layers)
    b <- do.call(leaf_batch, batch_drivers(p$psi_soil, ppfd = p$ppfd,
                                           vpd = p$vpd, layers = p$layers))
    auto <- leaf_gradient_batch(b, pars = pars,
                                model = "TF24", route = "collar")
    expect_identical(auto$status, "pinned", label = lab)
    ift <- leaf_gradient_batch(b, pars = pars, method = "ift",
                               model = "TF24", route = "collar")
    expect_identical(ift$status, "error", label = lab)
    expect_match(ift$message, "pinned optimum", label = lab)
    expect_true(all(is.na(ift$gradient)), label = lab)
  }
})

test_that("an unsolvable row costs that row and not the batch", {
  # ⚠️ THE REASON THE BATCH REPORTS A STATUS RATHER THAN THROWING. A proposal
  # during a fit will reach operating points the solve cannot handle -- that is
  # what a proposal distribution does -- and taking out the whole likelihood
  # evaluation for one observation would take out the draw.
  pars <- c("vcmax_25", "stem_P50")
  ok <- c(1.5, 2.0, 2.5)
  b <- leaf_batch(psi_soil = c(ok[[1]], 6.5, ok[[2]], 7.5, ok[[3]]), PPFD = 900,
                  supply = leaf_supply_multilayer())
  g <- leaf_gradient_batch(b, pars = pars, method = "ift",
                           model = "TF24", route = "collar")

  expect_identical(g$status, c("interior", "error", "interior", "error",
                               "interior"))
  # BIT-IDENTICAL to being differentiated alone, which is the real claim: a
  # failed row must not leave the leaf carrying state into its neighbours. The
  # displaced stem-vulnerability spline `perturb_stem_P50()` leaves behind is
  # exactly the kind of thing that would show up here and nowhere else.
  solo <- vapply(ok, function(p) {
    as.vector(leaf_gradient_batch(leaf_batch(psi_soil = p, PPFD = 900,
                                             supply = leaf_supply_multilayer()),
                                  pars = pars, method = "ift",
                                  model = "TF24",
                                  route = "collar")$gradient[1, , ])
  }, numeric(length(pars) * dim(g$gradient)[[3]]))
  for (i in seq_along(ok)) {
    expect_identical(as.vector(g$gradient[2 * i - 1, , ]), solo[, i],
                     label = paste("row", 2 * i - 1))
  }

  # A failed row's gradient is ALL NA rather than partially filled. The parameter
  # loop fails on the first parameter here, but even where it does not, returning
  # the entries computed before the failure alongside `status == "error"` would
  # invite a caller who checks the status per row to read a gradient with a hole
  # in it.
  expect_true(all(is.na(g$gradient[c(2, 4), , ])))
  expect_true(all(nzchar(g$message[c(2, 4)])))
  expect_identical(g$message[c(1, 3, 5)], rep("", 3L))
  # The point's own diagnostics survive, because they describe the operating
  # point rather than the derivative and were determined before the failure.
  expect_true(all(is.finite(g$stationarity[c(2, 4)]) |
                    is.infinite(g$stationarity[c(2, 4)])))
})

test_that("a batch that failed on the fast stem_b path leaves a usable leaf", {
  # ⚠️ `perturb_stem_P50()` leaves the splines built at a different stem_b ON
  # PURPOSE, and `set_traits()` is the way back. A row that threw between the two
  # skipped the restore, so the batch does it in the handler -- and the LAST row
  # is the case that needs it, because there is no next row whose own
  # `set_traits()` would have covered for it (hazard 8).
  b <- leaf_batch(psi_soil = c(2.0, 7.5), PPFD = 900,
                  supply = leaf_supply_multilayer())
  invisible(leaf_gradient_batch(b, pars = "stem_P50", method = "ift",
                                model = "TF24", route = "collar"))
  after <- leaf_gradient_batch(b, pars = "stem_P50",
                               model = "TF24", route = "collar")
  fresh <- leaf_gradient_batch(leaf_batch(psi_soil = c(2.0, 7.5), PPFD = 900,
                                          supply = leaf_supply_multilayer()),
                               pars = "stem_P50",
                               model = "TF24", route = "collar")
  expect_identical(after$gradient, fresh$gradient)
  expect_identical(after$value, fresh$value)

  # And the leaf itself, which the caller can still reach, solves as a fresh one
  # does -- the strongest form of the claim, since a rescaled curve would be
  # invisible in the gradient columns above if both routes carried it.
  set_drivers(b$leaf, psi_soil = 1.5, PPFD = 900)
  b$leaf$find_root_collar_psi()
  expect_identical(operating_point(b$leaf),
                   leaf_solve(psi_soil = 1.5, PPFD = 900,
                              supply = leaf_supply_multilayer(), model = "TF24",
                              route = "collar")[
                     , names(operating_point(b$leaf))])
})

test_that("per-observation theta differentiates each row at its own parameters", {
  # THE HIERARCHICAL CASE, which is what the batch is for: a calibration study maps
  # 40 fitted parameters onto 4 model ones, so the model parameters differ per
  # observation and a shared `theta` would be the wrong shape entirely.
  pars <- c("vcmax_25", "stem_P50")
  psv <- c(1.0, 1.5, 2.0)
  # A non-default stem curve, to keep this off the path every other test takes.
  # The direction no longer matters: psi_crit is a quantile of the same curve, so
  # it scales with stem_P50 and the solve cannot be walked outside its own domain.
  traits <- lapply(c(96, 110, 120), function(v) {
    leaf_traits(vcmax_25 = v, stem_P50 = 3.7)
  })

  b <- leaf_batch(psi_soil = psv, PPFD = 900, supply = leaf_supply_multilayer())
  # ⚠️ THE TRAIT COUNT IS DERIVED, not written down. It was `[1:15]` and a literal
  # `numeric(18)`, both of which had to be found by hand when a trait was added.
  n_traits <- length(leaf_traits())
  n_pars <- length(gradient_par_names())
  theta <- t(vapply(traits, function(tr) {
    unname(c(unlist(tr)[gradient_par_names()[seq_len(n_traits)]], 3.14e-5,
             rep(NA_real_, n_pars - n_traits - 1L)))
  }, numeric(n_pars)))
  g <- leaf_gradient_batch(b, theta = theta, pars = pars,
                           model = "TF24", route = "collar")

  for (i in seq_along(psv)) {
    ref <- leaf_gradient(psi_soil = psv[[i]], PPFD = 900, traits = traits[[i]],
                         pars = pars,
                         supply = leaf_supply_multilayer(), model = "TF24",
                         route = "collar")
    expect_identical(as.vector(g$gradient[i, , ]), as.vector(ref$gradient),
                     label = paste("row", i))
    expect_identical(unname(g$value[i, ]), unname(ref$value),
                     label = paste("row", i))
  }

  # A one-row theta is shared by every observation, which is the other half of the
  # recycling rule and the shape the default takes.
  shared <- leaf_gradient_batch(b, theta = theta[2, , drop = FALSE], pars = pars,
                                model = "TF24", route = "collar")
  all_at_2 <- leaf_gradient_batch(b, traits = traits[[2]], pars = pars,
                                  model = "TF24", route = "collar")
  expect_identical(shared$gradient, all_at_2$gradient)
})

test_that("the default theta describes the point the batch was built at", {
  # `leaf_specific_conductance_max` is both a driver and a fitted parameter, so
  # `theta` has to hold the value the solve will actually use -- the same thing
  # `.gradient_theta()` does for one observation. A default that quietly used the
  # package default instead would differentiate at a point nobody asked about and
  # report plausible numbers.
  kmax <- 5e-5
  b <- leaf_batch(psi_soil = 1.5, PPFD = 900,
                  leaf_specific_conductance_max = kmax,
                  supply = leaf_supply_multilayer())
  g <- leaf_gradient_batch(b, pars = c("vcmax_25",
                                       "leaf_specific_conductance_max"),
                                       model = "TF24", route = "collar")
  ref <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                       leaf_specific_conductance_max = kmax,
                       pars = c("vcmax_25", "leaf_specific_conductance_max"),
                       supply = leaf_supply_multilayer(), model = "TF24",
                       route = "collar")
  expect_identical(as.vector(g$gradient[1, , ]), as.vector(ref$gradient))

  # Non-shared conductances take the n-row path, which is a different branch in
  # `.gradient_theta_matrix()` and has to produce the same answers.
  kv <- c(3.14e-5, 5e-5)
  b2 <- leaf_batch(psi_soil = c(1.5, 1.5), PPFD = 900,
                   leaf_specific_conductance_max = kv,
                   supply = leaf_supply_multilayer())
  expect_identical(b2$theta_nrow, 2L)
  g2 <- leaf_gradient_batch(b2, pars = "leaf_specific_conductance_max",
                            model = "TF24", route = "collar")
  for (i in 1:2) {
    ref2 <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                          leaf_specific_conductance_max = kv[[i]],
                          pars = "leaf_specific_conductance_max",
                          supply = leaf_supply_multilayer(), model = "TF24",
                          route = "collar")
    expect_identical(as.vector(g2$gradient[i, , ]), as.vector(ref2$gradient),
                     label = paste("row", i))
  }
})

test_that("the result is shaped and named for a caller applying a Jacobian", {
  # The layering: C++ returns dY/dtheta for the model parameters and R applies the
  # parameterisation chain rule, vectorised over observations. That only works if
  # the array's dimensions are labelled, so it is asserted rather than assumed.
  pars <- c("vcmax_25", "stem_P50", "TF24_cost_scale")
  b <- leaf_batch(psi_soil = c(1.0, 1.5, 2.0, 2.5), PPFD = 900,
                  supply = leaf_supply_multilayer())
  g <- leaf_gradient_batch(b, pars = pars, model = "TF24", route = "collar")
  expect_identical(dim(g$gradient), c(4L, 3L, 5L))
  expect_identical(dimnames(g$gradient)[[2]], pars)
  expect_identical(dimnames(g$gradient)[[3]],
                   c("A", "gc", "psi_stem", "collar", "profit"))
  expect_identical(dim(g$value), c(4L, 5L))
  expect_identical(colnames(g$value),
                   c("A", "gc", "psi_stem", "collar", "profit"))
  for (f in c("status", "method", "message")) {
    expect_length(g[[f]], 4L)
  }
  for (f in c("H", "stationarity")) {
    expect_length(g[[f]], 4L)
  }
  # `pars` selects, and its ORDER is the array's order -- a fit reads columns by
  # position after this, so a permuted request has to come back permuted rather
  # than in the enumeration's order.
  rev_g <- leaf_gradient_batch(b, pars = rev(pars),
                               model = "TF24", route = "collar")
  expect_identical(dimnames(rev_g$gradient)[[2]], rev(pars))
  expect_identical(g$gradient[, "vcmax_25", ], rev_g$gradient[, "vcmax_25", ])
})

test_that("`pars` order does not change any gradient (#72)", {
  # ⚠️ THE REGRESSION GUARD FOR #72, AND IT IS THE PRIMARY ONE -- not the golden
  # file. The defect was worth up to 3.4e-5 relative, and the golden file is
  # compared with a 1e-3 tolerance off macOS/arm64, so on Linux it could not
  # catch a recurrence. This test compares two runs IN THE SAME PROCESS, so it is
  # exact everywhere and needs no tolerance at all.
  #
  # What went wrong, because the shape of it is what this asserts: the `stem_b`
  # shortcut calls `perturb_stem_P50()`, which rescales the spline and touches
  # nothing else, while the loops only restore base AFTER the whole parameter
  # loop. So a `stem_b` that was not first was differentiated at a point
  # displaced by one step in whichever parameter preceded it.
  #
  # ⚠️ Asserted as the GENERAL invariant -- order changes nothing, for any
  # parameter -- rather than as "stem_b is now reseated". `root_b` obeys the same
  # homogeneity identity and would get the same shortcut, and a stem_b-shaped
  # test would silently stop covering the case it was written for.
  b <- leaf_batch(psi_soil = 1.5, PPFD = 900, supply = leaf_supply_multilayer())
  pars <- c("vcmax_25", "stem_c", "phi_J", "stem_P50", "TF24_cost_scale", "root_P50")

  for (method in c("auto", "fd")) {
    for (fast in c(TRUE, FALSE)) {
      lab <- paste(method, "fast =", fast)
      ref <- leaf_gradient_batch(b, pars = pars, method = method,
                                 fast_stem_curve = fast,
                                 model = "TF24",
                                 route = "collar")$gradient[1, , ]
      # Reversed, and rotated, so the claim is not "one other ordering agrees".
      for (perm in list(rev(pars), pars[c(4, 1, 6, 2, 5, 3)],
                        c("stem_P50", setdiff(pars, "stem_P50")))) {
        got <- leaf_gradient_batch(b, pars = perm, method = method,
                                   fast_stem_curve = fast,
                                   model = "TF24",
                                   route = "collar")$gradient[1, , ]
        expect_identical(got[pars, ], ref[pars, ],
                         label = paste(lab, paste(perm, collapse = ",")))
      }
    }
  }

  # ⚠️ And the predecessor that made the defect WORST, which is not the one you
  # would pick by eye: `a` = 0.30 and `curv_fact_colim` = 0.99 take the absolute
  # 1e-6 step floor rather than a relative step, and `stem_c`'s perturbation
  # REBUILDS the very curve the `stem_b` shortcut then rescales. Every trait is
  # checked as a predecessor rather than a sample of them.
  solo <- leaf_gradient_batch(b, pars = "stem_P50",
                              model = "TF24",
                              route = "collar")$gradient[1, "stem_P50", ]
  for (p in setdiff(names(leaf_traits()), "stem_P50")) {
    got <- leaf_gradient_batch(b, pars = c(p, "stem_P50"),
                               model = "TF24",
                               route = "collar")$gradient[1, "stem_P50", ]
    expect_identical(got, solo, label = paste("after", p))
  }

  # And leaf_gradient(), the same composite, agrees.
  r_ref <- leaf_gradient(psi_soil = 1.5, PPFD = 900,
                         pars = c("TF24_cost_scale", "stem_P50"),
                         supply = leaf_supply_multilayer(), model = "TF24",
                         route = "collar")
  batch <- leaf_gradient_batch(b, pars = c("TF24_cost_scale", "stem_P50"),
                               model = "TF24", route = "collar")
  expect_identical(as.vector(batch$gradient[1, , ]), as.vector(r_ref$gradient))
})

test_that("leaf_batch() recycles its drivers the way leaf_solve() does", {
  # Same helpers, so the rules cannot drift -- but that is only true while both
  # actually call them, which is what this checks.
  b <- leaf_batch(psi_soil = c(1.0, 1.5, 2.0), PPFD = 900,
                  supply = leaf_supply_multilayer())
  expect_identical(b$n, 3L)
  g <- leaf_gradient_batch(b, pars = "vcmax_25",
                           model = "TF24", route = "collar")
  solved <- leaf_solve(psi_soil = c(1.0, 1.5, 2.0), PPFD = 900,
                       supply = leaf_supply_multilayer(), model = "TF24",
                       route = "collar")
  expect_identical(g$value[, "A"], solved$A)
  expect_identical(g$value[, "collar"], solved$collar)

  # A plain numeric psi_soil is N single-layer observations; a list is one
  # multi-layer observation each. Getting that backwards is the easiest mistake to
  # make with this function, and it fails loudly rather than silently, because the
  # network and the profile then disagree in length.
  expect_identical(leaf_batch(psi_soil = list(c(1, 1.5, 2)), PPFD = 900,
                              supply = leaf_supply_multilayer())$n, 1L)
  expect_length(leaf_gradient_batch(
    leaf_batch(psi_soil = list(c(1, 1.5, 2)), PPFD = 900,
               supply = leaf_supply_multilayer()),
    pars = "vcmax_25", model = "TF24", route = "collar")$status, 1L)
})

test_that("leaf_gradient_batch() rejects what it cannot do", {
  b <- leaf_batch(psi_soil = 1.5, PPFD = 900, supply = leaf_supply_multilayer())
  s <- leaf_batch(psi_soil = 1.5, supply = leaf_supply_singlelayer())

  expect_error(leaf_gradient_batch(b, pars = "not_a_trait",
                                   model = "TF24", route = "collar"),
               "cannot differentiate")
  # The same message `leaf_gradient()` gives, from the same function: `resistance`
  # is a real parameter on the other supply path, so saying "unknown" would be
  # misleading rather than merely unhelpful.
  expect_error(leaf_gradient_batch(b, pars = "resistance",
                                   model = "TF24", route = "collar"),
               "single-potential path only")
  expect_true(is.finite(leaf_gradient_batch(s, pars = "resistance",
                                            model = "TF24",
                                            route = "collar")$gradient[
    1, 1, "A"]))

  expect_error(leaf_gradient_batch(b, step = -1,
                                   model = "TF24",
                                   route = "collar"), "positive number")
  expect_error(leaf_gradient_batch(b, method = "magic",
                                   model = "TF24", route = "collar"))
  expect_error(leaf_gradient_batch(leaf_traits()), "must come from leaf_batch")
  expect_error(leaf_batch(psi_soil = 1.5, traits = list(vcmax_25 = 96),
                          supply = leaf_supply_multilayer()),
               "leaf_traits")

  # `traits` and `theta` both say where the gradient is taken, so passing both is
  # refused rather than resolved in favour of one of them.
  th <- matrix(1, nrow = 1, ncol = length(gradient_par_names()))
  expect_error(leaf_gradient_batch(b, traits = leaf_traits(), theta = th,
                                   model = "TF24", route = "collar"),
               "pass one")
  expect_error(leaf_gradient_batch(b, theta = matrix(1, 1, 13),
                                   model = "TF24", route = "collar"),
               paste(length(gradient_par_names()), "columns"))
  expect_error(leaf_gradient_batch(b, theta = th[rep(1, 3), , drop = FALSE],
                                   model = "TF24", route = "collar"),
               "1 row or one per observation")
  expect_error(leaf_gradient_batch(b, theta = as.data.frame(th),
                                   model = "TF24", route = "collar"),
               "numeric matrix")
  # ⚠️ Column names are checked against the enumeration ORDER, because C++
  # indexes by position: a reordered matrix would differentiate the wrong
  # parameters and return plausible numbers rather than fail.
  bad <- th
  colnames(bad) <- rev(gradient_par_names())
  expect_error(leaf_gradient_batch(b, theta = bad,
                                   model = "TF24",
                                   route = "collar"), "in that order")
})

test_that("a batch that has lost its C++ pointer says so rather than crashing", {
  # ⚠️ The prepared drivers and the leaf are both external pointers, so a batch
  # that has been through `serialize()` holds NULL ones -- and that is true within
  # the SAME session too, which is worth knowing: `readRDS()` does not restore a
  # pointer, it restores a null. Serialising a batch is a reasonable thing for
  # someone to try, and "external pointer is not valid" -- what the first read of
  # either would otherwise report, from whichever one Rcpp unwrapped first -- is
  # true and gives no hint that the fix is to rebuild.
  b <- leaf_batch(psi_soil = 1.5, PPFD = 900, supply = leaf_supply_multilayer())
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(b, f)
  expect_error(leaf_gradient_batch(readRDS(f), pars = "vcmax_25",
                                   model = "TF24", route = "collar"),
               "does not survive")
  # The live one still works afterwards: the check is on the handle, not a flag
  # set by serialising, so `saveRDS()` does not poison the object it copied.
  expect_true(is.finite(leaf_gradient_batch(b, pars = "vcmax_25",
                                            model = "TF24",
                                            route = "collar")$gradient[
    1, 1, "A"]))
  # Something that was never a handle at all, which the same check has to catch:
  # dereferencing a double as a pointer is not an error the session survives.
  expect_error(phylloptim:::gradient_batch_check(1.0), "not a prepared batch")
})

# --- a collar potential the caller supplies (#88) -----------------------------

test_that("the batch's prescribed psi agrees with leaf_gradient()", {
  # Both uses of the prescribed path -- the collar it solved for handed back
  # (which must reproduce the solve), and a collar deliberately off it -- and
  # both entry points' handling of `psi` and `dpsi_dtheta`.
  pars <- c("vcmax_25", "stem_P50", "TF24_cost_scale")
  soils <- c(1.0, 2.0, 3.0)
  b <- leaf_batch(psi_soil = soils, PPFD = 900,
                  supply = leaf_supply_multilayer())
  a <- leaf_gradient_batch(b, pars = pars, model = "TF24", route = "collar")
  expect_true(all(a$status == "interior"))

  g <- leaf_gradient_batch(b, pars = pars, psi = a$psi,
                           dpsi_dtheta = -a$M / a$H,
                           model = "TF24", route = "collar")
  expect_true(all(g$status == "prescribed"))
  expect_true(all(g$method == "prescribed"))
  expect_identical(g$gradient, a$gradient)
  expect_identical(g$M, a$M)
  expect_identical(g$dY_dpsi, a$dY_dpsi)

  off <- a$psi + 0.2
  h <- leaf_gradient_batch(b, pars = pars, psi = off,
                           model = "TF24", route = "collar")
  for (i in seq_along(soils)) {
    r <- leaf_gradient(psi_soil = soils[[i]], PPFD = 900, pars = pars,
                       psi = off[[i]],
                       supply = leaf_supply_multilayer(), model = "TF24",
                       route = "collar")
    expect_identical(as.vector(h$gradient[i, , ]), as.vector(r$gradient),
                     label = paste("row", i))
    expect_identical(as.vector(h$M[i, ]), unname(r$M), label = paste("M", i))
    expect_identical(as.vector(h$dY_dpsi[i, ]), unname(r$dY_dpsi),
                     label = paste("dY_dpsi", i))
    expect_identical(h$psi[[i]], r$psi)
  }
})

test_that("a refused row costs that row and not the batch", {
  # The same per-row discipline `status == "error"` has, for a case that is not
  # an error: a fit will hand over collar potentials that a perturbation has put
  # outside the feasible interval, and taking out the whole likelihood
  # evaluation for one of them would take out the draw.
  pars <- c("vcmax_25", "stem_P50")
  soils <- c(1.0, 2.0, 3.0)
  b <- leaf_batch(psi_soil = soils, PPFD = 900,
                  supply = leaf_supply_multilayer())
  a <- leaf_gradient_batch(b, pars = pars, model = "TF24", route = "collar")

  g <- leaf_gradient_batch(b, pars = pars, psi = c(99, a$psi[2:3]),
                           model = "TF24", route = "collar")
  expect_identical(g$status, c("refused", "prescribed", "prescribed"))
  expect_true(all(is.na(g$gradient[1, , ])))
  expect_true(all(is.na(g$M[1, ])))
  # The refused row keeps the psi it was given and has no operating point, and
  # the rows either side are BIT-IDENTICAL to being asked alone -- the leaf must
  # not carry the refused row's state into its neighbours.
  expect_identical(g$psi[[1]], 99)
  expect_true(all(is.na(g$value[1, ])))
  for (i in 2:3) {
    solo <- leaf_gradient_batch(leaf_batch(psi_soil = soils[[i]], PPFD = 900,
                                           supply = leaf_supply_multilayer()),
                                pars = pars, psi = a$psi[[i]],
                                model = "TF24", route = "collar")
    expect_identical(as.vector(g$gradient[i, , ]),
                     as.vector(solo$gradient[1, , ]),
                     label = paste("row", i))
  }
})

test_that("the batch validates psi and dpsi_dtheta by shape", {
  pars <- c("vcmax_25", "stem_P50")
  b <- leaf_batch(psi_soil = c(1.0, 2.0, 3.0), PPFD = 900,
                  supply = leaf_supply_multilayer())
  expect_error(leaf_gradient_batch(b, pars = pars, psi = 2.0,
                                   model = "TF24", route = "collar"),
               "one collar potential per observation")
  expect_error(leaf_gradient_batch(b, pars = pars, psi = rep(-1, 3),
                                   model = "TF24", route = "collar"),
               "finite and positive")
  expect_error(leaf_gradient_batch(b, pars = pars, psi = rep(2, 3),
                                   method = "fd",
                                   model = "TF24", route = "collar"),
               "cannot be given with `psi`")
  expect_error(leaf_gradient_batch(b, pars = pars, dpsi_dtheta = c(1, 2),
                                   model = "TF24", route = "collar"),
               "needs `psi`")
  expect_error(leaf_gradient_batch(b, pars = pars, psi = rep(2, 3),
                                   dpsi_dtheta = matrix(0, 2, 2),
                                   model = "TF24", route = "collar"),
               "3 x 2")

  # A vector is one value per parameter, shared by every observation -- and it
  # must equal the matrix that spells that out, or the recycling is wrong. Each
  # row is given its OWN psi* so none of them is refused; a shared psi across
  # three soil potentials could be, and then the columns compared here would be
  # NA.
  psi <- leaf_gradient_batch(b, pars = pars,
                             model = "TF24", route = "collar")$psi
  shared <- leaf_gradient_batch(b, pars = pars, psi = psi,
                                dpsi_dtheta = c(0.1, 0.2),
                                model = "TF24", route = "collar")
  spelt <- leaf_gradient_batch(b, pars = pars, psi = psi,
                               dpsi_dtheta = matrix(c(0.1, 0.1, 0.1,
                                                      0.2, 0.2, 0.2), 3, 2),
                                                      model = "TF24",
                                                      route = "collar")
  expect_identical(shared$gradient, spelt$gradient)
  expect_identical(unname(shared$gradient[, , "collar"]),
                   matrix(c(0.1, 0.1, 0.1, 0.2, 0.2, 0.2), 3, 2))
})

test_that("the batch refuses a price its model does not read, and a missing one, as leaf_solve() does", {
  net <- series_resistance(1e4)
  single <- leaf_supply_singlelayer()
  solve_error <- function(...) {
    tryCatch(leaf_solve(psi_soil = 1.5, PPFD = 1500, ...),
             error = conditionMessage)
  }
  batch_error <- function(model, ...) {
    b <- leaf_batch(psi_soil = c(1.0, 1.5), PPFD = 1500, ...)
    tryCatch(leaf_gradient_batch(b, pars = "vcmax_25", model = model,
                                 route = "collar"),
             error = conditionMessage)
  }

  # A price the model does not read: the same check, so the same message.
  expect_identical(
    batch_error("TF24_floor", supply = single, root_network = net,
                CF77_lambda = 1.5e5),
    solve_error(supply = single, root_network = net, model = "TF24_floor",
                CF77_lambda = 1.5e5))
  expect_identical(batch_error("TF24", TF24_floor_lambda_o = 1.5e5),
                   solve_error(model = "TF24", TF24_floor_lambda_o = 1.5e5))
  # given on one observation of several is still given
  expect_error(leaf_gradient_batch(
                 leaf_batch(psi_soil = c(1.0, 1.5), PPFD = 1500,
                            LeastCost_beta = c(NA, 146),
                            supply = leaf_supply_multilayer()),
                 pars = "vcmax_25", model = "TF24", route = "collar"),
               "model = \"TF24\" does not read it")

  # An unset price is the curve's default, on every observation it is unset on:
  # TF24_floor's is 0.
  floor_at <- function(price) {
    b <- leaf_batch(psi_soil = c(1.0, 1.5), PPFD = 1500,
                    TF24_floor_lambda_o = price,
                    supply = leaf_supply_multilayer())
    leaf_gradient_batch(b, pars = "vcmax_25", model = "TF24_floor",
                        route = "collar")$gradient
  }
  expect_identical(floor_at(NA), floor_at(0))
  expect_identical(floor_at(c(1.5e5, NA))[2, , ], floor_at(0)[2, , ])

  # A price that varies by observation reaches its own row while kmax and the
  # resistance are shared, so theta has a row per observation.
  run <- function(price) {
    b <- leaf_batch(psi_soil = 1.5, PPFD = 1500, supply = single,
                    root_network = net, CF77_lambda = price)
    leaf_gradient_batch(b, pars = "vcmax_25", model = "CF77", route = "stem")
  }
  both <- leaf_gradient_batch(
    leaf_batch(psi_soil = c(1.5, 1.5), PPFD = 1500, supply = single,
               root_network = net, CF77_lambda = c(1e5, 3e5)),
    pars = "vcmax_25", model = "CF77", route = "stem")
  expect_identical(both$value[, "A"],
                   unname(c(run(1e5)$value[, "A"], run(3e5)$value[, "A"])))
  expect_false(both$value[1, "A"] == both$value[2, "A"])

  # The check reads the columns the solve reads, so a supplied theta is held to
  # it too, and a batch built with one curve's price can be differentiated
  # under another once that column is NA.
  b <- leaf_batch(psi_soil = c(1.0, 1.5), PPFD = 1500, supply = single,
                  root_network = net, CF77_lambda = 1.5e5)
  th <- .gradient_theta_matrix(b, leaf_traits(), model = "CF77")
  colnames(th) <- gradient_par_names()
  expect_error(leaf_gradient_batch(b, pars = "vcmax_25",
                                   model = "TF24", route = "collar"),
               "CF77 curve's prescribed price")
  expect_error(leaf_gradient_batch(b, pars = "vcmax_25", model = "TF24",
                                   theta = th, route = "collar"),
               "CF77 curve's prescribed price")
  th[, "CF77_lambda_"] <- NA_real_
  expect_identical(
    leaf_gradient_batch(b, pars = "vcmax_25",
                        model = "TF24", theta = th, route = "collar")$status,
    c("interior", "interior"))

  # leaf_solve() seats one price for every row
  expect_error(leaf_solve(psi_soil = 1.5, PPFD = 1500, supply = single,
                          root_network = net, model = "CF77",
                          CF77_lambda = c(1e5, 3e5), route = "collar"),
               "must be a single number")
})

test_that("leaf_gradient() is the batch's one observation on every model", {
  # The owned prices reach theta by different routes in the two entry points --
  # off the caller's leaf in one, the batch's arguments in the other -- so a
  # slot one of them misses shows here as a difference.
  tr <- leaf_traits()
  net <- series_resistance(1e4)
  psv <- c(1.0, 1.5, 2.0)
  own <- c(TF24 = "TF24_cost_scale",
           CF77 = "CF77_lambda_", JS22 = "JS22_gamma", CMax = "CMax_a",
           SOX = "vcmax_25", JW26 = "vcmax_25", ProfitMax = "vcmax_25",
           TF24_floor = "TF24_floor_lambda_o",
           # ⚠️ ITS OWN SLOT, NOT A TRAIT: `gradient::apply()` must write it
           # from theta, or the column is an exact zero that no finiteness
           # check sees (test-gradient-referee.R).
           LeastCost = "LeastCost_beta")
  for (m in names(own)) for (rt in c("collar", "stem")) {
    p <- own[[m]]
    # Each model's batch carries its own price only: a price the model does not
    # read is refused, as it is by leaf_solve(). The single-observation leaf
    # holds all three, which covers each slot reaching the right member.
    priced <- function(field, value) {
      if (.owns_field(m, field)) value else NA_real_
    }
    b <- leaf_batch(psi_soil = psv, PPFD = 1500, root_network = net,
                    supply = leaf_supply_singlelayer(), traits = tr,
                    CF77_lambda = priced("CF77_lambda_", 1.5e5),
                    TF24_floor_lambda_o = priced("TF24_floor_lambda_o", 1.5e5),
                    LeastCost_beta = priced("LeastCost_beta", 146))
    gb <- leaf_gradient_batch(b, pars = p, model = m, route = rt)
    one <- vapply(psv, function(ps) {
      l <- leaf_model(tr, leaf_control(), leaf_supply_singlelayer())
      l$CF77_lambda_ <- 1.5e5
      l$TF24_floor_lambda_o <- 1.5e5
      l$LeastCost_beta <- 146
      leaf_gradient(psi_soil = ps, PPFD = 1500, root_network = net, x = l,
                    traits = tr, pars = p, model = m, route = rt)$gradient[1, "A"]
    }, numeric(1))
    expect_identical(unname(gb$gradient[, 1, "A"]), unname(one),
                     label = paste(m, rt))
  }
})
