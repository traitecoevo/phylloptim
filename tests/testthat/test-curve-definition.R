# Each cost curve is one definition, `Curve<K>` in leaf_cost_curves.hpp, and the
# R surface reads its parameters through cost_curve_parameters(). These hold the
# R side to the definition, so a parameter added there and not wired here fails
# by name.

cp <- cost_curve_parameters()

test_that("the parameter list is the curves' and every slot is a gradient slot", {
  expect_named(cp, c("curve", "parameter", "default", "trait", "aligned"))
  expect_true(all(cp$curve %in% cost_curve_names()))
  # In cost_curve_names() order, which is the order R's tables keep.
  expect_identical(unique(cp$curve), intersect(cost_curve_names(), cp$curve))
  expect_false(anyDuplicated(cp$parameter) > 0)
  expect_true(all(cp$parameter %in% gradient_par_names()),
              info = paste(setdiff(cp$parameter, gradient_par_names()), collapse = ", "))
})

test_that("every curve parameter is bound, in the R lists, and at its default", {
  leaf <- leaf_model()
  for (i in seq_len(nrow(cp))) {
    p <- cp$parameter[[i]]
    expect_true(p %in% names(leaf), info = p)
    if (is.na(cp$default[[i]])) {
      expect_true(is.na(leaf[[p]]), info = p)
    } else {
      expect_identical(leaf[[p]], cp$default[[i]], info = p)
    }
  }

  # A trait is a leaf_traits() argument whose default is the definition's.
  traits <- cp[cp$trait, ]
  expect_true(all(traits$parameter %in% names(formals(leaf_traits))))
  expect_identical(unlist(leaf_traits()[traits$parameter]),
                   stats::setNames(traits$default, traits$parameter))

  # An owned field is an argument of leaf_solve() and leaf_batch(), and the
  # table R routes it by is the definition's.
  owned <- cp[!cp$trait, ]
  for (f in owned$parameter) {
    expect_true(sub("_$", "", f) %in% names(formals(leaf_solve)), info = f)
    expect_true(sub("_$", "", f) %in% names(formals(leaf_batch)), info = f)
  }
  expect_identical(phylloptim:::.owned_fields(),
                   stats::setNames(owned$parameter, owned$curve))

  # An aligned default has a bracket for align_curves() to solve it in.
  expect_setequal(names(phylloptim:::.align_brackets),
                  cp$parameter[cp$aligned])
})

test_that("every curve parameter is registered", {
  skip_if_not_installed("yaml")
  fields <- model_registry()$parameters$field
  expect_true(all(cp$parameter %in% fields),
              info = paste(setdiff(cp$parameter, fields), collapse = ", "))
})

test_that("the by-name methods cover every curve, and are each curve's own", {
  l <- leaf_model(supply = leaf_supply_singlelayer())
  l$TF24_floor_lambda_o <- 1.5e5
  set_drivers(l, psi_soil = 1, PPFD = 900)
  l$prepare_profitmax()
  ps <- 1.6
  pu <- 1
  for (curve in cost_curve_names()) {
    v <- c(profit = l$curve_profit(curve, ps, pu),
           cost = l$curve_cost(curve, ps, pu),
           lambda = l$curve_lambda(curve, ps, pu))
    expect_true(all(is.finite(v)), info = curve)
  }
  expect_error(l$curve_cost("TF25", ps, pu), "unknown cost curve")

  # Bit for bit against the frozen per-curve methods. A product's or a ratio's
  # lambda reads the assimilation last evaluated, so the profit goes first.
  expect_identical(l$curve_cost("TF24", ps, pu), l$hydraulic_cost_TF(ps))
  for (curve in c("CF77", "JS22", "CMax", "TF24_floor")) {
    expect_identical(l$curve_cost(curve, ps, pu),
                     l[[paste0("hydraulic_cost_", curve)]](ps, pu), info = curve)
  }
  expect_identical(l$curve_cost("SOX", ps, pu), -log(l$sox_reduction(ps)))
  expect_identical(l$curve_cost("JW26", ps, pu), -log(l$jw26_reduction(ps)))
  for (curve in c("CF77", "JS22", "CMax", "TF24_floor", "SOX", "JW26")) {
    expect_identical(l$curve_profit(curve, ps, pu),
                     l[[paste0("profit_psi_stem_", curve)]](ps, pu), info = curve)
  }
  expect_identical(l$curve_profit("TF24", ps, pu), l$profit_psi_stem_TF(ps, pu))
  expect_identical(l$curve_lambda("TF24", ps, pu), l$lambda_TF24(ps))
  expect_identical(l$curve_lambda("TF24_floor", ps, pu), l$lambda_TF24_floor(ps))
  for (curve in c("JS22", "CMax", "SOX", "JW26")) {
    l$curve_profit(curve, ps, pu)
    expect_identical(l$curve_lambda(curve, ps, pu),
                     l[[paste0("lambda_", curve)]](ps, pu), info = curve)
  }
  expect_identical(l$curve_lambda("CF77", ps, pu), l$CF77_lambda_)
})
