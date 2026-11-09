# The overview figure's profit panel draws the objective the leaf route
# maximises, so the curve's maximum is the optimum it marks, at any soil
# potential; and the figure hands a caller's leaf back on its own model.

test_that("the profit panel's curve peaks at the optimum it marks", {
  for (psi_soil in c(1.0, 2.0)) {
    m <- leaf_model()
    set_drivers(m, psi_soil = psi_soil, PPFD = 1500)
    grDevices::pdf(NULL)
    panel <- phylloptim:::.fig_profit_panel(m)
    grDevices::dev.off()
    # A 300-point grid brackets the optimum, so its best point is within a
    # grid step's curvature of the solved profit and never above it.
    expect_lte(max(panel$profit), panel$profit_star + 1e-9)
    expect_equal(max(panel$profit), panel$profit_star, tolerance = 1e-4)
  }
})

test_that("plot_model_overview() leaves the caller's leaf on its own model", {
  x <- leaf_model()
  set_drivers(x, psi_soil = 1.5, PPFD = 900)
  x$set_model("CF77", "leaf")
  x$optimise()
  before <- x$profit_
  grDevices::pdf(NULL)
  plot_model_overview(x)
  grDevices::dev.off()
  expect_identical(x$model_curve(), "CF77")
  expect_identical(x$model_route(), "leaf")
  expect_identical(x$profit_, before)
})
