# The golden grid's drivers, shared by the two gradient test files.
#
# A row built with these is the same operating point tests/cpp/test_golden.cpp
# and test-golden.R use, which is what lets a gradient test say "the pinned rows"
# and mean the same rows as the rest of the suite. Taken from test_golden.cpp.
#
# In a helper rather than in test-gradient.R, where it used to live, because
# test-gradient-batch.R needs the same points to compare against and each test
# file gets its own environment -- a second copy of these drivers would be free
# to drift, and then the two files would be pinning different operating points
# while appearing to pin the same ones.
# The grid is the multi-layer supply, so it is named: the default is the single
# potential.
grid_drivers <- function(psi_soil, ppfd = 900, vpd = 2.0, layers = 1L) {
  theta <- 0.000157
  area_leaf <- 0.05
  list(supply = leaf_supply_multilayer(),
       psi_soil = psi_soil + 0.25 * (seq_len(layers) - 1),
       PPFD = ppfd,
       soil_depth = 1.0 * seq_len(layers),
       root_network = root_network_from_carbon(
         rep(1 / layers / area_leaf, layers),
         soil_depth = 1.0 * seq_len(layers)),
       leaf_specific_conductance_max = 1.0 * theta / 5.0,
       atm_vpd = vpd, ca = 40, leaf_temp = 25, atm_o2_kpa = 21,
       atm_kpa = 101.3)
}

# TF24 on the collar route, the grid's model, unless the call names one.
grid_gradient <- function(psi_soil, ppfd = 900, vpd = 2.0, layers = 1L, ...) {
  args <- c(grid_drivers(psi_soil, ppfd, vpd, layers), list(...))
  if (is.null(args$model)) args$model <- "TF24"
  if (is.null(args$route)) args$route <- "collar"
  do.call(leaf_gradient, args)
}

# The same drivers as ONE observation of a leaf_batch().
#
# ⚠️ `psi_soil` and `soil_depth` are wrapped in list(), which is the multi-layer
# form: a plain numeric vector is N single-layer observations, exactly as in
# leaf_solve(). Getting that backwards is the easiest mistake to make with
# leaf_batch(), and it fails loudly rather than silently -- the network and the
# profile then disagree in length.
batch_drivers <- function(psi_soil, ppfd = 900, vpd = 2.0, layers = 1L) {
  d <- grid_drivers(psi_soil, ppfd, vpd, layers)
  d$psi_soil <- list(d$psi_soil)
  d$soil_depth <- list(d$soil_depth)
  d
}

# The rows gradient_golden.tsv records: tools/gradient_golden.R writes them and
# test-gradient-batch.R reads them back, both from this one list. Each case is a
# one-observation leaf_batch() (`args`) differentiated in `pars` on (`model`,
# `route`), at a prescribed `psi` where one is given, after `setup` has been
# applied to the batch's leaf. Append, never edit: a case's rows are its record.
gradient_golden_cases <- function() {
  grid <- c("vcmax_25", "stem_P50", "stem_c", "R_d_25")
  single <- function(psi_soil) {
    list(psi_soil = psi_soil, PPFD = 900, atm_vpd = 2.0,
         supply = leaf_supply_singlelayer(),
         root_network = series_resistance(1e4))
  }
  case <- function(args, pars, model = "TF24", route = "collar", psi = NULL,
                   dpsi_dtheta = NULL, setup = NULL) {
    list(args = args, pars = pars, model = model, route = route, psi = psi,
         dpsi_dtheta = dpsi_dtheta, setup = setup)
  }
  cases <- list(
    "interior-1layer" = case(batch_drivers(2.0), grid),
    "interior-5layer" = case(batch_drivers(0.5, vpd = 0.5, layers = 5L), grid),
    # A dry pin: psi_crit bounds nothing under TF24 (#183), so the leaf pins dry
    # only where it transpires all its stem can carry -- a dry soil in humid
    # air, with the stem held at the end of its curve.
    "pinned-dry-1layer" = case(batch_drivers(6.0, ppfd = 500, vpd = 0.5), grid),
    # A shut leaf: a soil past the end of the stem curve (15.44 MPa at the
    # defaults), where the stem carries nothing.
    "shutdown-1layer" = case(c(batch_drivers(16.0),
                               list(traits = leaf_traits(root_P50 = 4.0))),
                             grid),
    "single-potential" = case(single(1.5),
                              c("vcmax_25", "leaf_specific_conductance_max",
                                "resistance", "R_d_25")))

  # Every curve on both routes, in its own parameters and three every curve
  # reads. The prices leave the collar route's optimum interior: at their
  # aligned defaults CF77 and LeastCost pin it at the dry end, where the rows
  # would come from the finite-difference arm the cases above already pin.
  cp <- cost_curve_parameters()
  prices <- list(CF77 = list(CF77_lambda = 3e5),
                 TF24_floor = list(TF24_floor_lambda_o = 1.5e5),
                 LeastCost = list(LeastCost_beta = 5))
  for (cv in cost_curve_names()) {
    own <- cp$parameter[cp$curve == cv]
    for (route in c("collar", "leaf")) {
      args <- c(if (route == "collar") batch_drivers(1.0) else single(1.0),
                prices[[cv]])
      pars <- c("vcmax_25", "stem_P50", "R_d_25", own)
      cases[[paste0(cv, "-", route)]] <- case(args, pars, cv, route)
    }
  }

  # The thermal-damage layer, on: its inputs are parameters only then.
  hot <- batch_drivers(1.0)
  hot$leaf_temp <- 38
  cases[["damage-collar"]] <- case(
    hot,
    c("vcmax_25", "A_", "D_", "topt_offset_", "c_repair_maint_"),
    setup = function(l) {
      l$use_thermal_damage_ <- TRUE
      l$A_ <- 0.1
      l$D_ <- 0.05
    })
  # A prescribed point on each route, off the optimum, with the caller's own
  # dpsi/dtheta; ProfitMax's collar is the one whose normalisers the prescribed
  # evaluation must seed.
  cases[["prescribed-collar"]] <- case(batch_drivers(1.0), grid, psi = 2.0,
                                       dpsi_dtheta = c(0.1, -0.2, 0.3, 0))
  cases[["prescribed-ProfitMax-collar"]] <- case(
    batch_drivers(1.0), c("vcmax_25", "stem_P50"), "ProfitMax", "collar",
    psi = 2.0)
  cases[["prescribed-JS22-leaf"]] <- case(
    single(1.0), c("vcmax_25", "stem_P50", "JS22_gamma"), "JS22", "leaf", psi = 2.0)

  # Every curve on both routes again with a cuticle: every case above runs at
  # g_min = 0, where the residual terms are exact zeros (hazard 18).
  cuticle <- list(traits = leaf_traits(g_min = 5e-4))
  for (cv in cost_curve_names()) {
    own <- cp$parameter[cp$curve == cv]
    for (route in c("collar", "leaf")) {
      args <- c(if (route == "collar") batch_drivers(1.0) else single(1.0),
                prices[[cv]], cuticle)
      pars <- c("vcmax_25", "stem_P50", "R_d_25", "g_min", own)
      cases[[paste0(cv, "-", route, "-cuticle")]] <- case(args, pars, cv, route)
    }
  }
  cases
}

# One golden case's gradient, through the batch.
gradient_golden_run <- function(cs) {
  b <- do.call(leaf_batch, cs$args)
  if (!is.null(cs$setup)) cs$setup(b$leaf)
  leaf_gradient_batch(b, pars = cs$pars, model = cs$model, route = cs$route,
                      psi = cs$psi, dpsi_dtheta = cs$dpsi_dtheta)
}
