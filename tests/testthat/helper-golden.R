# Helpers for the R-side golden check. See test-golden.R for why it exists.

# The C++ default constructor's trait values, restated. leaf::Leaf's default
# constructor is what tests/cpp/test_golden.cpp uses, and RcppR6 binds only the
# 19-argument one, so the R side has to name them. Stage 2 turns this into the
# package's own defaults; until then it lives here, and the first test in
# test-golden.R is the one that checks it still matches -- if the C++ defaults
# move and this list does not, every other test in the file starts measuring the
# wrong model and would otherwise just fail confusingly.
#
# root_c / root_b / root_psi_crit / beta_R_H / beta_R_V default in
# MultiLayerRoots rather than in Leaf: the root Weibull pair deliberately has no
# second copy in Leaf, because an unmarked duplicate of a parameter that exists
# in two versions is hazard 1 in the developer guide, and it has already cost
# this project a wrong exponent in a manuscript draft.
leaf_cpp_defaults <- list(
  vcmax_25 = 96,
  stem_c = 2.680147,
  stem_b = 3.898245,
  psi_crit = 5.870283,
  root_c = 2.680147,
  root_b = 3.898245,
  root_psi_crit = 5.870283,
  beta2 = 1.5,
  jmax_25 = 157.44,
  a = 0.30,
  curv_fact_elec_trans = 0.7,
  curv_fact_colim = 0.99,
  GSS_tol_abs = 1e-3,
  vulnerability_curve_ncontrol = 100,
  ci_abs_tol = 1e-3,
  ci_niter = 1000,
  cost_scale_TF24 = 7.5,
  beta_R_H = 3.4e2,
  beta_R_V = 9.4e3
)

default_leaf <- function() do.call(Leaf, leaf_cpp_defaults)

# The fixed drivers from tests/cpp/test_golden.cpp, which took them in turn from
# plant's tests/testthat/test-leaf.r.
golden_theta <- 0.000157
golden_ks <- 1.0
golden_h <- 5.0
golden_area_leaf <- 0.05
golden_ca <- 40.0
golden_o2 <- 21.0
golden_tleaf <- 25.0
golden_patm <- 101.3

# One grid point, set up exactly as test_golden.cpp's solve() does: the soil
# profile spread over `layers` equal 1 m layers drying with depth, and the root
# carbon split evenly. Returns the same nine outputs the golden file records.
golden_solve <- function(psi_soil, ppfd, vpd, layers) {
  l <- default_leaf()

  i <- seq_len(layers) - 1L
  ps <- psi_soil + 0.25 * i
  depth <- 1.0 * (i + 1)
  root <- rep(1.0 / layers / golden_area_leaf, layers)

  l$set_physiology(
    root_carbon_per_leaf_area = root,
    PPFD = ppfd,
    psi_soil = ps,
    soil_depth = depth,
    leaf_specific_conductance_max = golden_ks * golden_theta / golden_h,
    atm_vpd = vpd,
    ca = golden_ca,
    leaf_temp = golden_tleaf,
    atm_o2_kpa = golden_o2,
    atm_kpa = golden_patm
  )
  l$find_root_collar_psi()

  consumption <- l$soil_consumption_
  list(
    psi_stem = l$opt_psi_stem_,
    opt_root_psi = l$opt_root_psi_,
    ci = l$ci_,
    assim = l$assim_colimited_,
    transpiration = l$transpiration_,
    gc = l$stom_cond_CO2_,
    profit = l$profit_,
    e_up = l$E_up_,
    uptake = sum(consumption[is.finite(consumption)])
  )
}
