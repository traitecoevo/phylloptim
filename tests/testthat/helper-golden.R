# Helpers for the R-side golden check. See test-golden.R for why it exists.

# This used to restate the C++ default constructor's trait values, because stage
# 1 bound only the 19-argument constructor and the R side had to name them all.
# It no longer does, and that is the improvement: `leaf_model()` carries the
# package's own defaults now, so solving through it here is what CHECKS those
# defaults still equal `leaf::Leaf`'s. If the two ever drift apart the bit-exact
# comparisons in test-golden.R fail, rather than quietly measuring a different
# model than the one the golden file was generated from.
default_leaf <- function() leaf_model()

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
golden_kmax <- golden_ks * golden_theta / golden_h

# One grid point, set up exactly as test_golden.cpp's solve() does: the soil
# profile spread over `layers` equal 1 m layers drying with depth, and the root
# carbon split evenly. Returns the same nine outputs the golden file records.
#
# Deliberately calls $set_physiology() directly rather than going through
# set_drivers(). This is the tie-back to the C++ baseline, so it wants as little
# of our own R between it and the binding as possible; test-surface.R covers
# set_drivers() on its own and checks the two agree on the same point.
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
    leaf_specific_conductance_max = golden_kmax,
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
