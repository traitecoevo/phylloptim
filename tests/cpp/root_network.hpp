// -*-c++-*-
// The root-architecture model the C++ suite drives the leaf with.
//
// Before #33 these two constants were `Leaf` constructor arguments with defaults,
// and the carbon -> resistance map ran inside `set_physiology`. Now the map is a
// helper the CALLER invokes (`phylloptim::root_network_from_carbon`) and the
// constants belong to whoever owns that model -- in plant, `TF24_Strategy`; here,
// this file.
//
// It exists as a shared header rather than a copy per test program for one
// reason: the golden file's bit-exactness depends on these two numbers. Four
// copies would mean a change to one of them fails `test_golden` while the two
// bench programs quietly go on measuring a different root system.
#ifndef PHYLLOPTIM_TESTS_ROOT_NETWORK_HPP_
#define PHYLLOPTIM_TESTS_ROOT_NETWORK_HPP_

#include <phylloptim.hpp>

#include <vector>

namespace fixture {

// The values that used to be `Leaf`'s beta_R_H / beta_R_V defaults, and are
// still plant's TF24 defaults.
inline constexpr double beta_R_H = 3.4e2; // MPa s (mol C) / (mol H2O)
inline constexpr double beta_R_V = 9.4e3; // MPa (mol C) s / (mol H2O) / m^2

// Carbon -> resistances, exactly as set_physiology used to do it internally.
// The suite drives the leaf from a cumulative depth profile, so it differences
// that profile into per-layer widths here (#626) rather than holding widths of
// its own. Bit-exact for the golden grid, whose layers are exactly 1 m.
inline phylloptim::RootNetwork
root_network(const std::vector<double>& root_carbon_per_leaf_area,
             const std::vector<double>& soil_depth) {
  return phylloptim::root_network_from_carbon(
      root_carbon_per_leaf_area, phylloptim::layer_thicknesses(soil_depth),
      beta_R_H, beta_R_V);
}

// The single-potential path's driver: one series resistance, and no
// vulnerability-weighted horizontal term. The SAME RootNetwork the multi-layer
// path takes -- `r_R_V_sum` already means "series resistance to the collar", so
// this is that field with one layer, which is why one set_physiology argument
// serves both paths.
inline phylloptim::RootNetwork series_resistance(double r) {
  phylloptim::RootNetwork out;
  out.r_R_V_sum.assign(1, r);
  return out;
}

// A leaf seated as plant seats it: the multi-layer supply and TF24 on the
// collar route. A bare Leaf defaults to TF24_floor on the leaf route over the
// single-potential supply.
inline phylloptim::Leaf multilayer_leaf() {
  phylloptim::Leaf l;
  l.set_supply_multilayer();
  l.set_model(phylloptim::Leaf::CostCurve::TF24, true);
  return l;
}

// The operating-point kinds the derivative surface's tests must reach on each
// curve, which their coverage assertions require. A kind is left out only where
// the curve cannot reach it at plausible drivers (measured over soil 0.3 to 16
// MPa, light 1 to 1500, deficit 0.5 to 3.1 kPa and, for JS22 and CMax, their
// level at 1e-4 of its default): a price on the drop in potential keeps the
// leaf off the stem's end; a product's gain vanishes at psi_crit, so its
// optimum never sits there; the solve refuses a dark LeastCost leaf, whose
// log link has no slope at A <= 0, rather than shutting it on its wet bound;
// and a shut ProfitMax leaf has no stream to normalise over, so no objective,
// while none pinned dry over that grid with one.
inline std::vector<phylloptim::Leaf::OperatingPointKind>
surface_kinds(phylloptim::Leaf::CostCurve curve) {
  using CC = phylloptim::Leaf::CostCurve;
  using Kind = phylloptim::Leaf::OperatingPointKind;
  switch (curve) {
    case CC::TF24:
    case CC::TF24_floor:
    case CC::CF77:
      return {Kind::Interior, Kind::BoundarySoil, Kind::BoundaryDry,
              Kind::HydraulicShutdown};
    case CC::ProfitMax:
      return {Kind::Interior, Kind::BoundarySoil};
    case CC::JS22:
    case CC::CMax:
    case CC::SOX:
    case CC::JW26:
      return {Kind::Interior, Kind::BoundarySoil, Kind::HydraulicShutdown};
    case CC::LeastCost:
      return {Kind::Interior, Kind::BoundaryDry, Kind::HydraulicShutdown};
  }
  return {};
}

// Whether the surface has an objective to differentiate at a solved point:
// false only for a normalised curve with no stream to normalise over (|A|max or
// k_span not positive), which the leaf route meets with its degenerate exit and
// the surface refuses by name.
inline bool surface_has_objective(const phylloptim::Leaf& l) {
  return l.cost_curve_ != phylloptim::Leaf::CostCurve::ProfitMax ||
         (l.profitmax_A_max() > 0.0 && l.profitmax_k_span() > 0.0);
}

}  // namespace fixture

#endif
