// -*-c++-*-
#ifndef PHYLLOPTIM_SINGLE_POTENTIAL_HPP_
#define PHYLLOPTIM_SINGLE_POTENTIAL_HPP_

#include <phylloptim/constants.hpp>
#include <phylloptim/util.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace phylloptim {

// The other supply path: one soil water potential and a constant series
// resistance to the root collar (issue #2, stage 3).
//
// This is what a user arriving from plantecophys or tealeaves has -- a single
// psi_soil and no root-mass profile -- and it is what makes the item 7a model
// comparison possible at all: Medlyn, Prentice least-cost and Cowan-Farquhar are
// all formulated against a *single* soil water potential, so comparing them
// against TF24 under a multi-layer root network compares two things at once.
//
// It implements the same four-method contract as MultiLayerRoots:
//
//   begin_solve()   -> wettest (here: the only) suction
//   uptake(...)     -> E_up and the per-layer consumption buffer
//   duptake_dpsi(.) -> dE_up/dT_collar, a positive conductance
//   n_layers()      -> how many entries of the consumption buffer it writes
//
// SIGN CONVENTION matches MultiLayerRoots and the rest of the package: every
// water potential is a positive magnitude in MPa, and there is nothing to flip.
//
// NO KINKS, so duptake_dpsi never returns the NaN that MultiLayerRoots uses to
// request a finite-difference fallback. The contract is "NaN means fall back",
// not "NaN is expected"; a constant-conductance path simply never needs it.
class SinglePotential {
public:
  // Soil water potential as a POSITIVE magnitude, -MPa.
  double psi_soil_ = util::na_value;
  // Series soil-to-collar resistance PER UNIT LEAF AREA,
  // [MPa * s * (mol H2O)^-1 * m^2 leaf]. The whole path from soil to root collar
  // collapses to this one number -- that is the point. Per unit leaf area
  // because the leaf is purely intensive: MultiLayerRoots gets the same
  // normalisation for free, since its resistances come from root carbon that is
  // itself per unit leaf area. Zero means a perfectly-conducting path (see
  // uptake).
  double resistance_ = 0.0;
  // Gravitational head to lift water to the collar, MPa. Zero by default: the
  // bare-leaf user is not thinking about rooting depth.
  double grav_head_ = 0.0;

  // psi_soil_ as a one-element vector. Leaf threads "the current soil state"
  // through E_column / find_root_psi / find_psi_stem_from_psi_root as a vector;
  // carrying one here lets a single-potential Leaf reuse that machinery
  // unchanged, which is what keeps stage 2 off the R-facing signatures. One
  // double of duplication, and begin_solve is the only writer.
  std::vector<double> psi_soil_vec_{0.0};

  void clear() {
    psi_soil_ = util::na_value;
    psi_soil_vec_.assign(1, util::na_value);
  }

  // One layer, always: the consumption buffer gets exactly one entry.
  int n_layers() const { return 1; }

  void set_soil_state(double psi_soil) { psi_soil_ = psi_soil; }

  // Per-solve entry point. Mirrors MultiLayerRoots::begin_solve: publish the
  // soil state as a vector and report the wettest suction for bracketing. With
  // one layer the "wettest" is the only one.
  double begin_solve() {
    psi_soil_vec_.assign(1, psi_soil_);
    return psi_soil_;
  }

  // E_up at a collar suction. With a constant resistance this is Ohm's law:
  //
  //   E = (T_collar - T_soil - grav) / r
  //
  // A zero resistance means the collar equilibrates with the soil instantly, so
  // there is no potential drop to draw water across and the flux is whatever the
  // demand side asks for. That is not expressible here -- this function answers
  // "how much CAN the soil supply at this collar potential" -- so a zero
  // resistance is rejected at the point it would produce an infinity rather than
  // silently returning one.
  void uptake(double T_collar, std::vector<double>& soil_consumption,
              double& E_up) const {
    uptake_from(T_collar, psi_soil_, soil_consumption, E_up);
  }

  void uptake_from(double T_collar, double T_soil,
                   std::vector<double>& soil_consumption, double& E_up) const {
    if (!std::isfinite(T_collar)) {
      util::stop("SinglePotential::uptake invalid input; T_collar=" +
                 util::to_string(T_collar));
    }
    if (!(resistance_ > 0.0)) {
      util::stop("SinglePotential::uptake needs a positive resistance_; got " +
                 util::to_string(resistance_));
    }
    const double E_i = (T_collar - T_soil - grav_head_) / resistance_;
    if (!soil_consumption.empty()) {
      soil_consumption[0] = E_i;  // mol, as MultiLayerRoots leaves it
    }
    E_up = E_i * kg_per_mol_h2o;
    if (!std::isfinite(E_up)) {
      util::stop("SinglePotential::uptake non-finite E_up; T_collar=" +
                 util::to_string(T_collar));
    }
  }

  // Vector-taking forms, so Leaf can call both supply paths identically. The
  // suction is read from element 0 rather than from psi_soil_, so that
  // the R-facing Leaf::E_from_Soil_to_Root_Collar -- which may be handed any
  // vector at all -- means the same thing here as it does for MultiLayerRoots.
  void uptake_at(double T_collar, const std::vector<double>& psi_soil,
                 std::vector<double>& soil_consumption, double& E_up) const {
    if (psi_soil.empty()) {
      util::stop("SinglePotential::uptake_at needs at least one potential");
    }
    uptake_from(T_collar, psi_soil[0], soil_consumption, E_up);
  }

  double duptake_dpsi(double T_collar,
                      const std::vector<double>& psi_soil) const {
    static_cast<void>(T_collar);
    static_cast<void>(psi_soil);
    return duptake_dpsi();
  }

  // Constant, and exactly the analytic derivative of uptake above. POSITIVE: it
  // is a conductance, and uptake rises as the collar pulls harder. Same sign
  // convention as MultiLayerRoots::duptake_dpsi, and no longer negated by
  // marginal_cost_water_multilayer (#25).
  double duptake_dpsi() const { return kg_per_mol_h2o / resistance_; }
};

}  // namespace phylloptim

#endif
