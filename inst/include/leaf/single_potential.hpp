// -*-c++-*-
#ifndef LEAF_SINGLE_POTENTIAL_HPP_
#define LEAF_SINGLE_POTENTIAL_HPP_

#include <leaf/constants.hpp>
#include <leaf/util.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace leaf {

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
//   begin_solve()   -> wettest (here: the only) potential, SIGNED
//   uptake(...)     -> E_up and the per-layer consumption buffer
//   duptake_dpsi(.) -> dE_up/dP_collar
//   n_layers()      -> how many entries of the consumption buffer it writes
//
// SIGN CONVENTION matches MultiLayerRoots: psi_soil_ is supplied as a positive
// magnitude and flipped once, in begin_solve(), to the signed potential the
// transport works in.
//
// NO KINKS, so duptake_dpsi never returns the NaN that MultiLayerRoots uses to
// request a finite-difference fallback. The contract is "NaN means fall back",
// not "NaN is expected"; a constant-conductance path simply never needs it.
class SinglePotential {
public:
  // Soil water potential as a POSITIVE magnitude, -MPa.
  double psi_soil_ = util::na_value;
  // Series soil-to-collar resistance, [MPa * s * (mol H2O)^-1]. The whole path
  // from soil to root collar collapses to this one number -- that is the point.
  // Zero means a perfectly-conducting path (see uptake).
  double resistance_ = 0.0;
  // Gravitational head to lift water to the collar, MPa. Zero by default: the
  // bare-leaf user is not thinking about rooting depth.
  double grav_head_ = 0.0;

  // Signed (<= 0) form of psi_soil_, built by begin_solve.
  double psi_soil_inverted_ = util::na_value;
  // The same value as a one-element vector. Leaf threads "the current soil state,
  // signed" through E_column / find_root_psi / find_psi_stem_from_psi_root as a
  // vector; carrying one here lets a single-potential Leaf reuse that machinery
  // unchanged, which is what keeps stage 2 off the R-facing signatures. One
  // double of duplication, and begin_solve is the only writer of either.
  std::vector<double> psi_soil_inverted_vec_{0.0};

  void clear() {
    psi_soil_ = util::na_value;
    psi_soil_inverted_ = util::na_value;
    psi_soil_inverted_vec_.assign(1, util::na_value);
  }

  // One layer, always: the consumption buffer gets exactly one entry.
  int n_layers() const { return 1; }

  void set_soil_state(double psi_soil) { psi_soil_ = psi_soil; }

  // Per-solve entry point. Mirrors MultiLayerRoots::begin_solve: flip to the
  // signed convention and report the wettest potential for bracketing.
  double begin_solve() {
    psi_soil_inverted_ = -psi_soil_;
    psi_soil_inverted_vec_.assign(1, psi_soil_inverted_);
    return psi_soil_inverted_;
  }

  // E_up at a collar potential. With a constant resistance this is Ohm's law:
  //
  //   E = (psi_soil - P_collar - grav) / (r * area_leaf)
  //
  // A zero resistance means the collar equilibrates with the soil instantly, so
  // there is no potential drop to draw water across and the flux is whatever the
  // demand side asks for. That is not expressible here -- this function answers
  // "how much CAN the soil supply at this collar potential" -- so a zero
  // resistance is rejected at the point it would produce an infinity rather than
  // silently returning one.
  void uptake(double P_x_r, double area_leaf,
              std::vector<double>& soil_consumption, double& E_up) const {
    uptake_from(P_x_r, psi_soil_inverted_, area_leaf, soil_consumption, E_up);
  }

  void uptake_from(double P_x_r, double psi_soil_signed, double area_leaf,
                   std::vector<double>& soil_consumption, double& E_up) const {
    if (!std::isfinite(P_x_r) || !std::isfinite(area_leaf)) {
      util::stop("SinglePotential::uptake invalid input; P_x_r=" +
                 util::to_string(P_x_r) +
                 "; area_leaf=" + util::to_string(area_leaf));
    }
    if (!(resistance_ > 0.0)) {
      util::stop("SinglePotential::uptake needs a positive resistance_; got " +
                 util::to_string(resistance_));
    }
    const double E_i =
        (psi_soil_signed - P_x_r - grav_head_) / (resistance_ * area_leaf);
    if (!soil_consumption.empty()) {
      soil_consumption[0] = E_i;  // mol, as MultiLayerRoots leaves it
    }
    E_up = E_i * kg_per_mol_h2o;
    if (!std::isfinite(E_up)) {
      util::stop("SinglePotential::uptake non-finite E_up; P_x_r=" +
                 util::to_string(P_x_r));
    }
  }

  // Vector-taking forms, so Leaf can call both supply paths identically. The
  // potential is read from element 0 rather than from psi_soil_inverted_, so that
  // the R-facing Leaf::E_from_Soil_to_Root_Collar -- which may be handed any
  // vector at all -- means the same thing here as it does for MultiLayerRoots.
  void uptake_at(double P_x_r, const std::vector<double>& psi_soil,
                 double area_leaf, std::vector<double>& soil_consumption,
                 double& E_up) const {
    if (psi_soil.empty()) {
      util::stop("SinglePotential::uptake_at needs at least one potential");
    }
    uptake_from(P_x_r, psi_soil[0], area_leaf, soil_consumption, E_up);
  }

  double duptake_dpsi(double P_x_r, const std::vector<double>& psi_soil,
                      double area_leaf) const {
    static_cast<void>(P_x_r);
    static_cast<void>(psi_soil);
    return duptake_dpsi(area_leaf);
  }

  // Constant, and exactly the analytic derivative of uptake above. Negative,
  // because uptake rises as the collar gets MORE negative -- the same sign
  // convention MultiLayerRoots::duptake_dpsi uses, and the one
  // marginal_cost_water_multilayer negates to recover a conductance.
  double duptake_dpsi(double area_leaf) const {
    return -kg_per_mol_h2o / (resistance_ * area_leaf);
  }
};

}  // namespace leaf

#endif
