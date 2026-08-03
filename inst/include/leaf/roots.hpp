// -*-c++-*-
#ifndef LEAF_ROOTS_HPP_
#define LEAF_ROOTS_HPP_

#include <leaf/constants.hpp>
#include <leaf/potential.hpp>
#include <leaf/util.hpp>
#include <leaf/vulnerability.hpp>

#include <odelia/interpolator.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <string>
#include <vector>

namespace leaf {

// The per-layer root hydraulic resistances the supply solve actually consumes.
// Two of the five fields are load-bearing; see set_root_network on why the other
// three are here at all.
struct RootNetwork {
  // Minimum (fully-hydrated) horizontal, intra-layer soil->root resistance.
  // Divided by the vulnerability-weighted mean conductivity at the operating
  // potential to get the actual horizontal resistance.
  std::vector<double> r_R_H_min;
  // Cumulative vertical, inter-layer resistance from the surface down to layer i.
  std::vector<double> r_R_V_sum;
  // Diagnostics only -- read by nothing in the model.
  std::vector<double> c_r_V, c_r_H, r_R_V;
};

// The root-architecture model: how carbon invested in roots becomes hydraulic
// resistance. Each layer's root carbon is split 1/3 vertical : 2/3 horizontal,
// and
//   network_.r_R_H_min[i] = beta_R_H / c_r_h        (min horizontal resistance, i.e. the
//                                           reciprocal of max conductance)
//   r_R_V[i]     = beta_R_V * dz^2 / c_r_v (vertical; dz^2 because vertical
//                                           conductivity scales with root
//                                           cross-sectional area)
//   network_.r_R_V_sum[i] = cumulative vertical resistance from the surface to layer i.
//
// The returned vectors are sized to the deepest layer with non-zero root carbon,
// so the hot loop only iterates over layers that actually contain roots.
//
// This is deliberately a free function rather than a member: it is the one part
// of the supply path that is a *model of the plant*, not of water transport, and
// keeping it out of MultiLayerRoots is what lets an alternative supply path
// exist without inventing a root system. It stays in this package rather than
// moving to plant so that it remains covered by the golden file -- the moment
// this arithmetic crosses the package boundary it leaves the tested surface.
//
// NOTE on zero-carbon layers: a layer with no roots currently gets
// r_R_H_min = 0, i.e. *zero* horizontal resistance, which is infinite
// soil-to-root conductance in a layer with no roots -- backwards. It appears
// unreachable from plant today (max_soil_layer truncates at the last non-zero
// layer, and plant's Q() root distribution does not produce an exact interior
// zero), so this preserves the behaviour rather than changing it silently.
//
// Fills `out` in place so its buffers are reused: this runs once per
// set_physiology, which plant calls once per solve, and building five fresh
// vectors each time measured +0.074 us/call (0.061 -> 0.135), about +2% of a
// whole solve. The value-returning overload below is for tests and one-off
// callers, where that does not matter.
inline void root_network_from_carbon(
    const std::vector<double>& root_carbon_per_layer, double dz,
    double beta_R_H, double beta_R_V, RootNetwork& out) {
  const size_t n_layers = root_carbon_per_layer.size();

  // deepest layer with non-zero root carbon
  int max_soil_layer = 0;
  for (size_t i = 0; i < n_layers; ++i) {
    if (root_carbon_per_layer[i] != 0) {
      max_soil_layer = i + 1;
    }
  }
  out.c_r_V.assign(max_soil_layer, 0.0);
  out.c_r_H.assign(max_soil_layer, 0.0);
  out.r_R_H_min.resize(max_soil_layer);
  out.r_R_V.resize(max_soil_layer);
  out.r_R_V_sum.resize(max_soil_layer);

  const double dz_sq = dz * dz;
  double vertical_resistance_sum = 0.0;
  for (int i = 0; i < max_soil_layer; ++i) {
    if(root_carbon_per_layer[i] < 0){
            util::stop("Root mass lower than 0");
    }
    const double root_mass = root_carbon_per_layer[i];
    if (root_mass == 0.0) {
      out.r_R_H_min[i] = 0.0;
      out.r_R_V[i] = 0.0;
      out.r_R_V_sum[i] = vertical_resistance_sum;
      continue;
    }

    const double c_r_v = root_mass / 3.0;
    const double c_r_h = root_mass * 2.0 / 3.0;
    out.c_r_V[i] = c_r_v;
    out.c_r_H[i] = c_r_h;

    // Set horizantal minimum resistance per soil layer (i.e. reciprocal of maximum conductance).
    out.r_R_H_min[i] = beta_R_H / c_r_h;
    // The vertical conductivity is likely linearly proportional to the root area projected onto the horizontal plane, hence dz^2.
    out.r_R_V[i] = beta_R_V * dz_sq / c_r_v;
    vertical_resistance_sum += out.r_R_V[i];
    out.r_R_V_sum[i] = vertical_resistance_sum;
  }
}

// Same, returning a fresh network. Convenience for tests and standalone callers.
inline RootNetwork root_network_from_carbon(
    const std::vector<double>& root_carbon_per_layer, double dz,
    double beta_R_H, double beta_R_V) {
  RootNetwork out;
  root_network_from_carbon(root_carbon_per_layer, dz, beta_R_H, beta_R_V, out);
  return out;
}

// ---------------------------------------------------------------------------
// SOIL -> ROOT-COLLAR WATER SUPPLY
// ---------------------------------------------------------------------------
// The water supply side of the model, lifted out of Leaf (issue #2). Everything
// here answers one question: given a collar potential, how much water can the
// root system deliver, and how fast does that change? The gas-exchange core does
// not need to know that soil has layers -- it only ever consumes
//
//   E_up  = uptake(P_collar)        and       dE_up/dP = duptake_dpsi(P_collar)
//
// so this is where a single-potential alternative will plug in (stage 3).
//
// Scientific model (after Potkay et al. 2021; prototyped in plant's
// vignettes/models/root_water_uptake.Rmd as E_from_Soil_to_Root_Collar):
//
// The root system is represented as a set of parallel soil layers, each
// connected to a single root collar (the point where roots join the stem).
// Within each layer i, water flows from soil to collar driven by the water
// potential gradient (psi_soil[i] - P_x_r), corrected for the gravitational
// head needed to lift water to the layer midpoint (gravity_head * z_soil_mid).
//
// The hydraulic resistance of each layer is the sum of two terms:
//   * r_R_H : horizontal (intra-layer, soil->root) resistance. Set during
//             set_root_network as network_.r_R_H_min[i] / f_r, where r_R_H_min scales
//             with the carbon invested in horizontal roots and f_r is the
//             fractional loss of conductivity from the root vulnerability
//             curve at the operating potential.
//   * r_R_V : vertical (inter-layer, along the root axis to the collar)
//             resistance, accumulated from the surface down to layer i
//             (r_R_V_sum). It scales with dz^2 / carbon-in-vertical-roots.
//
// Because the root vulnerability curve f_r is non-linear in psi, the
// horizontal resistance is evaluated using the *average* fractional
// conductivity over the potential interval spanned between the soil and the
// collar (P_src_min..P_src_max). This mean is obtained as
// (1/(b-a)) * integral_a^b f_r dpsi from a pre-integrated curve
// (root_vuln_integral_from_psi) with two spline evals, the same technique used
// for stem transpiration in Leaf::setup_transpiration.
//
// SIGN CONVENTION, NOW IN THE TYPE (issue #8): everything in this class works in
// SIGNED potentials, spelled `Psi` (<= 0, MPa). psi_soil_ arrives from the
// caller as positive magnitudes -- it is plant-bound as std::vector<double>, so it
// stays that -- and is flipped once, in begin_solve(), into psi_soil_inverted_,
// which is std::vector<Psi>. The vulnerability splines take a magnitude, so
// those sites ask for one by name: root_vuln_from_psi.eval(P_src_min.abs()
// .value) rather than the old unmarked `eval(-P_src_min)`. See leaf/potential.hpp
// for why, and the sign map in leaf_model.hpp for the whole picture.
//
// TWO OUTPUTS ARE NOT OWNED HERE. E_up and the per-layer soil_consumption
// buffer are passed in by reference rather than stored, because plant reaches
// into `leaf.E_up_` / `leaf.soil_consumption_[a]` by name and *writes back*
// into them after crown integration (tf24_strategy.cpp:501-508). They are
// plant's buffers, not this object's state. Note the deliberate unit split:
// E_up is kg H2O m^-2 s^-1, soil_consumption[i] is mol, converted downstream.
class MultiLayerRoots {
public:
  // --- root vulnerability trait pair (hazard 1: NOT the stem's b/c) ---------
  double root_c = 2.680147;       // unitless
  // The characteristic TENSION for the ROOT curve (hazard 1: not the stem's).
  AbsPsi root_b{3.898245};        // -MPa
  // The driest collar potential this supply path can be asked about. A published
  // trait, quoted as a magnitude, so it is a AbsPsi -- which is what stops it
  // being compared against a signed potential (plant #584; see potential.hpp).
  AbsPsi root_psi_crit{5.870283};

  // proportionality constant between minimum horizontal (intralayer) root
  // hydraulic resistance and C_r^-1, [MPa * s * (mol C) / (mol H2O)]
  double beta_R_H = 3.4e2;
  // proportionality constant between minimum vertical (interlayer) root
  // hydraulic resistance and dz^2/C_r, [MPa * (mol C) * s / (mol H2O) / m^2]
  double beta_R_V = 9.4e3;

  // pre-computed root vulnerability curve f_r(m) = exp(-(m/root_b)^root_c)
  odelia::interpolator::Interpolator root_vuln_from_psi;
  // cumulative integral of it, G(m) = int_0^m f_r(s) ds, indexed by magnitude
  // m = -psi. Lets uptake() obtain the mean conductivity over a potential
  // interval from 2 evals instead of (n+1).
  odelia::interpolator::Interpolator root_vuln_integral_from_psi;

  // --- soil geometry -------------------------------------------------------
  // The four scalars carry the same unset sentinels clear() assigns, so a bare
  // MultiLayerRoots is never indeterminate. Leaf reaches them only after
  // set_soil_state / set_root_network, but the class is public now.
  double soil_number_of_depths_ = util::na_value_int;
  int max_soil_layer = util::na_value_int;  // deepest layer with non-zero root mass
  std::vector<double> soil_depth_;
  std::vector<double> z_soil_mid_;
  // Per-layer gravitational head gravity_head * z_soil_mid_[i], precomputed once
  // per set_soil_state (z_soil_mid_ is fixed across a collar solve). Used three
  // times per layer in uptake()'s hot loop; caching it removes a redundant
  // multiply per layer per (re)evaluation.
  std::vector<double> grav_head_z_;
  bool use_precomputed_z_soil_mid_ = false;
  double dz_ = util::na_value;

  // --- soil state ----------------------------------------------------------
  // Positive magnitudes, as supplied. Stays std::vector<double> rather than
  // becoming std::vector<AbsPsi>: plant's RcppR6 binds this by name as
  // `roots_.psi_soil_` with type std::vector<double> (hazard 7), and it is an
  // input that is never mixed with a signed potential -- begin_solve is its only
  // reader.
  std::vector<double> psi_soil_;
  // The same state in the signed convention, built by begin_solve. This is the
  // member the `_inverted` suffix used to be the only guard on.
  std::vector<Psi> psi_soil_inverted_;
  // Per-layer cache of root_vuln_integral_from_psi.eval(psi_soil_inverted_[i].abs().value).
  // psi_soil_inverted_ is fixed for the whole collar solve, so the soil-side
  // endpoint of the cumulative-integral lookup in uptake() is constant across
  // every (re)evaluation of the nested root-finders. Precomputing it once per
  // solve (alongside the P_x_r-side eval, hoisted out of the layer loop)
  // collapses ~2 spline evals per layer to ~1 per call. Rebuilt in begin_solve.
  std::vector<double> root_vuln_integral_soil_;

  // --- root resistance network --------------------------------------------
  // Held as one object rather than five loose vectors so the carbon -> resistance
  // map can fill it in place and reuse its buffers. The solve reads exactly two
  // of its fields.
  RootNetwork network_;

  // -------------------------------------------------------------------------

  // Reset every state member to the unset sentinel. Mirrors Leaf::setup_clean_leaf.
  void clear() {
    psi_soil_.clear();
    soil_depth_.clear();
    z_soil_mid_.clear();
    grav_head_z_.clear();
    use_precomputed_z_soil_mid_ = false;
    network_.c_r_V.clear();
    network_.c_r_H.clear();
    network_.r_R_H_min.clear();
    network_.r_R_V.clear();
    network_.r_R_V_sum.clear();
    soil_number_of_depths_ = util::na_value_int;
    max_soil_layer = util::na_value_int;
  }

  // Pre-compute the root vulnerability curve and its cumulative integral over
  // [0, psi_max_root], the range where conductivity drops to 1%. Avoids repeated
  // exp(pow(...)) inside uptake().
  void setup_vulnerability(double resolution) {
    std::vector<double> x_psi_root, y_integral;
    cumulative_vulnerability_integral(root_b.value, root_c, resolution, x_psi_root,
                                      y_integral);

    // f_r conductivity knots on the same grid. f_r(0) = exp(-pow(0,root_c)) = 1.
    std::vector<double> y_f_r(x_psi_root.size());
    for (size_t i = 0; i < x_psi_root.size(); ++i) {
      y_f_r[i] = exp(-pow(AbsPsi{x_psi_root[i]} / root_b, root_c));
    }
    root_vuln_from_psi.init(x_psi_root, y_f_r);
    root_vuln_from_psi.set_extrapolate(true); // clamp to last value beyond range

    root_vuln_integral_from_psi.init(x_psi_root, y_integral);
    // linear extrapolation beyond range: slope ~= f_r at the tail (~1%), so the
    // integral keeps growing consistently with the clamped-conductivity tail.
    root_vuln_integral_from_psi.set_extrapolate(true);
  }

  // Per-timestep soil state: the layer potentials, the layer depths, and the
  // gravitational head that follows from them.
  void set_soil_state(const std::vector<double>& psi_soil,
                      const std::vector<double>& soil_depth) {
    psi_soil_ = psi_soil;
    soil_depth_ = soil_depth;
    soil_number_of_depths_ = soil_depth_.size();

    if (!(use_precomputed_z_soil_mid_ &&
          z_soil_mid_.size() == static_cast<size_t>(soil_number_of_depths_))) {
      // Fallback for paths that do not provide environment-precomputed midpoints.
      z_soil_mid_.resize(soil_number_of_depths_);
      for (size_t i = 0; i < soil_number_of_depths_; ++i) {
        if (i == 0) {
          z_soil_mid_[i] = (soil_depth_[i] / 2.0);
        } else {
          z_soil_mid_[i] = ((soil_depth_[i - 1] + soil_depth_[i]) / 2.0);
        }
      }
    }

    use_precomputed_z_soil_mid_ = false;

    grav_head_z_.resize(soil_number_of_depths_);
    for (size_t i = 0; i < soil_number_of_depths_; ++i) {
      grav_head_z_[i] = gravity_head * z_soil_mid_[i];
    }

    // Layer thickness is soil geometry, not root architecture, so it is set here
    // rather than alongside the resistance network that consumes it.
    dz_ = soil_depth_.back()/soil_number_of_depths_;
  }

  // Per-timestep root resistance network. Takes the resistances themselves, not
  // the root carbon they are derived from.
  //
  // WHY THIS TAKES RESISTANCES. The solve reads exactly two of these vectors --
  // r_R_H_min and r_R_V_sum -- plus grav_head_z_ and max_soil_layer. Nothing in
  // uptake() or duptake_dpsi() touches root carbon, the 1/3 : 2/3 split, dz, or
  // either beta_R_* constant; those are inputs to a *root architecture* model
  // that happens to run just before. Splitting them out is the same move
  // leaf_specific_conductance_max already makes: plant computes
  // kmax = K_s*theta/(h*eta_c) and hands over a scalar, so which
  // conductance-versus-height model is in force is not this package's business.
  // The carbon -> resistance map lives in root_network_from_carbon below, so an
  // alternative supply path can supply resistances any way it likes.
  //
  // c_r_V_, c_r_H_ and r_R_V ride along as diagnostics: nothing in the model
  // reads them, but plant exposes them through RcppR6, so they are carried
  // rather than dropped. They are removal candidates with item 6.
  // Takes by value and moves: the caller always passes a freshly-built
  // RootNetwork, so this hands over its buffers rather than copying five vectors
  // on every set_physiology.
  void set_root_network(RootNetwork network) {
    network_ = std::move(network);
    max_soil_layer = static_cast<int>(network_.r_R_H_min.size());
  }

  // The carbon-based caller. Kept as an entry point (rather than making everyone
  // go through root_network_from_carbon) so the fill happens straight into
  // network_ and its buffers are reused across calls.
  void set_root_network_from_carbon(const std::vector<double>& root_carbon_per_layer) {
    root_network_from_carbon(root_carbon_per_layer, dz_, beta_R_H, beta_R_V,
                             network_);
    max_soil_layer = static_cast<int>(network_.r_R_H_min.size());
  }

  // Per-solve entry point. Flips psi_soil_ into the signed convention, builds
  // the soil-side cumulative-integral cache, and returns the wettest (least
  // negative) layer potential, which is the bracket endpoint the collar solve
  // needs. Both caches are valid until the next call.
  //
  // This is a single pass on purpose: the cache is the measured hot-path
  // optimisation described on root_vuln_integral_soil_, and the wettest layer
  // falls out of the same loop.
  Psi begin_solve() {
    psi_soil_inverted_.resize(max_soil_layer);
    root_vuln_integral_soil_.resize(max_soil_layer);
    Psi wettest_soil_layer{-std::numeric_limits<double>::infinity()};
    for (int i = 0; i < max_soil_layer; ++i) {
      // The one flip. psi_soil_ is supplied as a magnitude, so say so and ask for
      // the signed form by name.
      const Psi psi_inverted = AbsPsi{psi_soil_[i]}.as_psi();
      psi_soil_inverted_[i] = psi_inverted;
      root_vuln_integral_soil_[i] =
          root_vuln_integral_from_psi.eval(psi_inverted.abs().value);
      wettest_soil_layer = std::max(wettest_soil_layer, psi_inverted);
    }
    return wettest_soil_layer;
  }

  // Uptake at a collar potential, against the soil state begin_solve() cached.
  // This is the hot path: ~10^3 calls per collar solve.
  void uptake(Psi P_x_r, std::vector<double>& soil_consumption,
              double& E_up) const {
    uptake_impl(P_x_r, psi_soil_inverted_,
                root_vuln_integral_soil_.size() ==
                    static_cast<size_t>(max_soil_layer),
                soil_consumption, E_up);
  }

  // Uptake against an arbitrary vector of layer potentials, for callers that
  // want to probe the supply function away from the current soil state (the
  // R-facing Leaf::E_from_Soil_to_Root_Collar).
  //
  // The `&psi_soil == &psi_soil_inverted_` test is what used to select the
  // cached path for every caller, including the hot one. It is kept here only so
  // this entry point cannot change behaviour for a caller that happens to hand
  // back psi_soil_inverted_ itself; the hot path above no longer depends on
  // address identity to be fast. PLAN 7b-ii trap 3.
  void uptake_at(Psi P_x_r, const std::vector<Psi>& psi_soil,
                 std::vector<double>& soil_consumption, double& E_up) const {
    uptake_impl(P_x_r, psi_soil,
                (&psi_soil == &psi_soil_inverted_) &&
                    root_vuln_integral_soil_.size() ==
                        static_cast<size_t>(max_soil_layer),
                soil_consumption, E_up);
  }

  // Analytic d(E_up)/d(P_x_r): the signed-collar-potential derivative of the
  // uptake, mirroring the general branch of uptake_impl. Per layer, with
  // span = |psi_soil[i] - P_x_r| and integral = \int f_r over
  // [P_src_min, P_src_max] (root_vuln_integral_from_psi, whose integrand is
  // root_vuln_from_psi):
  //   E_i        = (psi_soil[i] - P_x_r - grav) / r_R,
  //   r_R        = r_R_H_min[i] * span / integral + r_R_V_sum[i],
  //   dspan/dP   = sign_var   (+1 if P_x_r is the upper bound, else -1),
  //   dinteg/dP  = sign_var * f_r(-P_x_r)  for P_x_r<0  (else sign_var, f_r==1),
  // and dE_i/dP follows by the quotient rule.
  //
  // CONTRACT: returns NaN when any layer sits on a branch kink (P_x_r ==
  // psi_soil[i], the gravity-balance point, or P_x_r == 0). That is deliberate,
  // not a failure -- the analytic general-branch derivative is not valid across
  // those, and the caller falls back to a central difference. An implementation
  // that threw, or returned 0, would silently degrade TF24f's acclimation
  // gradient. Any alternative supply path must keep this contract.
  double duptake_dpsi(Psi P_x_r,
                      const std::vector<Psi>& psi_soil) const {
    const double kink_tol = 1e-8;
    const Psi zero{0.0};
    double dEup_dr_mol = 0.0;

    for (int i = 0; i < max_soil_layer; i++) {
      if (std::abs(P_x_r - psi_soil[i]) < kink_tol ||
          std::abs((psi_soil[i] - P_x_r) - grav_head_z_[i]) < kink_tol ||
          std::abs(P_x_r.value) < kink_tol) {
        return std::numeric_limits<double>::quiet_NaN();
      }

      const Psi P_src_min = std::min(psi_soil[i], P_x_r);
      const Psi P_src_max = std::max(psi_soil[i], P_x_r);
      const double span = P_src_max - P_src_min;
      const double sign_var = (P_x_r > psi_soil[i]) ? 1.0 : -1.0;  // = dspan/dP_x_r

      // integral, replicated bit-for-bit from uptake_impl.
      const Psi hi_neg = std::min(P_src_max, zero);
      const Psi lo_pos = std::max(P_src_min, zero);
      double integral = 0.0;
      if (hi_neg > P_src_min) {
        integral += root_vuln_integral_from_psi.eval(P_src_min.abs().value) -
                    root_vuln_integral_from_psi.eval(hi_neg.abs().value);
      }
      if (P_src_max > lo_pos) {
        integral += (P_src_max - lo_pos);
      }

      // d(integral)/d(P_x_r): for P_x_r<0 the moving bound is in the vulnerable
      // region. The integrand is the derivative of the *same* cumulative spline
      // that produced `integral` (root_vuln_integral_from_psi.deriv), NOT the
      // separate root_vuln_from_psi spline: the two agree on the knot domain but
      // extrapolate independently (both clamp-to-last-value, #527), so beyond the
      // domain only the integral spline's own derivative stays consistent with its
      // value. For P_x_r>0 the moving bound is in the above-atmospheric part
      // (f_r==1), contributed linearly, so the slope is 1.
      const double fr_at = (P_x_r < zero)
                               ? root_vuln_integral_from_psi.deriv(
                                     P_x_r.abs().value)
                               : 1.0;
      const double dinteg_dr = sign_var * fr_at;

      const double r_R_H = network_.r_R_H_min[i] * span / integral;
      const double r_R = r_R_H + network_.r_R_V_sum[i];
      const double dr_R_H_dr =
          network_.r_R_H_min[i] * (sign_var * integral - span * dinteg_dr) / (integral * integral);
      const double dr_R_dr = dr_R_H_dr;

      const double num = (psi_soil[i] - P_x_r) - grav_head_z_[i];
      const double dnum_dr = -1.0;
      // E_i = num / r_R  ->  quotient rule.
      dEup_dr_mol += (dnum_dr * r_R - num * dr_R_dr) / (r_R * r_R);
    }

    return dEup_dr_mol * kg_per_mol_h2o;  // match E_up's kg units
  }

private:
  // Total water drawn from all layers to the collar. Writes E_up (kg H2O m^-2
  // leaf s^-1) and soil_consumption[i] (mol H2O m^-2 leaf s^-1, note the unit
  // split); a negative E_i in a layer means that layer is *gaining* water
  // (hydraulic redistribution).
  //
  // Implementation decisions:
  //   * f_r and its running integral are read from pre-computed splines
  //     (root_vuln_from_psi, root_vuln_integral_from_psi) instead of repeatedly
  //     evaluating exp(-(psi/b)^c); see setup_vulnerability.
  //   * Two special cases are handled exactly to avoid division/round-off
  //     issues: (a) collar potential equals layer potential, and (b) the
  //     gradient exactly balances gravity (E_i = 0).
  //   * The isfinite() guards are present because this is called from within
  //     nested root-finders where bad brackets can produce NaNs; they fail fast
  //     with diagnostic context rather than propagating NaN.
  void uptake_impl(Psi P_x_r, const std::vector<Psi>& psi_soil,
                   bool use_integral_cache,
                   std::vector<double>& soil_consumption, double& E_up) const {

    if (!std::isfinite(P_x_r.value)) {
      util::stop("E_from_Soil_to_Root_Collar invalid input; P_x_r=" + util::to_string(P_x_r.value));
    }

    E_up = 0;

    // Cumulative-integral spline caching (bit-identical fast path). The only two
    // arguments ever passed to root_vuln_integral_from_psi in the loop below are
    // -P_src_min and -hi_neg, each of which resolves to exactly one of
    // {-psi_soil[i], -P_x_r, 0}. -P_x_r is constant across all layers (compute
    // once), and -psi_soil[i] is constant across the whole solve (precomputed in
    // begin_solve).
    const double neg_P_x_r = P_x_r.abs().value;
    const double G_at_P_x_r =
        use_integral_cache ? root_vuln_integral_from_psi.eval(neg_P_x_r) : 0.0;

    // GUARD POLICY (the per-layer isfinite/stop guards here were added while
    // debugging the #485 drought-NaN, now fixed at source by the soil residual-
    // moisture floor). Most were defensive and redundant, so they have been
    // removed from this hot loop; the remaining two are load-bearing:
    //   * the equal-potentials f_ri <= 0 check below: root_vuln_from_psi
    //     LINEARLY extrapolates NEGATIVE beyond its domain, so a deep-drought
    //     layer can produce negative conductivity -> negative-but-FINITE r_R ->
    //     wrong-sign E_i that the post-loop isfinite(E_up) net would NOT catch.
    //   * the post-loop isfinite(E_up) check: any non-finite produced anywhere
    //     in the loop propagates into the sum and is caught there once per call.
    // Everything else is provably safe to drop on the valid path: psi_soil is
    // validated in Leaf::set_physiology; P_src_min<=P_src_max by construction;
    // the general-branch integral comes from a monotone-increasing spline so it
    // is strictly > 0 (span>0), giving r_R>0 and finite E_i; and any stray
    // NaN/Inf still reaches the post-loop net.
    const Psi zero{0.0};

    for(int i = 0; i < max_soil_layer; i++){

    // Find the most negative soil potential out of the given soil layer and the root collar
    Psi P_src_min = std::min(psi_soil[i], P_x_r);

    // Find the least negative soil potential out of the given soil layer and the root collar
    Psi P_src_max = std::max(psi_soil[i], P_x_r);

     // If root collar soil water potential equals the soil water potential in a given layer
    if(std::abs(P_x_r - psi_soil[i]) < 1e-8){

      // Fraction of conductance in roots in a given layer at most negative soil water potential (but actually is equal to root collar)
      // root_vuln_from_psi is a pre-built spline of exp(-(|psi|/b_root)^c_root), so
      // it wants the magnitude, and now has to be asked for it by name.
      double f_ri = root_vuln_from_psi.eval(P_src_min.abs().value);
      if (!std::isfinite(f_ri) || f_ri <= 0.0) {
        util::stop("E_from_Soil_to_Root_Collar invalid f_ri; layer=" + std::to_string(i) +
                   "; f_ri=" + util::to_string(f_ri) +
                   "; P_src_min=" + util::to_string(P_src_min.value) +
                   "; P_x_r=" + util::to_string(P_x_r.value));
      }

      // Fraction of conductance in roots in a given layer at most negative soil water potential
      double r_R_H = network_.r_R_H_min[i] / f_ri; // [MPa * s * (mol H2O)^-1]

      // Total root resistance (horizantal plus vertical)
      double r_R = r_R_H + network_.r_R_V_sum[i];

      // Transpiration is equivalent to gravitational water loss (i.e. layer gains water)
      double E_i = -grav_head_z_[i] / r_R ;

      soil_consumption[i] = E_i;
      E_up += E_i;

    }
    else if(std::abs((psi_soil[i] - P_x_r) - grav_head_z_[i]) < 1e-8){
      // If pressure difference perfectly balances gravity transpiration is equal to zero
      double E_i = 0.0; // [mol H2O / m^2 / s]

      soil_consumption[i] = E_i;

      E_up += E_i;

    } else{

      // Mean fractional root conductivity over the potential interval
      // [P_src_min, P_src_max], i.e. (1/(b-a)) * integral_a^b f_r dpsi.
      // Computed from the pre-integrated curve G(m) = integral_0^m f_r(s) ds
      // (root_vuln_integral_from_psi, indexed by magnitude m = -psi) with 2
      // evals instead of the old (n+1)-point sample mean. The interval is split
      // at psi = 0: for psi > 0 (above-atmospheric) vulnerability is 1.
      Psi hi_neg = std::min(P_src_max, zero); // boundary of the psi<=0 part
      Psi lo_pos = std::max(P_src_min, zero); // boundary of the psi>0 part

      // Memoised cumulative-integral lookup. Returns the exact same double the
      // spline would (same input -> same output); the comparisons select the
      // precomputed value because -P_src_min / -hi_neg are bit-for-bit equal to
      // one of the cached arguments in the common (psi<=0) case.
      const double neg_psi_soil_i = psi_soil[i].abs().value;
      auto G_integral = [&](double arg) -> double {
        if (use_integral_cache) {
          if (arg == neg_P_x_r) return G_at_P_x_r;
          if (arg == neg_psi_soil_i) return root_vuln_integral_soil_[i];
        }
        return root_vuln_integral_from_psi.eval(arg);
      };

      double integral = 0.0;
      if (hi_neg > P_src_min) {
        // psi<=0 part: magnitude m runs from hi_neg's magnitude up to P_src_min's
        integral += G_integral(P_src_min.abs().value) -
                    G_integral(hi_neg.abs().value);
      }
      if (P_src_max > lo_pos) {
        // psi>0 part: f_r == 1 over its length
        integral += (P_src_max - lo_pos);
      }

    // span = P_src_max - P_src_min > 0 here (the equal-potentials case is
    // handled in the branch above). integral comes from the monotone-increasing
    // cumulative-vulnerability spline so it is strictly > 0 over a span>0
    // interval; forming r_R_H as r_R_H_min * span / integral is one division
    // (vs the old f_r_average = integral/span then r_R_H_min/f_r_average two),
    // and needs no per-layer finiteness guard (any stray NaN/Inf propagates to
    // the post-loop isfinite(E_up) net).
    const double span = P_src_max - P_src_min;

    // Find the horizantal resistance in a given layer by dividing the minimum resistance (i.e. maximum conductivity) by the fractional loss of conductivity
    double r_R_H = network_.r_R_H_min[i] * span / integral; // [MPa * s * (mol H2O)^-1]

    // Find the total resistance in a given layer by adding the vertical resistance in that layer
    double r_R = r_R_H + network_.r_R_V_sum[i]; // [MPa * s * (mol H2O)^-1]

    // Transpiration is equal to the potentail gradient between the root collar and the soil, accounting for gravitational potential
    double E_i = ((psi_soil[i] - P_x_r) - grav_head_z_[i]) / r_R; // [mol H2O / m^2 / s]

    soil_consumption[i] = E_i;
    E_up += E_i;

    }
  }
  // Convert the summed uptake to kg H2O m^-2 s^-1, consistent with the rest of
  // the leaf model and environment. NOTE (review #10): only the aggregate E_up
  // is converted to kg here; the per-layer soil_consumption[i] above is left in
  // mol H2O m^-2 s^-1 and converted downstream in TF24_Strategy::compute_rates.
  // The two siblings therefore carry different units by design.
  E_up = E_up * kg_per_mol_h2o;
  if (!std::isfinite(E_up)) {
    util::stop("E_from_Soil_to_Root_Collar non-finite E_up_; P_x_r=" + util::to_string(P_x_r.value) +
               "; max_soil_layer=" + std::to_string(max_soil_layer));
  }
  }
};

}  // namespace leaf

#endif
