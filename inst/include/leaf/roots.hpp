// -*-c++-*-
#ifndef LEAF_ROOTS_HPP_
#define LEAF_ROOTS_HPP_

#include <leaf/constants.hpp>
#include <leaf/util.hpp>
#include <leaf/vulnerability.hpp>

#include <odelia/interpolator.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace leaf {

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
//             set_root_network as r_R_H_min[i] / f_r, where r_R_H_min scales
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
// SIGN CONVENTION: everything in this class works in SIGNED (negative)
// potentials. psi_soil_ arrives from the caller as positive magnitudes and is
// flipped once, in begin_solve(), into psi_soil_inverted_ (<= 0). The
// vulnerability splines take a magnitude, so those sites flip back with a
// leading `-` (e.g. root_vuln_from_psi.eval(-P_src_min)). See the sign map in
// leaf_model.hpp.
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
  double root_b = 3.898245;       // -MPa
  double root_psi_crit = 5.870283; // -MPa

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
  std::vector<double> psi_soil_;           // positive magnitudes, as supplied
  std::vector<double> psi_soil_inverted_;  // signed (<= 0), built by begin_solve
  // Per-layer cache of root_vuln_integral_from_psi.eval(-psi_soil_inverted_[i]).
  // psi_soil_inverted_ is fixed for the whole collar solve, so the soil-side
  // endpoint of the cumulative-integral lookup in uptake() is constant across
  // every (re)evaluation of the nested root-finders. Precomputing it once per
  // solve (alongside the P_x_r-side eval, hoisted out of the layer loop)
  // collapses ~2 spline evals per layer to ~1 per call. Rebuilt in begin_solve.
  std::vector<double> root_vuln_integral_soil_;

  // --- root resistance network --------------------------------------------
  std::vector<double> c_r_V_;  // carbon per layer in vertical transport (kg m^-2)
  std::vector<double> c_r_H_;  // carbon per layer in horizontal transport (kg m^-2)
  std::vector<double> r_R_H_min;  // min horizontal resistance per layer
  std::vector<double> r_R_V;      // vertical resistance per layer
  std::vector<double> r_R_V_sum;  // cumulative vertical resistance with depth

  // -------------------------------------------------------------------------

  // Reset every state member to the unset sentinel. Mirrors Leaf::setup_clean_leaf.
  void clear() {
    psi_soil_.clear();
    soil_depth_.clear();
    z_soil_mid_.clear();
    grav_head_z_.clear();
    use_precomputed_z_soil_mid_ = false;
    c_r_V_.clear();
    c_r_H_.clear();
    r_R_H_min.clear();
    r_R_V.clear();
    r_R_V_sum.clear();
    soil_number_of_depths_ = util::na_value_int;
    max_soil_layer = util::na_value_int;
  }

  // Pre-compute the root vulnerability curve and its cumulative integral over
  // [0, psi_max_root], the range where conductivity drops to 1%. Avoids repeated
  // exp(pow(...)) inside uptake().
  void setup_vulnerability(double resolution) {
    std::vector<double> x_psi_root, y_integral;
    cumulative_vulnerability_integral(root_b, root_c, resolution, x_psi_root,
                                      y_integral);

    // f_r conductivity knots on the same grid. f_r(0) = exp(-pow(0,root_c)) = 1.
    std::vector<double> y_f_r(x_psi_root.size());
    for (size_t i = 0; i < x_psi_root.size(); ++i) {
      y_f_r[i] = exp(-pow(x_psi_root[i] / root_b, root_c));
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
  }

  // Per-timestep root resistance network. Each layer's root carbon
  // (mass_root_prop[i], kg) is split 1/3 vertical : 2/3 horizontal (c_r_V_,
  // c_r_H_). From these:
  //   r_R_H_min[i] = beta_R_H / c_r_h        (min horizontal resistance, i.e.
  //                                           reciprocal of max conductance)
  //   r_R_V[i]     = beta_R_V * dz^2 / c_r_v (vertical; dz^2 because vertical
  //                                           conductivity scales with root
  //                                           cross-sectional area)
  //   r_R_V_sum[i] = cumulative vertical resistance from surface to layer i.
  // max_soil_layer is the deepest layer with non-zero root mass; the resistance
  // vectors are sized to it so uptake()'s hot loop only iterates over layers
  // that actually contain roots.
  //
  // Requires set_soil_state to have run this step (reads soil_depth_ and
  // soil_number_of_depths_).
  void set_root_network(const std::vector<double>& mass_root_prop) {
    dz_ = soil_depth_.back()/soil_number_of_depths_;

    // find max soil layer as last iteration with mass_root_prop greater than 0
    max_soil_layer = 0;
    for (size_t i = 0; i < soil_number_of_depths_; ++i) {
      if (mass_root_prop[i] != 0) {
        max_soil_layer = i + 1;
      }
    }
    c_r_V_.assign(max_soil_layer, 0.0);
    c_r_H_.assign(max_soil_layer, 0.0);
    r_R_H_min.resize(max_soil_layer);
    r_R_V.resize(max_soil_layer);
    r_R_V_sum.resize(max_soil_layer);

    const double dz_sq = dz_ * dz_;
    double vertical_resistance_sum = 0.0;
    for (int i = 0; i < max_soil_layer; ++i) {
      if(mass_root_prop[i] < 0){
              util::stop("Root mass lower than 0");
      }
      const double root_mass = mass_root_prop[i];
      if (root_mass == 0.0) {
        r_R_H_min[i] = 0.0;
        r_R_V[i] = 0.0;
        r_R_V_sum[i] = vertical_resistance_sum;
        continue;
      }

      const double c_r_v = root_mass / 3.0;
      const double c_r_h = root_mass * 2.0 / 3.0;
      c_r_V_[i] = c_r_v;
      c_r_H_[i] = c_r_h;

      // Set horizantal minimum resistance per soil layer (i.e. reciprocal of maximum conductance).
      r_R_H_min[i] = beta_R_H / c_r_h;
      // The vertical conductivity is likely linearly proportional to the root area projected onto the horizontal plane, hence dz^2.
      r_R_V[i] = beta_R_V * dz_sq / c_r_v;
      vertical_resistance_sum += r_R_V[i];
      r_R_V_sum[i] = vertical_resistance_sum;
    }
  }

  // Per-solve entry point. Flips psi_soil_ into the signed convention, builds
  // the soil-side cumulative-integral cache, and returns the wettest (least
  // negative) layer potential, which is the bracket endpoint the collar solve
  // needs. Both caches are valid until the next call.
  //
  // This is a single pass on purpose: the cache is the measured hot-path
  // optimisation described on root_vuln_integral_soil_, and the wettest layer
  // falls out of the same loop.
  double begin_solve() {
    psi_soil_inverted_.resize(max_soil_layer);
    root_vuln_integral_soil_.resize(max_soil_layer);
    double wettest_soil_layer = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < max_soil_layer; ++i) {
      const double psi_inverted = -psi_soil_[i];
      psi_soil_inverted_[i] = psi_inverted;
      root_vuln_integral_soil_[i] =
          root_vuln_integral_from_psi.eval(-psi_inverted);
      wettest_soil_layer = std::max(wettest_soil_layer, psi_inverted);
    }
    return wettest_soil_layer;
  }

  // Uptake at a collar potential, against the soil state begin_solve() cached.
  // This is the hot path: ~10^3 calls per collar solve.
  void uptake(double P_x_r, double area_leaf,
              std::vector<double>& soil_consumption, double& E_up) const {
    uptake_impl(P_x_r, psi_soil_inverted_,
                root_vuln_integral_soil_.size() ==
                    static_cast<size_t>(max_soil_layer),
                area_leaf, soil_consumption, E_up);
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
  void uptake_at(double P_x_r, const std::vector<double>& psi_soil,
                 double area_leaf, std::vector<double>& soil_consumption,
                 double& E_up) const {
    uptake_impl(P_x_r, psi_soil,
                (&psi_soil == &psi_soil_inverted_) &&
                    root_vuln_integral_soil_.size() ==
                        static_cast<size_t>(max_soil_layer),
                area_leaf, soil_consumption, E_up);
  }

  // Analytic d(E_up)/d(P_x_r): the signed-collar-potential derivative of the
  // uptake, mirroring the general branch of uptake_impl. Per layer, with
  // span = |psi_soil[i] - P_x_r| and integral = \int f_r over
  // [P_src_min, P_src_max] (root_vuln_integral_from_psi, whose integrand is
  // root_vuln_from_psi):
  //   E_i        = (psi_soil[i] - P_x_r - grav) / area_leaf / r_R,
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
  double duptake_dpsi(double P_x_r, const std::vector<double>& psi_soil,
                      double area_leaf) const {
    const double inv_area_leaf = 1.0 / area_leaf;
    const double kink_tol = 1e-8;
    double dEup_dr_mol = 0.0;

    for (int i = 0; i < max_soil_layer; i++) {
      if (std::abs(P_x_r - psi_soil[i]) < kink_tol ||
          std::abs((psi_soil[i] - P_x_r) - grav_head_z_[i]) < kink_tol ||
          std::abs(P_x_r) < kink_tol) {
        return std::numeric_limits<double>::quiet_NaN();
      }

      const double P_src_min = std::min(psi_soil[i], P_x_r);
      const double P_src_max = std::max(psi_soil[i], P_x_r);
      const double span = P_src_max - P_src_min;
      const double sign_var = (P_x_r > psi_soil[i]) ? 1.0 : -1.0;  // = dspan/dP_x_r

      // integral, replicated bit-for-bit from uptake_impl.
      const double hi_neg = std::min(P_src_max, 0.0);
      const double lo_pos = std::max(P_src_min, 0.0);
      double integral = 0.0;
      if (hi_neg > P_src_min) {
        integral += root_vuln_integral_from_psi.eval(-P_src_min) -
                    root_vuln_integral_from_psi.eval(-hi_neg);
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
      const double fr_at =
          (P_x_r < 0.0) ? root_vuln_integral_from_psi.deriv(-P_x_r) : 1.0;
      const double dinteg_dr = sign_var * fr_at;

      const double r_R_H = r_R_H_min[i] * span / integral;
      const double r_R = r_R_H + r_R_V_sum[i];
      const double dr_R_H_dr =
          r_R_H_min[i] * (sign_var * integral - span * dinteg_dr) / (integral * integral);
      const double dr_R_dr = dr_R_H_dr;

      const double num = (psi_soil[i] - P_x_r - grav_head_z_[i]) * inv_area_leaf;
      const double dnum_dr = -inv_area_leaf;
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
  void uptake_impl(double P_x_r, const std::vector<double>& psi_soil,
                   bool use_integral_cache, double area_leaf,
                   std::vector<double>& soil_consumption, double& E_up) const {

    if (!std::isfinite(P_x_r) || !std::isfinite(area_leaf)) {
      util::stop("E_from_Soil_to_Root_Collar invalid input; P_x_r=" + util::to_string(P_x_r) +
                 "; area_leaf_=" + util::to_string(area_leaf));
    }

    E_up = 0;

    // area_leaf is constant across the whole solve; fold its reciprocal into a
    // per-layer multiply instead of a per-layer division (1 fdiv/call vs 15).
    const double inv_area_leaf = 1.0 / area_leaf;

    // Cumulative-integral spline caching (bit-identical fast path). The only two
    // arguments ever passed to root_vuln_integral_from_psi in the loop below are
    // -P_src_min and -hi_neg, each of which resolves to exactly one of
    // {-psi_soil[i], -P_x_r, 0}. -P_x_r is constant across all layers (compute
    // once), and -psi_soil[i] is constant across the whole solve (precomputed in
    // begin_solve).
    const double neg_P_x_r = -P_x_r;
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
    for(int i = 0; i < max_soil_layer; i++){

    // Find the most negative soil potential out of the given soil layer and the root collar
    double P_src_min = std::min(psi_soil[i], P_x_r);

    // Find the least negative soil potential out of the given soil layer and the root collar
    double P_src_max = std::max(psi_soil[i], P_x_r);

     // If root collar soil water potential equals the soil water potential in a given layer
    if(std::abs(P_x_r - psi_soil[i]) < 1e-8){

      // Fraction of conductance in roots in a given layer at most negative soil water potential (but actually is equal to root collar)
      // root_vuln_from_psi is a pre-built spline of exp(-(|psi|/b_root)^c_root)
      double f_ri = root_vuln_from_psi.eval(-P_src_min);
      if (!std::isfinite(f_ri) || f_ri <= 0.0) {
        util::stop("E_from_Soil_to_Root_Collar invalid f_ri; layer=" + std::to_string(i) +
                   "; f_ri=" + util::to_string(f_ri) +
                   "; P_src_min=" + util::to_string(P_src_min) +
                   "; P_x_r=" + util::to_string(P_x_r));
      }

      // Fraction of conductance in roots in a given layer at most negative soil water potential
      double r_R_H = r_R_H_min[i] / f_ri; // [MPa * s * (mol H2O)^-1]

      // Total root resistance (horizantal plus vertical)
      double r_R = r_R_H + r_R_V_sum[i];

      // Transpiration is equivalent to gravitational water loss (i.e. layer gains water)
      double E_i = -grav_head_z_[i] * inv_area_leaf / r_R ;

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
      double hi_neg = std::min(P_src_max, 0.0); // boundary of the psi<=0 part
      double lo_pos = std::max(P_src_min, 0.0); // boundary of the psi>0 part

      // Memoised cumulative-integral lookup. Returns the exact same double the
      // spline would (same input -> same output); the comparisons select the
      // precomputed value because -P_src_min / -hi_neg are bit-for-bit equal to
      // one of the cached arguments in the common (psi<=0) case.
      const double neg_psi_soil_i = -psi_soil[i];
      auto G_integral = [&](double arg) -> double {
        if (use_integral_cache) {
          if (arg == neg_P_x_r) return G_at_P_x_r;
          if (arg == neg_psi_soil_i) return root_vuln_integral_soil_[i];
        }
        return root_vuln_integral_from_psi.eval(arg);
      };

      double integral = 0.0;
      if (hi_neg > P_src_min) {
        // psi<=0 part: magnitude m runs from -hi_neg up to -P_src_min
        integral += G_integral(-P_src_min) - G_integral(-hi_neg);
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
    double r_R_H = r_R_H_min[i] * span / integral; // [MPa * s * (mol H2O)^-1]

    // Find the total resistance in a given layer by adding the vertical resistance in that layer
    double r_R = r_R_H + r_R_V_sum[i]; // [MPa * s * (mol H2O)^-1]

    // Transpiration is equal to the potentail gradient between the root collar and the soil, accounting for gravitational potential
    double E_i = (psi_soil[i] - P_x_r - grav_head_z_[i]) * inv_area_leaf / r_R; // [mol H2O / m^2 / s]

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
    util::stop("E_from_Soil_to_Root_Collar non-finite E_up_; P_x_r=" + util::to_string(P_x_r) +
               "; max_soil_layer=" + std::to_string(max_soil_layer) +
               "; area_leaf_=" + util::to_string(area_leaf));
  }
  }
};

}  // namespace leaf

#endif
