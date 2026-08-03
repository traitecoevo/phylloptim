// -*-c++-*-
#ifndef LEAF_CLOSED_FORM_HPP_
#define LEAF_CLOSED_FORM_HPP_

// A fast approximate solver for the TF24 leaf optimum, as an alternative to the
// exact golden-section search in Leaf::optimise_psi_stem_TF.
//
// WHY. The exact solve is the model's whole cost. Measured in the companion
// analysis (`Falster-stomatal_analytical_analysis`, notes/tf24_closed_form_bench.cpp,
// benchmarked against the genuine plant::Leaf):
//
//     exact optimise_psi_stem_TF          2.611 us      1x
//     power law + 1 Newton step           0.241 us     10.8x
//     explicit form at beta2 = 1/stem_c   0.056 us     47x
//
// At beta2 = 1/stem_c, 0.051 of the 0.056 us is set_physiology -- the leaf solve itself
// has essentially vanished. For scale: making this library header-only, and
// enabling LTO, each measured at *zero*. This is where the speed is.
//
// HOW. Every model in this family satisfies dA/dE = lambda at the optimum, and
// given lambda the solution collapses to the Medlyn USO form
//
//     ci/ca = xi/(xi + sqrt(D)),     xi = sqrt(Q/lambda)
//
// with Q below. Because lambda_TF24 is a power law in psi to leading order, that
// can be inverted explicitly for a starting psi, and one Newton step on the
// unapproximated supply-minus-demand residual finishes the job. When
// beta2 = 1/stem_c the power cancels and there is nothing left to solve at all.
//
// STATUS: this is NOT wired into Leaf, and nothing calls it by default. It is a
// selectable alternative whose accuracy has to be established per use (see
// within_guard below and PLAN.md item 9) before it can replace the exact solve on
// a production path.
//
// FOUR THINGS NOT TO GET WRONG, all learned the expensive way over there:
//
//  1. `newton_steps = 1` is deliberate. Two steps are *worse* in the tail.
//  2. The validity guard tests an OUTPUT (ci/ca > 0.5), so it can only be applied
//     after the fact, with a fallback to the exact solve. You cannot branch on it
//     up front.
//  3. Quote the realised speedup, not the ceiling. With a fallback fraction phi
//     it is 1/[phi + (1-phi)/10.8] -- 3.6x at phi = 0.2, not 10.8x. Measuring phi
//     on a real water-limited scenario is still an open item.
//  4. Smoothness of the argmax is a hard constraint, not a nicety: plant chose
//     golden-section over Brent precisely because its argmax varies smoothly with
//     inputs, which the demographic growth-rate gradient depends on (see the
//     comment in leaf/optimize.hpp). Any switch to this solver must re-measure
//     that smoothness -- the reference reports roughness in dA/dh of 0.0015 for
//     the closed form against 0.0011 for the exact solve.
//
// LIMITATION: single-layer only, with the collar held at zero
// (`l.transpiration(psi, 0.0)`), matching the reference. The multi-layer case is a
// three-level nest, so the prize there is larger -- but the closed form for it does
// not exist yet. See PLAN.md item 7b.

#include <leaf/constants.hpp>
#include <leaf/leaf_model.hpp>

#include <cmath>

namespace leaf {
namespace closed_form {

// The group that converts a marginal cost of water into the USO slope xi:
//     xi = sqrt(Q/lambda),   Q = 3*Gstar*kg_to_mol_h2o/(1.67e-3)
// with Gstar the CO2 compensation point in Pa and lambda in umol CO2 (kg H2O)^-1.
//
// The 1.67e-3 is the H2O:CO2 diffusion ratio carrying the same 1e-3 scale factor
// that appears in the E-from-A relation below; it is kept in this one place rather
// than spread through the expressions. NOTE that H2O_CO2_stom_diff_ratio is 1.67
// here while the g1 literature uses 1.6 -- see PLAN.md item 8 -- so a g1 compared
// against fitted values carries that 2.2% offset.
inline double uso_group(const Leaf &l) {
  const double gstar_Pa = l.gamma_ * l.umol_per_mol_to_Pa_;
  return 3.0 * gstar_Pa * kg_to_mol_h2o / (H2O_CO2_stom_diff_ratio * 1e-3);
}

// Transpiration implied by assimilation on the DEMAND side, kg H2O m-2 s-1.
// The supply-side counterpart is Leaf::transpiration; the closed form works by
// driving their difference to zero.
inline double transpiration_from_assim(const Leaf &l, double assim, double ci) {
  return H2O_CO2_stom_diff_ratio * 1e-3 * assim * l.atm_vpd_ /
         ((l.ca_ - ci) * kg_to_mol_h2o);
}

// NOTE ON stem_b.value BELOW. This solver uses b as a bare SCALE -- inside a sqrt,
// inside a pow's base, and as a divisor -- not only as the tension ratio psi/b. Those
// uses need b's numeric value and need it positive, which is why b is typed AbsPsi (a
// tension) and NOT as a signed potential: signing it would put a negative under a
// square root at :197 and under a fractional power at :158. See the retraction in
// leaf/potential.hpp.
//
// d(lambda)/d(psi), analytic. Differentiates Leaf::lambda_TF24, whose form is
//     lambda = K * (1-f)^(beta2-1) * p^(stem_c-1),  p = psi/stem_b,  f = exp(-p^stem_c)
inline double dlambda_TF24(const Leaf &l, double psi) {
  const double K =
      l.cost_scale_TF24 * l.beta2 * l.stem_c / (l.stem_b.value * l.leaf_specific_conductance_max_);
  const double p = AbsPsi{psi} / l.stem_b;
  const double f = std::exp(-std::pow(p, l.stem_c));
  return K * std::pow(p, l.stem_c - 2.0) *
         ((l.beta2 - 1.0) * std::pow(1.0 - f, l.beta2 - 2.0) * f * l.stem_c *
              std::pow(p, l.stem_c) +
          std::pow(1.0 - f, l.beta2 - 1.0) * (l.stem_c - 1.0)) /
         l.stem_b.value;
}

// dA/dci for the colimitation quadratic, analytic. Replaces the central
// difference an earlier prototype used, saving two assimilation evaluations per
// Newton step.
inline double dassim_dci(const Leaf &l, double ci, double electron_transport) {
  const double gstar = l.gamma_ * l.umol_per_mol_to_Pa_;
  const double ar = l.vcmax_ * (ci - gstar) / (ci + l.km_);
  const double ae =
      electron_transport / 4.0 * (ci - gstar) / (ci + 2.0 * gstar);
  const double dar =
      l.vcmax_ * (l.km_ + gstar) / ((ci + l.km_) * (ci + l.km_));
  const double dae = electron_transport / 4.0 * 3.0 * gstar /
                     ((ci + 2.0 * gstar) * (ci + 2.0 * gstar));
  const double s = ar + ae, ds = dar + dae;
  const double cv = l.curv_fact_colim;
  const double disc = std::sqrt(s * s - 4.0 * cv * ar * ae);
  const double ddisc =
      (2.0 * s * ds - 4.0 * cv * (dar * ae + ar * dae)) / (2.0 * disc);
  return (ds - ddisc) / (2.0 * cv);
}

struct Solution {
  double psi_stem;      // MPa, positive magnitude (NaN for the explicit form)
  double ci;            // Pa
  double assim;         // umol CO2 m-2 s-1
  double transpiration; // kg H2O m-2 s-1
  double stom_cond_CO2; // mol CO2 m-2 s-1
  double g1_eff;        // the USO slope xi, kPa^0.5
};

// Assemble the outputs implied by a stem potential.
inline Solution evaluate_at(Leaf &l, double psi, double Q, double sqrt_D) {
  const double lambda = l.lambda_TF24(psi);
  const double xi = std::sqrt(Q / lambda);
  const double ci = l.ca_ * xi / (xi + sqrt_D);
  const double assim = l.assim_colimited(ci);
  const double E = transpiration_from_assim(l, assim, ci);
  return Solution{psi, ci, assim, E,
                  l.atm_kpa_ * E * kg_to_mol_h2o / l.atm_vpd_ /
                      H2O_CO2_stom_diff_ratio,
                  xi};
}

// General beta2. Explicit power-law leading order, then `newton_steps` Newton
// steps on the full supply-minus-demand residual. Requires set_physiology to have
// run. See note 1 above: leave newton_steps at 1.
inline Solution solve(Leaf &l, int newton_steps = 1) {
  const double kmax = l.leaf_specific_conductance_max_;
  const double sqrt_D = std::sqrt(l.atm_vpd_);
  const double Q = uso_group(l);
  const double K_lambda = l.cost_scale_TF24 * l.beta2 * l.stem_c / (l.stem_b.value * kmax);
  const double Xi = std::sqrt(Q / K_lambda);
  const double n = l.stem_c * l.beta2 - 1.0;
  const double electron_transport = l.electron_transport();

  // Leading order, taking kappa from the wet-end limit A(ci -> ca).
  const double assim_wet = l.assim_colimited(l.ca_ * (1.0 - 1e-9));
  const double kappa = H2O_CO2_stom_diff_ratio * 1e-3 * assim_wet /
                       (l.ca_ * kg_to_mol_h2o);
  double p = std::pow(kappa * sqrt_D * Xi / (kmax * l.stem_b.value), 2.0 / (n + 2.0));

  for (int k = 0; k < newton_steps; ++k) {
    const double psi = p * l.stem_b.value;
    if (!(psi > 0.0) || psi >= l.psi_crit.value) {
      break;
    }
    const double lambda = l.lambda_TF24(psi);
    const double xi = std::sqrt(Q / lambda);
    const double ci = l.ca_ * xi / (xi + sqrt_D);
    const double assim = l.assim_colimited(ci);
    const double dassim = dassim_dci(l, ci, electron_transport);
    const double u = l.ca_ - ci;
    const double E = transpiration_from_assim(l, assim, ci);
    const double dE_dci = H2O_CO2_stom_diff_ratio * 1e-3 * l.atm_vpd_ /
                          kg_to_mol_h2o * (dassim * u + assim) / (u * u);
    const double dci_dxi = l.ca_ * sqrt_D / ((xi + sqrt_D) * (xi + sqrt_D));
    const double dxi_dpsi = -0.5 * xi / lambda * dlambda_TF24(l, psi);
    // Residual: hydraulic supply minus stomatal demand, both kg H2O m-2 s-1.
    const double R = l.transpiration(psi, 0.0) - E;
    const double dR =
        kmax * std::exp(-std::pow(p, l.stem_c)) - dE_dci * dci_dxi * dxi_dpsi;
    if (dR == 0.0 || !std::isfinite(R) || !std::isfinite(dR)) {
      break;
    }
    double psi_next = psi - R / dR;
    psi_next = std::max(1e-6, std::min(psi_next, l.psi_crit.value * (1.0 - 1e-9)));
    p = AbsPsi{psi_next} / l.stem_b;
  }

  return evaluate_at(l, p * l.stem_b.value, Q, sqrt_D);
}

// beta2 == 1/stem_c. The psi dependence cancels out of lambda entirely, so xi is
// constant and there is nothing to solve -- no power law, no Newton step. This is
// the 47x case. Only correct when beta2 == 1/stem_c; check beta2_is_exact below.
inline Solution solve_exact_beta2(Leaf &l) {
  const double sqrt_D = std::sqrt(l.atm_vpd_);
  const double Q = uso_group(l);
  const double xi = std::sqrt(Q * l.stem_b.value * l.leaf_specific_conductance_max_ /
                              (l.cost_scale_TF24 * l.beta2 * l.stem_c));
  const double ci = l.ca_ * xi / (xi + sqrt_D);
  const double assim = l.assim_colimited(ci);
  const double E = transpiration_from_assim(l, assim, ci);
  return Solution{std::nan(""), ci, assim, E,
                  l.atm_kpa_ * E * kg_to_mol_h2o / l.atm_vpd_ /
                      H2O_CO2_stom_diff_ratio,
                  xi};
}

inline bool beta2_is_exact(const Leaf &l, double tol = 1e-12) {
  return std::abs(l.beta2 - 1.0 / l.stem_c) <= tol * std::max(1.0, 1.0 / l.stem_c);
}

// The validity guard. The closed form degrades where the leaf is far from the
// wet-end limit its leading order is expanded about, and ci/ca is the diagnostic:
// the reference reports good agreement while ci/ca > 0.5 and does not claim it
// below. Tests an OUTPUT, so it can only be applied after solving -- see note 2.
inline bool within_guard(const Leaf &l, const Solution &s) {
  return std::isfinite(s.ci) && std::isfinite(s.assim) && s.ci / l.ca_ > 0.5;
}

} // namespace closed_form
} // namespace leaf

#endif
