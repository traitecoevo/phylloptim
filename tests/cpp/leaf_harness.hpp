// The shared harness of test_leaf: the counters, the assertions, the drivers and
// the leaf builders more than one part uses. One program, several sources: each
// part defines its tests and a run_<part>() that main calls.
#pragma once

#include "root_network.hpp"
#include <odelia/tangent.hpp>
#include <cmath>
#include <functional>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <type_traits>
#include <vector>

namespace leaf_test {

// Pure-C++ test suite for the leaf model. No R, no test framework, no linking:
//   make -C tests/cpp && ./tests/cpp/test_leaf
//
// The trait values and drivers below are lifted from plant's
// tests/testthat/test-leaf.r so the two suites exercise the same operating
// point. The expected values here were produced BY this implementation and are
// regression guards rather than independent references -- but the implementation
// has been cross-checked against plant's compiled build and found bit-identical,
// so they are guarding a verified model.

#include <phylloptim.hpp>

inline int failures = 0;

inline int checks = 0;

inline void ok(bool pass, const std::string &what) {
  ++checks;
  if (!pass) {
    ++failures;
    printf("  FAIL  %s\n", what.c_str());
  }
}

// A refusal, checked for its REASON and not only for throwing. The message is
// half the point of every `util::stop` in this package -- several of them exist to
// tell the caller which of two similar-looking configurations they are in -- so a
// bare "it threw" would pass on a refusal that had drifted onto the wrong branch.
// Returns rather than asserting, so the call site names the case.
template <typename F>
bool throws_with(F &&f, const std::string &fragment) {
  try {
    f();
  } catch (const std::exception &e) {
    return std::string(e.what()).find(fragment) != std::string::npos;
  }
  return false;
}

inline void near(double got, double want, double tol, const std::string &what) {
  ++checks;
  const double err = std::abs(got - want);
  const double scale = std::max(1.0, std::abs(want));
  if (!(err / scale <= tol)) {
    ++failures;
    printf("  FAIL  %s: got %.12g, want %.12g (rel err %.3g > %.3g)\n",
           what.c_str(), got, want, err / scale, tol);
  }
}

// Trait values from plant's test-leaf.r.
struct Drivers {
  double theta = 0.000157; // Huber value, m2 sapwood m-2 leaf
  double K_s = 1.0;        // stem-specific conductivity
  double h = 5.0;          // path length, m
  double PPFD = 900.0;
  double atm_vpd = 2.0;
  double ca = 40.0;
  double atm_o2_kpa = 21.0;
  double leaf_temp = 25.0;
  double atm_kpa = 101.3;
  double area_leaf = 0.05;
};

// set_traits' fourteenth argument, R_d_25: dark respiration at 25 C. The class
// default, so a call that is not about respiration can pass it and change nothing.
const double kRd25 = 1.44;

// JS22 s hydraulic unit cost. Inert on every curve these tests exercise -- it is
// read only by the JS22 arm -- but set_traits is positional, so it is named here
// rather than repeated as a literal at a dozen call sites.
const double kGammaJS22 = 1.0;

const double kCMaxA = 0.6;

const double kCMaxB = 0.0;

// The TF24_floor curve's ONLY parameter: its price of water as transpiration goes to
// zero. The hydraulic half is TF24's cost at TF24's own traits, so there is nothing
// else to name. Mirrors kLambdaCF77 in test_golden.cpp and carries the same caveat:
// its scale is set by the leaf, whose own marginal cost of water runs 9e4 to 3e5 at
// these drivers.
const double kFloorLambdaO = 1.5e5;

// The price CF77 and TF24_floor are given where a test only needs the price's own
// term to be visible, not a particular operating point. Both prices are caller
// inputs with no default (NA until set), so any test of a priced curve must name
// one. At 1e4 the upstream term in the collar FOC is ~0.2 umol m^-2 s^-1 MPa^-1,
// three orders above the difference's floor.
const double kPriceForFOC = 1e4;

inline phylloptim::Leaf make_leaf(const Drivers &d, std::vector<double> psi_soil,
                     std::vector<double> soil_depth) {
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  // root carbon per unit leaf area: the old absolute carbon divided by area_leaf
  std::vector<double> mass_root_prop(psi_soil.size(),
                                     1.0 / double(psi_soil.size()) / d.area_leaf);
  l.set_physiology(fixture::root_network(mass_root_prop, soil_depth), d.PPFD, psi_soil, soil_depth,
                   d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                   d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// Build a leaf with the energy-balance gate (and the wind model) configured BEFORE
// set_physiology runs, which is what the PM path needs: `ra_`, `Rn_` and `Tair_` are
// all derived there, and the non-finite-wind check is made there. make_leaf() calls
// set_physiology itself, so setting the gate on its return value is too late for
// anything set_physiology decides.
inline phylloptim::Leaf make_pm_leaf(const Drivers &d, std::vector<double> psi_soil,
                        std::vector<double> soil_depth, bool gate,
                        double wind_speed = 2.0, double leaf_dim = 0.05,
                        bool thermal_cost = false) {
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  l.use_energy_balance_ = gate;
  l.use_thermal_cost_ = thermal_cost;
  l.wind_speed_ = wind_speed;
  l.d_ = leaf_dim;
  // Daylight: the NIR band carries as much energy as the PAR band.
  l.NIR_int_ = d.PPFD / phylloptim::umol_par_per_joule;
  std::vector<double> root(psi_soil.size(),
                           1.0 / double(psi_soil.size()) / d.area_leaf);
  l.set_physiology(fixture::root_network(root, soil_depth), d.PPFD, psi_soil, soil_depth,
                   d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                   d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// The dry end has TWO bounds, closed by different conditions, and which one bound
// is a classification a tally has to be able to report.
//
// ⚠️ AT SHIPPED DEFAULTS ONE ARM NEVER WINS, so it needs a fixture built for it
// rather than a sweep that waits for one. root_psi_crit is 5.87 MPa against a
// continuity root near 1.3, so the root's own limit is never the min until the
// trait that sets it is lowered deliberately. Without this the kind exists and
// nothing reaches it, which is a classification no run can produce -- and while
// the distinction was a private bool, no tally could report it either.
//
// The two arms place the collar by DIFFERENT expressions -- a residual on one and
// a closed form on the other -- so the check is the placement as well as the
// name: at the root-limit arm the collar is root_psi_crit exactly.
// d(profit)/d(pars[par]) through the three calls, at whatever point the leaf
// holds, on the tangent scalar.
inline double profit_row_at_current_point(phylloptim::Leaf &l, int par,
                                   double *residual_row = nullptr) {
  using T = odelia::ode::tangent_scalar<double>;
  phylloptim::leaf_pars<T> pars;
  const phylloptim::leaf_pars<double> seated = l.passive_pars();
  for (std::size_t i = 0; i < pars.size(); ++i) pars[i] = T(seated[i]);
  odelia::ode::set_direction(pars[par], 1.0);
  std::vector<T> soil, r_h, r_v;
  for (int i = 0; i < l.supply_n_layers(); ++i) {
    soil.push_back(T(l.roots_.psi_soil_[std::size_t(i)]));
    r_h.push_back(T(l.roots_.network_.r_R_H_min[std::size_t(i)]));
    r_v.push_back(T(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
  }
  const phylloptim::SupplyInputs<T> supply{soil, r_h, r_v,
                                       pars[phylloptim::par_root_P50],
                                       pars[phylloptim::par_root_c]};
  const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), supply);
  const double curvature =
      l.operating_point_kind() == phylloptim::Leaf::OperatingPointKind::Interior
          ? l.profit_curvature_at_collar<phylloptim::Leaf::CostCurve::TF24>()
          : std::numeric_limits<double>::quiet_NaN();
  phylloptim::Leaf::HeldCoords<T> held;
  const T collar = l.operating_collar<phylloptim::Leaf::CostCurve::TF24, T>(
      draw, pars, curvature, &held);
  const auto out = l.billed_outputs_at<phylloptim::Leaf::CostCurve::TF24, T>(
      collar, draw, pars, &held);
  if (residual_row != nullptr) {
    *residual_row = odelia::ode::derivative_along(out.residual_flux);
  }
  return odelia::ode::derivative_along(out.profit);
}

// Capture the message a throwing call produced, or "" if it did not throw. Used
// by every test below that asserts on what an error SAYS rather than that one
// happened -- see test_out_of_domain_names_the_spline for why the wording is
// load-bearing here.
inline std::string message_of(const std::function<void()>& f) {
  try {
    f();
  } catch (const std::exception& e) {
    return std::string(e.what());
  }
  return std::string();
}

inline bool mentions(const std::string& haystack, const char* needle) {
  return haystack.find(needle) != std::string::npos;
}

// Build a leaf on the single-potential path, which is what the ProfitMax entry
// points require and what Sicangco's model is (one soil potential, no root
// resistance network).
inline phylloptim::Leaf make_single_leaf(const Drivers &d, double psi_soil,
                                  bool gate = false) {
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  l.use_energy_balance_ = gate;
  l.set_supply_single(0.0);
  l.NIR_int_ = d.PPFD / phylloptim::umol_par_per_joule;  // daylight
  phylloptim::RootNetwork rn;
  rn.r_R_V_sum = std::vector<double>{1.0e3};
  rn.r_R_H_min = std::vector<double>{0.0};
  rn.r_R_V = std::vector<double>{1.0e3};
  rn.c_r_V = std::vector<double>{0.0};
  rn.c_r_H = std::vector<double>{0.0};
  l.set_physiology(rn, d.PPFD, {psi_soil}, {1.0}, d.K_s * d.theta / d.h,
                   d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// ===========================================================================
// The benefit link with the energy balance on (#124)
// ---------------------------------------------------------------------------
// With the gate on, dA/dpsi has a second path, through Tleaf(E(psi)), and a
// non-identity link takes the same single h'(A) on it. The two costs that read
// Tleaf, ProfitMax's thermal cost and LeastCost's Vcmax, move with psi the same
// way and take none. Four checks, on the linked curves with the gate on:
//
//   1. dprofit_dpsi_leaf against a swept-step difference of each curve's own
//      objective -- the log of the product or ratio for the three log links,
//      the normalised profit for ProfitMax;
//   2. the same derivative as h'(A) times the identity link's dA/dpsi minus
//      dC/dpsi, with dA/dpsi read off a zero-cost curve: the vignette's
//      link-identity chunk with the gate on, so the h' factor is checked against
//      something other than itself. SOX, JW26 and ProfitMax only: LeastCost's
//      temperature-dependent cost has no independent slope to compare with, so
//      check 1 is its guard;
//   3. the solved leaf optimum is a ROOT of dJ/dpsi, which a bracket-width
//      search does not give;
//   4. the collar route's first-order condition against a difference of the
//      profit it reports, both gate states. A temperature term added without
//      its h' converges there to a point that is not the maximum, with every
//      reported number plausible, so only this comparison sees it.

// Thresholds that put the thermal cost on its steep flank at 25 C air. At the
// class defaults (Sicangco's heatwave pair, T50 50.4) the cost is ~1e-5 wherever
// assimilation is positive under these drivers -- every hotter row has ci
// assigned at ca -- so the term would be compiled and never exercised.
const double kLiveT50 = 32.0, kLiveTcrit = 28.0;

// A single-potential leaf with the energy balance and the wind model configured
// before set_physiology, as make_pm_leaf does for the layered supply. With `tc`
// the thermal cost is on at the live thresholds above.
inline phylloptim::Leaf make_linked_eb_leaf(const Drivers &d, double psi_soil, bool tc) {
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  l.use_energy_balance_ = true;
  l.use_thermal_cost_ = tc;
  if (tc) {
    l.T50_ = kLiveT50;
    l.Tcrit_ = kLiveTcrit;
  }
  l.wind_speed_ = 2.0;
  l.d_ = 0.05;
  l.LeastCost_beta = 146.0;
  l.set_supply_single(0.0);
  l.NIR_int_ = d.PPFD / phylloptim::umol_par_per_joule;
  phylloptim::RootNetwork rn;
  rn.r_R_V_sum = std::vector<double>{1.0e3};
  rn.r_R_H_min = std::vector<double>{0.0};
  rn.r_R_V = std::vector<double>{1.0e3};
  rn.c_r_V = std::vector<double>{0.0};
  rn.c_r_H = std::vector<double>{0.0};
  l.set_physiology(rn, d.PPFD, {psi_soil}, {1.0}, d.K_s * d.theta / d.h,
                   d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// EVERY curve's returned point is the actual maximum of ITS OWN objective, checked
// against a derivative-free scan of that objective. This is the optimality test the
// suite was missing for four of the seven: before it, `TF24`, `CF77` and `ProfitMax`
// had scan oracles and `JS22`, `CMax`, `SOX` and `JW26` had only
// `psi_leaf_optima.tsv` -- which is a REGRESSION baseline and would freeze a wrong
// optimum exactly as happily as a right one.
//
// A scan is the right oracle rather than a first-order check, because a maximum here
// is often at a bound: `dJ/dpsi != 0` at a pinned optimum, so an FOC test would
// report a correct answer as wrong (hazard 11). Comparing OBJECTIVE VALUES asks the
// question that is true in both regimes -- is there any feasible point that beats
// the one returned.
//
// ⚠️ THE SCAN CAN ONLY EVER TIE OR LOSE, which is what makes the direction of the
// assertion meaningful: the solver includes both endpoints and refines within a
// cell, so it should never come back below a grid it contains. A scan that WINS is
// the failure -- it means the search missed a basin.
// ============================================================================
// THE TWO-TERM COST: `Theta(E) = Theta~(psi) + lambda_o*E`
// ----------------------------------------------------------------------------
// The decomposition is canonical -- any cost splits into a part depending on the
// potential alone and a linear price of water -- and the content of this curve is
// that the second part is not zero. Every conductance-loss curve in this package
// has `lambda -> 0` as the drop closes, i.e. prices water as free precisely when
// it is abundant; `lambda_o` is the floor that removes that.
//
// `Theta~` here is TF24's OWN cost, read at TF24's own traits, so this curve is
// TF24 plus CF77's price and has exactly ONE parameter of its own. That is what
// makes the first test below the strongest available: two reductions, asserted
// bit-for-bit, needing no fixture and no tolerance -- and the first of them is a
// ONE-restriction statement about the production model.
// ============================================================================

// Bit-equality that survives NaN, which several reported fields legitimately are.
inline bool same_bits(double a, double b) {
  return (std::isnan(a) && std::isnan(b)) || a == b;
}

// The reported operating point, as one vector, so a reduction can be compared
// field by field rather than on `profit_` alone -- a curve wired to the right
// objective and the wrong output assembly would pass the latter.
inline std::vector<double> reported(const phylloptim::Leaf &l) {
  return {l.opt_psi_leaf_, l.profit_, l.hydraulic_cost_, l.ci_,
          l.assim_colimited_, l.transpiration_, l.stom_cond_CO2_,
          l.lambda_emergent()};
}

inline const char *const reported_names[] = {"opt_psi_leaf", "profit", "hydraulic_cost",
                                "ci", "assim", "transpiration", "gc",
                                "lambda_emergent"};

// --- the ATLS thermal-damage layer (use_thermal_damage_) ----------------------
//
// tests/testthat/test-thermal-damage.R checks what the layer computes. These check
// what only C++ reaches: that the gate-off path reads none of it, that the solve's
// own dprofit carries it, that the two finite-difference helpers put its outputs
// back, and that every derivative path which would not carry it refuses.

// Every input of the layer at a value other than its default.
inline void set_atls_inputs(phylloptim::Leaf &l) {
  l.topt_offset_ = 2.5;
  l.k_i_ = 5.0;
  l.k_i_Ea_ = 2.5e5;
  l.k_i_Tiso_ = 47.0;
  l.k_i_baseline_ = 8.0;
  l.k_rec_ = 0.4;
  l.k_rec_Ea_ = 3e4;
  l.dTopt_max_ = 4.0;
  l.dTimpair_max_ = 2.5;
  l.K_A_ = 2.0;
  l.A_ = 1.5;
  l.I_ = 0.2;
  l.D_ = 0.15;
  l.k_mat_ = 0.05;
  l.k_mat_Ea_ = 2e4;
  l.acclim_alpha_ = 0.03;
  l.acclim_beta_ = 0.2;
  l.acclim_T_ = 28.0;
  l.acclim_s_ = 2.0;
  l.c_acclim_maint_ = 0.05;
  l.c_repair_maint_ = 0.3;
  l.c_protect_maint_ = 0.4;
  l.c_repair_flux_ = 0.8;
}

}  // namespace leaf_test

using namespace leaf_test;

namespace residual {

// The drivers of the manuscript's figure-g0 and of the #183 sweeps.
constexpr double kGmin = 5e-4;   // 0.5 mmol H2O m^-2 s^-1, De Kauwe et al.'s mid value

// A leaf with the residual conductance seated through set_traits BEFORE the
// drivers, as a caller has to seat a trait. `single` puts it on the single
// potential path with the manuscript's series resistance; otherwise the
// layered path of the rest of this file.
inline phylloptim::Leaf make(const Drivers &d, const std::vector<double> &psi_soil,
                      const std::vector<double> &depth, double g_min,
                      bool eb = false, bool single = false,
                      double rho = std::numeric_limits<double>::quiet_NaN(),
                      void (*more)(phylloptim::Leaf &) = nullptr) {
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(phylloptim::Leaf::ncontrol_default);
  l.setup_root_vulnerability(phylloptim::Leaf::ncontrol_default);
  l.set_traits(l.vcmax_25, l.stem_c, l.stem_P50, l.roots_.root_c,
               l.roots_.root_P50, l.TF24_beta2, l.jmax_25, l.phi_J,
               l.curv_fact_elec_trans, l.curv_fact_colim, l.TF24_cost_scale,
               l.R_d_25, l.JS22_gamma, l.CMax_a, l.CMax_b, g_min);
  if (std::isfinite(rho)) l.g_min_co2_ratio_ = rho;
  l.use_energy_balance_ = eb;
  l.NIR_int_ = d.PPFD / phylloptim::umol_par_per_joule;
  if (more != nullptr) more(l);
  if (single) {
    l.set_supply_single(0.0);
    phylloptim::RootNetwork rn;
    rn.r_R_V_sum = std::vector<double>{1.0e3};
    rn.r_R_H_min = std::vector<double>{0.0};
    rn.r_R_V = std::vector<double>{1.0e3};
    rn.c_r_V = std::vector<double>{0.0};
    rn.c_r_H = std::vector<double>{0.0};
    l.set_physiology(rn, d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  } else {
    std::vector<double> root(psi_soil.size(),
                             1.0 / double(psi_soil.size()) / d.area_leaf);
    l.set_physiology(fixture::root_network(root, depth), d.PPFD, psi_soil,
                     depth, d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  }
  return l;
}

}  // namespace residual
