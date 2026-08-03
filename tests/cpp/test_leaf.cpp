// Pure-C++ test suite for the leaf model. No R, no test framework, no linking:
//   make -C tests/cpp && ./tests/cpp/test_leaf
//
// The trait values and drivers below are lifted from plant's
// tests/testthat/test-leaf.r so the two suites exercise the same operating
// point. See PLAN.md step 1: the expected values here were produced BY this
// implementation and are regression guards only -- they have not yet been
// cross-checked against plant's compiled build, which is the first job.

#include <leaf.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

int failures = 0;
int checks = 0;

void ok(bool pass, const std::string &what) {
  ++checks;
  if (!pass) {
    ++failures;
    printf("  FAIL  %s\n", what.c_str());
  }
}

void near(double got, double want, double tol, const std::string &what) {
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
  double rho = 608.0;
  double a_bio = 0.0245;
};

leaf::Leaf make_leaf(const Drivers &d, std::vector<double> psi_soil,
                     std::vector<double> soil_depth) {
  leaf::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  std::vector<double> mass_root_prop(psi_soil.size(),
                                     1.0 / double(psi_soil.size()));
  l.set_physiology(d.area_leaf, mass_root_prop, d.rho, d.a_bio, d.PPFD, psi_soil,
                   soil_depth, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                   d.theta * d.h, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// ---------------------------------------------------------------------------

void test_defaults_are_unset() {
  printf("defaults are unset until set_physiology\n");
  leaf::Leaf l;
  ok(!std::isfinite(l.ci_), "ci_ starts unset");
  ok(!std::isfinite(l.assim_colimited_), "assim_colimited_ starts unset");
  ok(!std::isfinite(l.opt_psi_stem_), "opt_psi_stem_ starts unset");
  ok(l.roots_.psi_soil_.empty(), "psi_soil_ starts empty");
  ok(l.use_energy_balance_ == false, "energy balance defaults off");
}

void test_vulnerability_curve() {
  printf("xylem vulnerability curve\n");
  leaf::Leaf l;
  near(l.proportion_of_conductivity(0.0), 1.0, 1e-12,
       "full conductivity at zero potential");
  // psi_crit is NOT the 1%-conductivity point -- with the default b and c it sits
  // at ~5%. (The 1% point is where setup_transpiration ends its knot grid, which
  // is a different and larger potential.) Asserting the closed form rather than a
  // round number so this stays honest if the defaults change.
  near(l.proportion_of_conductivity(l.psi_crit),
       std::exp(-std::pow(l.psi_crit / l.b, l.c)), 1e-12,
       "conductivity at psi_crit matches the Weibull closed form");
  ok(l.proportion_of_conductivity(l.psi_crit) < 0.06,
     "conductivity at psi_crit is a few percent");
  ok(l.proportion_of_conductivity(1.0) > l.proportion_of_conductivity(2.0),
     "conductivity declines monotonically");
}

void test_spline_matches_direct_integration() {
  printf("pre-integrated spline vs direct quadrature\n");
  Drivers d;
  leaf::Leaf l = make_leaf(d, {2.0}, {1.0});
  // This is the check plant makes at test-leaf.r:214. The spline is what the hot
  // path reads; adaptive Simpson integrates the curve directly.
  for (double psi_stem : {2.5, 3.0, 4.0, 5.0}) {
    near(l.transpiration(psi_stem, 2.0),
         l.transpiration_full_integration(psi_stem, 2.0), 1e-6,
         "transpiration at psi_stem=" + std::to_string(psi_stem));
  }
}

void test_arrhenius() {
  printf("temperature response\n");
  leaf::Leaf l;
  near(l.arrh_curve(leaf::vcmax_ha, 100.0, 25.0), 100.0, 1e-12,
       "Arrhenius is the identity at the 25 C reference");
  ok(l.arrh_curve(leaf::vcmax_ha, 100.0, 35.0) > 100.0,
     "Arrhenius rises above the reference temperature");
  ok(l.peak_arrh_curve(leaf::jmax_ha, 100.0, 60.0, leaf::jmax_H_d,
                       leaf::jmax_d_S) <
         l.peak_arrh_curve(leaf::jmax_ha, 100.0, 30.0, leaf::jmax_H_d,
                           leaf::jmax_d_S),
     "peaked Arrhenius declines past its optimum");
}

void test_saturation_vapour_pressure() {
  printf("saturation vapour pressure\n");
  leaf::Leaf l;
  // Tetens at 25 C is ~3.167 kPa.
  near(l.saturation_vapour_pressure(25.0), 3.167, 1e-3, "es(25 C)");
  // Delta ~ 0.189 kPa/K at 25 C.
  near(l.saturation_vapour_pressure_slope(25.0), 0.189, 5e-3, "Delta(25 C)");
  ok(l.saturation_vapour_pressure(30.0) > l.saturation_vapour_pressure(20.0),
     "es increases with temperature");
}

void test_solve_single_layer() {
  printf("root-collar solve, single soil layer\n");
  Drivers d;
  leaf::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  ok(std::isfinite(l.opt_psi_stem_), "psi_stem is finite");
  ok(std::isfinite(l.profit_), "profit is finite");
  ok(l.opt_psi_stem_ > 2.0, "stem is drier than the soil");
  ok(l.opt_psi_stem_ <= l.psi_crit, "stem stays within psi_crit");
  ok(l.transpiration_ > 0.0, "transpiration is positive");
  ok(l.assim_colimited_ > 0.0, "assimilation is positive");
  // Regression guards -- see the note at the top of this file.
  near(l.opt_psi_stem_, 3.595247, 1e-5, "opt_psi_stem_");
  near(l.assim_colimited_, 5.599511, 1e-5, "assim_colimited_");
  near(l.transpiration_, 1.141941e-05, 1e-5, "transpiration_");
  near(l.profit_, 2.515843, 1e-5, "profit_");
}

void test_solve_is_deterministic() {
  printf("repeated solves are bit-identical\n");
  Drivers d;
  leaf::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double first = l.profit_, psi = l.opt_psi_stem_;
  for (int i = 0; i < 50; ++i) {
    l.find_root_collar_psi();
  }
  ok(l.profit_ == first, "profit is unchanged after 50 re-solves");
  ok(l.opt_psi_stem_ == psi, "psi_stem is unchanged after 50 re-solves");
}

void test_drier_soil_costs_carbon() {
  printf("response to drying soil\n");
  Drivers d;
  double prev_profit = 1e9, prev_E = 1e9;
  for (double psi : {0.5, 1.0, 2.0, 3.0, 4.0}) {
    leaf::Leaf l = make_leaf(d, {psi}, {1.0});
    l.find_root_collar_psi();
    ok(l.profit_ <= prev_profit,
       "profit does not rise as soil dries to " + std::to_string(psi));
    ok(l.transpiration_ <= prev_E,
       "transpiration does not rise as soil dries to " + std::to_string(psi));
    prev_profit = l.profit_;
    prev_E = l.transpiration_;
  }
}

void test_light_response() {
  printf("light response\n");
  Drivers dim = Drivers(), bright = Drivers();
  dim.PPFD = 200;
  bright.PPFD = 1800;
  leaf::Leaf a = make_leaf(dim, {2.0}, {1.0});
  leaf::Leaf b = make_leaf(bright, {2.0}, {1.0});
  a.find_root_collar_psi();
  b.find_root_collar_psi();
  ok(b.assim_colimited_ > a.assim_colimited_, "brighter light assimilates more");
  ok(b.transpiration_ > a.transpiration_, "brighter light transpires more");
}

void test_multi_layer_soil() {
  printf("multi-layer soil\n");
  Drivers d;
  leaf::Leaf l = make_leaf(d, {1.0, 2.0, 3.0}, {0.5, 0.5, 0.5});
  l.find_root_collar_psi();
  ok(std::isfinite(l.profit_), "profit is finite with three layers");
  ok(l.soil_consumption_.size() == 3u, "one consumption term per layer");
  double total = 0.0;
  for (double s : l.soil_consumption_) {
    ok(std::isfinite(s), "per-layer consumption is finite");
    total += s;
  }
  ok(total > 0.0, "total soil water uptake is positive");
}

void test_shutdown_when_soil_is_drier_than_psi_crit() {
  printf("shutdown past psi_crit\n");
  Drivers d;
  leaf::Leaf l = make_leaf(d, {12.0}, {1.0});
  l.find_root_collar_psi();
  ok(std::isfinite(l.profit_), "profit stays finite past psi_crit");
  ok(l.profit_ <= 0.0, "profit is non-positive when shut down");
  near(l.opt_psi_stem_, l.psi_crit, 1e-12, "stem is held at psi_crit");
}

// KNOWN DEFECT, inherited from plant unchanged -- see PLAN.md, "Fix the
// shutdown-state leak". set_shutdown_state() writes only root_collar_psi_,
// opt_psi_stem_ and profit_. It does NOT reset transpiration_,
// assim_colimited_, stom_cond_CO2_, ci_, E_up_ or soil_consumption_, and
// set_physiology() does not either (setup_clean_leaf() runs from the
// constructors only). So a Leaf object reused across solves -- which is exactly
// how plant uses it, one persistent Leaf per TF24_Strategy driving every node
// and timestep -- reports the PREVIOUS solve's water and carbon fluxes after a
// shutdown.
//
// This test pins the broken behaviour deliberately, so that the fix is a
// visible, deliberate change rather than a silent one. When the leak is fixed,
// these assertions flip to expecting zeros.
void test_shutdown_leaves_stale_state_known_defect() {
  printf("shutdown state leak (known defect, pinned)\n");
  Drivers d;
  leaf::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  std::vector<double> mrp{1.0}, depth{1.0};
  const auto solve = [&](double psi) {
    std::vector<double> ps{psi};
    l.set_physiology(d.area_leaf, mrp, d.rho, d.a_bio, d.PPFD, ps, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.theta * d.h,
                     d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
  };
  solve(4.0); // wet enough to transpire
  const double E_wet = l.transpiration_;
  const double A_wet = l.assim_colimited_;
  const double S_wet = l.soil_consumption_[0];
  ok(E_wet > 0.0, "the wet solve transpires");

  solve(20.0); // far drier than psi_crit: the leaf shuts down
  ok(l.profit_ < 0.0, "the dry solve is a shutdown (profit < 0)");
  ok(l.transpiration_ == E_wet,
     "DEFECT: transpiration_ still holds the wet value");
  ok(l.assim_colimited_ == A_wet,
     "DEFECT: assim_colimited_ still holds the wet value");
  ok(l.soil_consumption_[0] == S_wet,
     "DEFECT: soil_consumption_ still holds the wet value");

  // A freshly constructed leaf shows the same hole from the other side: the
  // fluxes are never written at all, so they stay at the NA sentinel.
  leaf::Leaf fresh = make_leaf(d, {20.0}, {1.0});
  fresh.find_root_collar_psi();
  ok(!std::isfinite(fresh.transpiration_),
     "DEFECT: a never-transpired shutdown leaves transpiration_ unset");
}

void test_analytic_gradient_matches_finite_difference() {
  printf("analytic dprofit/dpsi_collar vs central difference\n");
  Drivers d;
  leaf::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double p0 = l.root_collar_psi_ < 0 ? -l.root_collar_psi_ : l.root_collar_psi_;
  const double target = std::max(2.2, std::min(p0, l.psi_crit - 0.5));
  const double eps = 1e-5;
  const double analytic = l.dprofit_droot_collar_psi(target);
  const double up = l.evaluate_root_collar_psi(target + eps);
  const double dn = l.evaluate_root_collar_psi(target - eps);
  const double fd = (up - dn) / (2 * eps);
  ok(std::isfinite(analytic), "analytic gradient is finite");
  near(analytic, fd, 2e-3, "analytic gradient matches central difference");
}

// lambda = dA/dE is the first-order condition every model in this family shares,
// so checking the analytic lambda against a finite-difference dA/dE is a check on
// the whole optimisation, not just on one formula.
void test_lambda_equals_dA_dE_single_layer() {
  printf("marginal cost of water: analytic lambda vs dA/dE (stem free)\n");
  Drivers d;
  for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
    leaf::Leaf l = make_leaf(d, {psi_soil}, {1.0});
    // optimise_psi_stem_TF holds the collar fixed at psi_soil_[0] and optimises
    // the stem, so the single-layer lambda is the one that applies here.
    l.optimise_psi_stem_TF();
    const double psi = l.opt_psi_stem_;
    const double eps = 1e-6;
    l.set_leaf_states_rates_from_psi_stem(psi + eps, psi_soil);
    const double A1 = l.assim_colimited_, E1 = l.transpiration_;
    l.set_leaf_states_rates_from_psi_stem(psi - eps, psi_soil);
    const double A0 = l.assim_colimited_, E0 = l.transpiration_;
    const double fd = (A1 - A0) / (E1 - E0);
    l.optimise_psi_stem_TF();
    // Tolerance is set by the optimiser: GSS_tol_abs is 1e-3 on psi, so the
    // first-order condition only holds to about that accuracy.
    near(l.marginal_cost_water() / fd, 1.0, 1e-3,
         "lambda/(dA/dE) at psi_soil=" + std::to_string(psi_soil));
    ok(l.marginal_cost_water() > 0.0, "lambda is positive");
  }
}

// The multi-layer identity from the companion manuscript:
//   lambda_multi = lambda_single * [1 + kmax*f(psi_r)/S],   S = dE_up/dpsi_r
// find_root_collar_psi optimises the COLLAR, so lambda_multi is what it
// equalises. The single-layer lambda should be badly wrong here -- that is the
// point of the correction, and this test pins the size of the error.
void test_multilayer_lambda_identity() {
  printf("multi-layer lambda: series-resistance correction vs dA/dE (collar free)\n");
  Drivers d;
  for (int layers : {1, 3, 5}) {
    std::vector<double> ps(layers), depth(layers), root(layers);
    for (int i = 0; i < layers; ++i) {
      ps[i] = 1.0 + 0.25 * i;
      depth[i] = 1.0 * (i + 1);
      root[i] = 1.0 / layers;
    }
    leaf::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const double single = l.marginal_cost_water();
    const double multi = l.marginal_cost_water_multilayer();

    const double target = -l.root_collar_psi_;
    const double eps = 1e-6;
    l.evaluate_root_collar_psi(target + eps);
    const double A1 = l.assim_colimited_, E1 = l.transpiration_;
    l.evaluate_root_collar_psi(target - eps);
    const double A0 = l.assim_colimited_, E0 = l.transpiration_;
    const double fd = (A1 - A0) / (E1 - E0);

    const std::string tag = std::to_string(layers) + " layers";
    ok(std::isfinite(multi), "lambda_multi is finite, " + tag);
    near(multi / fd, 1.0, 1e-3, "lambda_multi/(dA/dE), " + tag);
    // The correction bracket is >= 1, so the single-layer value must understate.
    ok(multi > single, "lambda_multi exceeds lambda_single, " + tag);
    ok(multi / single > 2.0 && multi / single < 12.0,
       "the correction factor is in the reported 2-12 range, " + tag);
    // And the single-layer value should be badly wrong when the collar is free,
    // which is exactly why the correction matters.
    ok(std::abs(single - fd) / fd > 0.5,
       "lambda_single is >50% wrong when the collar is free, " + tag);
  }
}

void test_g1_eff() {
  printf("equivalent Medlyn slope\n");
  Drivers d;
  leaf::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double g1 = l.g1_eff();
  ok(std::isfinite(g1) && g1 > 0.0, "g1_eff is finite and positive");
  // g1_eff is defined by inverting chi = g1/(g1 + sqrt(D)), so that must hold.
  const double chi = l.ci_ / l.ca_;
  near(g1 / (g1 + std::sqrt(l.atm_vpd_)), chi, 1e-12,
       "g1_eff inverts the USO relation exactly");
  // Drier soil closes stomata, lowering chi and therefore g1_eff.
  leaf::Leaf dry = make_leaf(d, {4.0}, {1.0});
  dry.find_root_collar_psi();
  ok(dry.g1_eff() < g1, "g1_eff falls as the soil dries");
  // And a higher marginal cost of water goes with a lower g1_eff.
  ok(dry.marginal_cost_water() > l.marginal_cost_water(),
     "lambda rises as the soil dries");
}

void test_energy_balance_path_runs() {
  printf("Penman-Monteith energy-balance path\n");
  Drivers d;
  leaf::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double A_prescribed = l.assim_colimited_;

  leaf::Leaf eb = make_leaf(d, {2.0}, {1.0});
  eb.use_energy_balance_ = true;
  eb.wind_speed_ = 2.0;
  eb.d_ = 0.05;
  // ra and Rn are derived in set_physiology, so re-run it with the gate on.
  eb = make_leaf(d, {2.0}, {1.0});
  eb.use_energy_balance_ = true;
  eb.find_root_collar_psi();
  ok(std::isfinite(eb.profit_), "energy-balance profit is finite");
  ok(std::isfinite(eb.assim_colimited_), "energy-balance assimilation is finite");
  ok(eb.assim_colimited_ != A_prescribed,
     "energy balance changes the operating point");
  const double Tleaf = eb.leaf_temp_from_E(eb.transpiration_);
  ok(Tleaf >= leaf::leaf_temp_min && Tleaf <= leaf::leaf_temp_max,
     "leaf temperature stays inside the physical clamp");
}

// The closed-form fast path (leaf/closed_form.hpp). Two things matter: that it is
// actually fast, and that its error is characterised honestly rather than asserted
// to be small.
void test_closed_form() {
  printf("closed-form fast path\n");
  // Reference geometry, matching the companion analysis: collar held at zero,
  // single layer, kmax from height.
  const double eta = 12.0, eta_c = 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta);
  const double theta = 1.0 / 4669.0;
  const auto setp = [&](leaf::Leaf &l, double h, double vpd) {
    std::vector<double> ps{0.0}, dp{1.0}, rt{1.0};
    l.set_physiology(1.0, rt, 608.0, 0.0245, 900.0, ps, dp,
                     1.0 * theta / (h * eta_c), vpd, 40.0, theta * h, 25.0, 21.0,
                     101.3);
  };
  leaf::Leaf l;

  // Near the wet end, where the leading-order expansion is centred, it should be
  // very accurate.
  setp(l, 1.0, 2.0);
  l.optimise_psi_stem_TF();
  const double A_wet = l.assim_colimited_;
  setp(l, 1.0, 2.0);
  const leaf::closed_form::Solution wet = leaf::closed_form::solve(l, 1);
  ok(std::abs(wet.assim / A_wet - 1.0) < 2e-3,
     "closed form is within 0.2% of the exact solve at h=1 m");
  ok(leaf::closed_form::within_guard(l, wet), "h=1 m passes the guard");

  // Error grows steeply as the leaf moves away from the wet end. These bounds
  // record measured behaviour -- they are deliberately loose enough to be stable
  // and tight enough to catch a regression.
  struct Case {
    double h, max_err;
  };
  for (const Case &cs : {Case{3.0, 2e-3}, Case{8.0, 1.5e-2}, Case{12.0, 4e-2}}) {
    setp(l, cs.h, 2.0);
    l.optimise_psi_stem_TF();
    const double A_ex = l.assim_colimited_;
    setp(l, cs.h, 2.0);
    const double A_cf = leaf::closed_form::solve(l, 1).assim;
    ok(std::abs(A_cf / A_ex - 1.0) < cs.max_err,
       "closed-form error is bounded at h=" + std::to_string(cs.h) + " m");
  }

  // The guard is coarse, and saying so is the point. It admits ~2% error in A at
  // h = 12 m (ci/ca ~ 0.55) and only rejects once the error is ~8% (h = 20 m).
  // So it is a filter on gross failure, not an error bound -- and note that the
  // heights it rejects are the dominant canopy trees.
  setp(l, 20.0, 2.0);
  l.optimise_psi_stem_TF();
  const double A_tall = l.assim_colimited_;
  setp(l, 20.0, 2.0);
  const leaf::closed_form::Solution tall = leaf::closed_form::solve(l, 1);
  ok(!leaf::closed_form::within_guard(l, tall), "h=20 m is rejected by the guard");
  ok(std::abs(tall.assim / A_tall - 1.0) > 3e-2,
     "and it is rejected because the error really is large there");

  // The beta2 = 1/c leaf, where xi is constant and nothing needs solving.
  leaf::Leaf exact_leaf(96.0, 2.680147, 3.898245, 5.870283, 2.680147, 3.898245,
                        5.870283, 1.0 / 2.680147, 157.44, 0.30, 0.7, 0.99, 1e-3,
                        100, 1e-3, 1000, 7.5, 3.4e2, 9.4e3);
  ok(leaf::closed_form::beta2_is_exact(exact_leaf),
     "beta2_is_exact recognises beta2 = 1/c");
  ok(!leaf::closed_form::beta2_is_exact(l), "and rejects the default beta2 = 1.5");
  setp(exact_leaf, 5.0, 1.5);
  exact_leaf.optimise_psi_stem_TF();
  const double A_ref = exact_leaf.assim_colimited_;
  setp(exact_leaf, 5.0, 1.5);
  const leaf::closed_form::Solution ex =
      leaf::closed_form::solve_exact_beta2(exact_leaf);
  ok(std::isnan(ex.psi_stem),
     "the explicit form reports no psi_stem -- it never solves for one");
  ok(std::abs(ex.assim / A_ref - 1.0) < 4e-2,
     "the explicit form is within a few percent of the exact solve");

  // Timing, reported rather than asserted: absolute microseconds are
  // machine-dependent, so a hard threshold would be a flaky test.
  const std::vector<double> hs{1, 2, 3, 5, 8, 12}, ds{0.8, 1.0, 1.5, 2.0};
  const int reps = 20000;
  double sink = 0;
  const auto time_it = [&](leaf::Leaf &leaf_ref, auto fn) {
    const auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < reps; ++r) {
      setp(leaf_ref, hs[r % 6], ds[r % 4]);
      sink += fn();
    }
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::micro>(t1 - t0).count() / reps;
  };
  const double t_setp = time_it(l, [&] { return l.ca_; });
  const double t_exact = time_it(l, [&] {
    l.optimise_psi_stem_TF();
    return l.assim_colimited_;
  });
  const double t_cf =
      time_it(l, [&] { return leaf::closed_form::solve(l, 1).assim; });
  const double t_expl = time_it(
      exact_leaf, [&] { return leaf::closed_form::solve_exact_beta2(exact_leaf).assim; });
  printf("    set_physiology %.3f us | exact %.3f us | 1-Newton %.3f us (%.1fx) |"
         " explicit %.3f us (%.1fx)\n",
         t_setp, t_exact, t_cf, t_exact / t_cf, t_expl, t_exact / t_expl);
  ok(t_cf < t_exact, "the closed form is faster than the exact solve");
  ok(sink != 0.0, "timing loop was not optimised away");
}

// The carbon -> resistance map (root_network_from_carbon) is the one piece of
// root *architecture* left in this package; the supply solve itself only ever
// reads r_R_H_min and r_R_V_sum. Testing it directly is the point of having
// pulled it out of MultiLayerRoots -- and it is why the map stayed here rather
// than moving to plant, where the golden file could not reach it.
// The second supply path (issue #2 stage 3). Not wired into Leaf yet -- it exists
// so the concept in stage 2 has two real alternatives to dispatch between, and so
// the dispatch measurement was made against a genuine second type rather than a
// stub the optimiser could see through.
void test_single_potential() {
  printf("single-potential supply path\n");
  leaf::SinglePotential sp;
  sp.set_soil_state(1.5);        // positive magnitude, -MPa
  sp.resistance_ = 2.0e4;
  const double area_leaf = 0.05;

  // begin_solve flips to the signed convention and reports the only potential.
  near(sp.begin_solve(), -1.5, 1e-14, "begin_solve returns the signed potential");
  ok(sp.n_layers() == 1, "single potential writes exactly one layer");

  // Ohm's law, and the sign that matters: a collar drier than the soil draws
  // water UP (positive uptake).
  std::vector<double> consumption(1, 0.0);
  double E_up = 0.0;
  sp.uptake(-2.5, area_leaf, consumption, E_up);
  ok(E_up > 0.0, "a collar drier than the soil draws water up");
  near(E_up, (-1.5 - -2.5) / (sp.resistance_ * area_leaf) * leaf::kg_per_mol_h2o,
       1e-14, "uptake is the Ohm's-law flux");
  ok(consumption[0] > 0.0, "per-layer consumption is filled");

  // A collar WETTER than the soil pushes water back into it. Losing this sign is
  // how hydraulic redistribution silently becomes extra uptake.
  sp.uptake(-0.5, area_leaf, consumption, E_up);
  ok(E_up < 0.0, "a collar wetter than the soil loses water to it");

  // The analytic derivative must match a central difference on uptake, and stay
  // finite everywhere -- unlike MultiLayerRoots there are no branch kinks, so it
  // never asks the caller for a finite-difference fallback.
  const double h = 1e-6, p0 = -2.5;
  double up = 0.0, dn = 0.0;
  sp.uptake(p0 + h, area_leaf, consumption, up);
  sp.uptake(p0 - h, area_leaf, consumption, dn);
  const double fd = (up - dn) / (2.0 * h);
  near(sp.duptake_dpsi(area_leaf), fd, 1e-8, "analytic duptake_dpsi matches FD");
  ok(sp.duptake_dpsi(area_leaf) < 0.0,
     "duptake_dpsi is negative: uptake rises as the collar gets more negative");

  // A zero resistance would be an infinite flux; it is rejected, not returned.
  leaf::SinglePotential bad;
  bad.set_soil_state(1.0);
  bad.begin_solve();
  bool threw = false;
  try {
    bad.uptake(-2.0, area_leaf, consumption, E_up);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "zero resistance throws rather than returning an infinity");
}

void test_root_network_from_carbon() {
  printf("root architecture: carbon -> resistance\n");
  const double beta_H = 3.4e2, beta_V = 9.4e3, dz = 0.5;

  // Closed form, straight from the documented model: carbon splits 1/3 vertical
  // : 2/3 horizontal, r_R_H_min = beta_H/c_r_h, r_R_V = beta_V*dz^2/c_r_v.
  const std::vector<double> carbon{3.0, 6.0, 1.5};
  const auto n = leaf::root_network_from_carbon(carbon, dz, beta_H, beta_V);

  ok(n.r_R_H_min.size() == 3u, "one resistance per rooted layer");
  near(n.r_R_H_min[0], beta_H / (3.0 * 2.0 / 3.0), 1e-14, "r_R_H_min layer 0");
  near(n.r_R_H_min[1], beta_H / (6.0 * 2.0 / 3.0), 1e-14, "r_R_H_min layer 1");
  near(n.r_R_V[2], beta_V * dz * dz / (1.5 / 3.0), 1e-14, "r_R_V layer 2");

  // r_R_V_sum is a running total down the profile, so it is monotone and its
  // last entry is the whole column's vertical resistance.
  ok(n.r_R_V_sum[0] < n.r_R_V_sum[1] && n.r_R_V_sum[1] < n.r_R_V_sum[2],
     "cumulative vertical resistance increases with depth");
  near(n.r_R_V_sum[2], n.r_R_V[0] + n.r_R_V[1] + n.r_R_V[2], 1e-14,
       "r_R_V_sum is the running sum of r_R_V");

  // More carbon is less resistance, in both directions. This is the sign that
  // matters: getting it backwards would make investment in roots harmful.
  const auto rich = leaf::root_network_from_carbon({12.0}, dz, beta_H, beta_V);
  const auto poor = leaf::root_network_from_carbon({3.0}, dz, beta_H, beta_V);
  ok(rich.r_R_H_min[0] < poor.r_R_H_min[0], "more root carbon -> less horizontal resistance");
  ok(rich.r_R_V_sum[0] < poor.r_R_V_sum[0], "more root carbon -> less vertical resistance");

  // Trailing zero-carbon layers are dropped, so the hot loop never visits them.
  const auto trailing = leaf::root_network_from_carbon({3.0, 6.0, 0.0, 0.0}, dz,
                                                      beta_H, beta_V);
  ok(trailing.r_R_H_min.size() == 2u, "trailing rootless layers are dropped");

  // Negative carbon is rejected rather than producing a negative resistance.
  bool threw = false;
  try {
    leaf::root_network_from_carbon({3.0, -1.0}, dz, beta_H, beta_V);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "negative root carbon throws");
}

void test_bad_input_throws() {
  printf("input validation\n");
  Drivers d;
  leaf::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  bool threw = false;
  try {
    std::vector<double> psi_soil{2.0}, depth{1.0, 2.0}, mrp{1.0};
    l.set_physiology(d.area_leaf, mrp, d.rho, d.a_bio, d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.theta * d.h,
                     d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  ok(threw, "mismatched soil vector lengths throw std::runtime_error");
}

void benchmark() {
  printf("\ntiming\n");
  Drivers d;
  leaf::Leaf l = make_leaf(d, {2.0}, {1.0});
  const int N = 20000;
  const auto t0 = std::chrono::steady_clock::now();
  double acc = 0;
  for (int i = 0; i < N; ++i) {
    l.find_root_collar_psi();
    acc += l.profit_;
  }
  const auto t1 = std::chrono::steady_clock::now();
  printf("  %.2f us per find_root_collar_psi() (checksum %.6f)\n",
         std::chrono::duration<double, std::micro>(t1 - t0).count() / N,
         acc / N);
}

} // namespace

int main() {
  test_defaults_are_unset();
  test_vulnerability_curve();
  test_spline_matches_direct_integration();
  test_arrhenius();
  test_saturation_vapour_pressure();
  test_solve_single_layer();
  test_solve_is_deterministic();
  test_drier_soil_costs_carbon();
  test_light_response();
  test_multi_layer_soil();
  test_shutdown_when_soil_is_drier_than_psi_crit();
  test_shutdown_leaves_stale_state_known_defect();
  test_analytic_gradient_matches_finite_difference();
  test_lambda_equals_dA_dE_single_layer();
  test_multilayer_lambda_identity();
  test_g1_eff();
  test_energy_balance_path_runs();
  test_closed_form();
  test_single_potential();
  test_root_network_from_carbon();
  test_bad_input_throws();
  benchmark();

  printf("\n%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
