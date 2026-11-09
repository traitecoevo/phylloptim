// The cost curves: each curve's definition, its solve on both routes, its lambda, and its reductions to a parent.

#include "leaf_harness.hpp"

namespace {

// lambda_emergent on the collar route, for EVERY curve. At an interior optimum
// the first-order condition makes the marginal cost of water equal the marginal
// carbon of water along the path the solve moves on, whatever the curve and its
// link: h'(A) dA/dpsi = dC/dpsi gives lambda = dA/dE for the identity, log and
// scaled links alike. So the reported lambda is checked against dA/dE taken by a
// central difference ALONG THE COLLAR PATH: assimilation and transpiration at
// the collar either side of the optimum.
//
// ⚠️ This is what catches a lambda of the wrong curve: a collar solve reporting
// TF24's series-corrected lambda whatever it had maximised read 1.05e6 for CF77
// at a price of 1e5. A check on TF24 alone, or on the leaf route, cannot see
// that.
template <phylloptim::Leaf::CostCurve K>
int check_collar_lambda_emergent(const Drivers &d, double *worst) {
  const struct {
    std::vector<double> psi;
    std::vector<double> depth;
  } soils[] = {{{0.5}, {1.0}}, {{1.0}, {1.0}}, {{0.5, 1.0, 1.5}, {0.5, 1.0, 1.5}}};
  int compared = 0;
  for (const auto &s : soils) {
    phylloptim::Leaf l = make_leaf(d, s.psi, s.depth);
    // 3e5, not kPriceForFOC: at 1e5 and below CF77 pins at psi_crit on all three
    // soils here, and a pinned optimum carries no dA/dE identity to check.
    l.CF77_lambda_ = 3e5;
    l.TF24_floor_lambda_o = kPriceForFOC;
    l.LeastCost_beta = 2.0;
    l.JS22_gamma = kGammaJS22;
    l.CMax_a = kCMaxA;
    l.CMax_b = kCMaxB;
    l.set_model(K, true);
    l.optimise();
    if (l.operating_point_kind() !=
        phylloptim::Leaf::OperatingPointKind::Interior) {
      continue;  // pinned: the gradient is genuinely non-zero there
    }
    const double lambda = l.lambda_emergent();
    const double p0 = l.opt_root_psi_;
    const double h = 1e-4;
    l.template evaluate_root_collar_psi_for<K>(p0 + h);
    const double A_up = l.assim_colimited_, E_up = l.transpiration_;
    l.template evaluate_root_collar_psi_for<K>(p0 - h);
    const double A_dn = l.assim_colimited_, E_dn = l.transpiration_;
    const double dA_dE = (A_up - A_dn) / (E_up - E_dn);
    const double err = std::abs(lambda / dA_dE - 1.0);
    *worst = std::max(*worst, err);
    near(lambda / dA_dE, 1.0, 1e-4,
         phylloptim::Leaf::curve_name(static_cast<int>(K)) +
             " collar lambda_emergent against dA/dE, " +
             std::to_string(s.psi.size()) + " layer(s), psi_soil " +
             std::to_string(s.psi[0]));
    ++compared;
  }
  return compared;
}

void test_collar_lambda_emergent_is_dA_dE() {
  printf("collar lambda_emergent is dA/dE along the collar path, every curve\n");
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  d.PPFD = 1500.0;
  double worst = 0.0;
  int per_curve[9] = {
      check_collar_lambda_emergent<CC::TF24>(d, &worst),
      check_collar_lambda_emergent<CC::CF77>(d, &worst),
      check_collar_lambda_emergent<CC::JS22>(d, &worst),
      check_collar_lambda_emergent<CC::CMax>(d, &worst),
      check_collar_lambda_emergent<CC::SOX>(d, &worst),
      check_collar_lambda_emergent<CC::JW26>(d, &worst),
      check_collar_lambda_emergent<CC::ProfitMax>(d, &worst),
      check_collar_lambda_emergent<CC::TF24_floor>(d, &worst),
      check_collar_lambda_emergent<CC::LeastCost>(d, &worst)};
  int total = 0;
  for (int i = 0; i < 9; ++i) {
    total += per_curve[i];
    // Hazard 16: a curve whose every solve pinned would pass by never comparing.
    ok(per_curve[i] > 0, "curve " + std::to_string(i) + " reached an interior optimum");
  }
  printf("    %d interior optima | worst relative error %.3e\n", total, worst);
}

// Cowan & Farquhar (1977): the objective is `A - lambda*E`, so its first-order
// condition in psi_leaf is dA/dE == lambda EXACTLY -- not approximately, and not
// mediated by a vulnerability curve. That identity IS the model, so this is its
// defining assertion rather than a consistency check on it.
//
// Contrast test_lambda_equals_dA_dE_single_layer above, which asks the same
// question of the TF24 cost and has to go through `lambda_emergent()` to get
// an answer: there lambda is a DERIVED property of the operating point, here it
// is the prescribed input. The two tests look alike and are asking opposite
// questions.
void test_cowan_farquhar_equates_dA_dE_to_lambda() {
  printf("Cowan-Farquhar: dA/dE == the prescribed lambda at the optimum\n");
  Drivers d;
  d.PPFD = 1500.0;
  // Reported so the tolerance below can be set from the measurement rather than
  // guessed, and so a solver change that loosens the condition shows up as a
  // number moving rather than only as a pass/fail.
  double worst = 0.0;
  int interior_rows = 0, bound_rows = 0;
  // ⚠️ lambda's SCALE is set by the model, not chosen freely: it is umol C per kg
  // of water, and at these drivers the leaf's own marginal cost of water runs
  // 9e4 to 3e5 (TF24's `lambda_emergent()`). A lambda far below that band makes
  // water free and every optimum pins wide open at psi_crit; far above it and
  // they all pin shut. Either way the identity below is never exercised, and the
  // test passes while asserting nothing -- so keep these values inside the band
  // and keep the interior-row count in the printed line.
  for (double lambda : {8.0e4, 1.5e5, 3.0e5}) {
    for (double psi_soil : {0.5, 1.0, 2.0}) {
      phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
      l.CF77_lambda_ = lambda;
      l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::CF77>();
      const double psi = l.opt_psi_leaf_;

      // An interior optimum is where the condition applies. At a bound the
      // gradient is genuinely non-zero and dA/dE != lambda is the right answer,
      // so assert the bound instead of loosening the tolerance.
      const bool interior = psi > psi_soil + 1e-6 && psi < l.psi_crit - 1e-6;
      if (!interior) {
        ++bound_rows;
        ok(psi == psi_soil || psi == l.psi_crit,
           "a non-interior optimum sits exactly on a bound, lambda=" +
               std::to_string(lambda));
        continue;
      }
      ++interior_rows;

      const double eps = 1e-6;
      l.set_leaf_states_rates_from_psi_leaf(psi + eps, psi_soil);
      const double A1 = l.assim_colimited_, E1 = l.transpiration_;
      l.set_leaf_states_rates_from_psi_leaf(psi - eps, psi_soil);
      const double A0 = l.assim_colimited_, E0 = l.transpiration_;
      const double fd = (A1 - A0) / (E1 - E0);

      worst = std::max(worst, std::abs(fd / lambda - 1.0));
      // Measured 2.3e-06 on the generating platform. The tolerance is set well
      // above that but far below the 1e-3 an argmax-resolution argument would
      // give, because the condition is STATIONARY: a displaced argmax changes
      // dA/dE only to second order, so this ratio is much better conditioned
      // than the potential it is evaluated at.
      near(fd / lambda, 1.0, 1e-4,
           "dA/dE == lambda at lambda=" + std::to_string(lambda) +
               " psi_soil=" + std::to_string(psi_soil));
    }
  }
  printf("    %d interior rows, %d on a bound | worst |dA/dE / lambda - 1|: %.3e\n",
         interior_rows, bound_rows, worst);
  ok(interior_rows > 0,
     "at least one row has an interior optimum, so the identity was exercised");
}

// THE STRONGER FORM OF THE SAME CLAIM, and the one that makes lambda a common
// currency rather than a per-model parameter: the TF24 cost and the
// Cowan-Farquhar cost are different functions of psi, but both optima satisfy
// dA/dE == lambda. So handing Cowan-Farquhar the lambda that TF24 IMPLIES at its
// own optimum must land it on that same optimum.
//
// This is a much sharper instrument than a finite-difference ratio: it compares
// two argmaxes computed by two objectives, with no derivative approximation in
// between, and it fails if `lambda_emergent()` reports the wrong quantity for
// either curve.
void test_cowan_farquhar_reproduces_the_TF_optimum() {
  printf("Cowan-Farquhar at TF24's implied lambda finds TF24's optimum\n");
  Drivers d;
  d.PPFD = 1500.0;
  double worst = 0.0;
  for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
    phylloptim::Leaf ref = make_leaf(d, {psi_soil}, {1.0});
    ref.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
    const double psi_tf = ref.opt_psi_leaf_;
    const double lambda_implied = ref.lambda_emergent();

    phylloptim::Leaf cf = make_leaf(d, {psi_soil}, {1.0});
    cf.CF77_lambda_ = lambda_implied;
    cf.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::CF77>();

    worst = std::max(worst, std::abs(cf.opt_psi_leaf_ - psi_tf));
    // Measured 1.2e-05 MPa on the generating platform. The tolerance is 5e-4
    // rather than that, because this is an ARGMAX difference and argmax-derived
    // quantities are the sqrt-amplified class: they disagree across platforms by
    // up to 1.4e-04. Tightening to the local measurement would make this fail off
    // macOS/arm64 for a reason that is not a regression.
    near(cf.opt_psi_leaf_, psi_tf, 5e-4,
         "the two optima coincide at psi_soil=" + std::to_string(psi_soil));
  }
  printf("    worst |psi_CF - psi_TF| at TF24's implied lambda: %.3e MPa\n", worst);
}

// The cost is lambda*E and E is exactly zero when the two potentials coincide, so
// the no-flow profit is exactly the no-flow assimilation -- zero in light --
// whatever lambda is. That makes this the one degenerate branch in the file whose
// value can be predicted rather than recorded.
void test_cowan_farquhar_closed_state() {
  printf("Cowan-Farquhar: the closed-stomata profit is zero, for any lambda\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {6.0}, {1.0}); // drier than psi_crit
  l.CF77_lambda_ = 800.0;
  l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::CF77>();
  near(l.opt_psi_leaf_, 6.0, 1e-12, "the stem is held at the soil potential");
  near(l.transpiration_, 0.0, 1e-300, "transpiration is exactly zero");
  near(l.hydraulic_cost_, 0.0, 1e-300, "so the cost is exactly zero");
  // No gas moves, so ci sits at the compensation point and A is zero to the
  // ci root-find's tolerance.
  near(l.profit_, 0.0, 1e-9, "and profit is zero: no gas moves, so A = 0");
  // Every reported field describes that point rather than a search probe -- the
  // property the evaluate-at-the-closed-point convention exists to give.
  near(l.profit_, l.assim_colimited_ - l.hydraulic_cost_, 1e-12,
       "profit == assim - cost holds in the degenerate branch too");
}

void test_cowan_farquhar_refuses_an_unset_lambda() {
  printf("Cowan-Farquhar refuses a cleared lambda\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {1.0}, {1.0});
  ok(std::isfinite(l.CF77_lambda_), "CF77_lambda_ starts at its aligned default");
  l.CF77_lambda_ = std::numeric_limits<double>::quiet_NaN();
  bool threw = false;
  std::string what;
  try {
    l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::CF77>();
  } catch (const std::exception &e) {
    threw = true;
    what = e.what();
  }
  ok(threw, "it refuses rather than maximising a NaN objective");
  ok(what.find("CF77_lambda_") != std::string::npos, "and the message names CF77_lambda_");

  // So do the collar route's solve, its prescribed collar and its replay, which
  // would otherwise report a NaN objective's value at the point asked for.
  phylloptim::Leaf c = make_leaf(d, {1.0}, {1.0});
  c.set_model(phylloptim::Leaf::CostCurve::CF77, true);
  c.optimise();
  const double collar = c.opt_root_psi_;
  const auto kind = c.operating_point_kind();
  c.CF77_lambda_ = std::numeric_limits<double>::quiet_NaN();
  const auto refuses = [&](auto call) {
    try {
      call();
    } catch (const std::exception &e) {
      return std::string(e.what()).find("CF77_lambda_") != std::string::npos;
    }
    return false;
  };
  ok(refuses([&] { c.optimise(); }), "the collar solve refuses it, naming it");
  ok(refuses([&] { c.evaluate_root_collar_psi(collar); }),
     "a prescribed collar refuses it, naming it");
  ok(refuses([&] { c.replay_operating_point(collar, kind); }),
     "a replay refuses it, naming it");
}

// ProfitMax's solve on a reused leaf writes hydraulic_cost_, as NA, on every
// exit: its objective is unitless, so the field holds no cost of its own, and
// left alone it reports the last curve's (hazard 8). The dark leaf takes the
// degenerate exit, which writes its outputs by hand.
void test_profitmax_writes_hydraulic_cost_on_every_exit() {
  printf("ProfitMax writes hydraulic_cost_ on a reused leaf, on every exit\n");
  for (double ppfd : {0.0, 900.0}) {
    for (bool collar : {false, true}) {
      Drivers d;
      d.PPFD = ppfd;
      phylloptim::Leaf l = make_leaf(d, {0.5}, {1.0});
      l.set_model(phylloptim::Leaf::CostCurve::TF24, collar);
      l.optimise();
      ok(std::isfinite(l.hydraulic_cost_), "TF24 leaves a finite cost behind");
      l.set_model(phylloptim::Leaf::CostCurve::ProfitMax, collar);
      l.optimise();
      const std::string at = std::string(collar ? " (collar" : " (leaf") +
                             " route, PPFD " + std::to_string(int(ppfd)) + ")";
      if (ppfd == 0.0 && !collar) {
        ok(l.operating_point_kind() == phylloptim::Leaf::OperatingPointKind::Unsolved,
           "the dark leaf takes the degenerate exit" + at);
      }
      ok(std::isnan(l.hydraulic_cost_), "ProfitMax reports hydraulic_cost_ as NA" + at);
    }
  }
}

// EVERY cost curve reports its emergent lambda on the same axis, which is what
// makes `lambda_emergent_` a shared output rather than three fields wearing one
// name. Two independent checks, one per curve:
//
//   * Cowan-Farquhar prices water at a constant, so its emergent lambda IS the
//     prescribed one -- exactly, since dC/dpsi = lambda*(dE/dpsi) identically;
//   * TF24's must equal a finite difference of its own cost over its own E.
//
// The ProfitMax case is checked the same way in
// test_profitmax_emergent_lambda_matches_a_finite_difference.
// ProfitMax's emergent lambda has a derivation rather than a closed form lying
// around, so it gets its own check: against a finite difference of its OWN cost
// over its own transpiration, scaled by |A|max to restore carbon units. Run with
// both gates on, which is the arm where the thermal term is not zero.
void test_profitmax_emergent_lambda_matches_a_finite_difference() {
  printf("ProfitMax's emergent lambda against a finite difference of dC/dE\n");
  Drivers d;
  d.PPFD = 1500.0;
  const double h = 1e-3;
  double worst = 0.0;
  int rows = 0;

  for (double t : {25.0, 40.0}) {
    for (bool eb : {false, true}) {
      d.leaf_temp = t;
      phylloptim::Leaf l = make_leaf(d, {1.0}, {1.0});
      l.use_energy_balance_ = eb;
      l.use_thermal_cost_ = true;
      if (eb) { l.Rn_ = 300.0; l.d_ = 0.05; l.wind_speed_ = 2.0; }
      l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::ProfitMax>();
      const double psi = l.opt_psi_leaf_;
      if (psi <= 1.0 + h || psi >= l.psi_crit - h) continue;

      // dC/dpsi of the NORMALISED cost, and dE/dpsi, both by difference.
      auto cost_and_E = [&](double p) {
        l.profit_psi_leaf_ProfitMax(p, 1.0);
        return std::pair<double, double>(l.hydraulic_cost_norm_ + l.thermal_cost_,
                                        l.transpiration_);
      };
      const auto up = cost_and_E(psi + h);
      const auto dn = cost_and_E(psi - h);
      const double fd = l.profitmax_A_max() *
                        ((up.first - dn.first) / (up.second - dn.second));

      l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::ProfitMax>();
      ++rows;
      worst = std::max(worst, std::abs(l.lambda_emergent() / fd - 1.0));
      near(l.lambda_emergent() / fd, 1.0, 5e-3,
           "ProfitMax emergent lambda at T=" + std::to_string(t) +
               " eb=" + std::to_string(int(eb)));
    }
  }
  printf("    %d rows | worst relative error: %.3e\n", rows, worst);
  ok(rows > 0, "at least one row had an interior optimum to check");
}

void test_every_curve_reports_an_emergent_lambda() {
  printf("every cost curve reports its emergent lambda on one axis\n");
  Drivers d;
  d.PPFD = 1500.0;
  const double h = 1e-3;   // measured step; see the dprofit test for why

  for (double psi_soil : {0.5, 1.0, 2.0}) {
    // --- Cowan-Farquhar: its price, to the stem table's slope accuracy --------
    {
      phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
      l.CF77_lambda_ = 1.5e5;
      l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::CF77>();
      near(l.lambda_emergent(), l.CF77_lambda_, 1e-7,
           "Cowan-Farquhar's emergent lambda is its prescribed one");
    }
    // --- TF24: against a finite difference of dC/dE -------------------------
    {
      phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
      l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
      const double psi = l.opt_psi_leaf_;
      if (psi <= psi_soil + h || psi >= l.psi_crit - h) continue;
      const double dC =
          (l.hydraulic_cost_TF24(psi + h) - l.hydraulic_cost_TF24(psi - h)) / (2 * h);
      l.profit_psi_leaf_TF24(psi + h, psi_soil);
      const double E1 = l.transpiration_;
      l.profit_psi_leaf_TF24(psi - h, psi_soil);
      const double E0 = l.transpiration_;
      const double fd = dC / ((E1 - E0) / (2 * h));
      l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
      near(l.lambda_emergent() / fd, 1.0, 1e-4,
           "TF24's emergent lambda is dC/dE at psi_soil=" +
               std::to_string(psi_soil));
    }
  }
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
  // The root network is built ONCE, outside the closure. It does not depend on h
  // or vpd, and the timing loop below reports a "set_physiology" figure -- if the
  // architecture model were rebuilt per call that figure would be measuring the
  // caller's helper as much as the leaf's setter (measured: 0.106 -> 0.243 us,
  // i.e. more than half the reported cost would have been the helper).
  const std::vector<double> depth_1m{1.0};
  const phylloptim::RootNetwork net_1m = fixture::root_network({1.0}, depth_1m);
  const auto setp = [&](phylloptim::Leaf &l, double h, double vpd) {
    std::vector<double> ps{0.0}, dp{1.0};
    l.set_physiology(net_1m, 900.0, ps, dp, 1.0 * theta / (h * eta_c), vpd, 40.0,
                     25.0, 21.0, 101.3);
  };
  phylloptim::Leaf l = fixture::multilayer_leaf();

  // Near the wet end, where the leading-order expansion is centred, it should be
  // very accurate.
  setp(l, 1.0, 2.0);
  l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
  const double A_wet = l.assim_colimited_;
  setp(l, 1.0, 2.0);
  // ⚠️ `solve` TAKES psi_upstream NOW, and the 0.0 is the reference geometry's
  // collar-at-zero rather than a placeholder. It used to take only the Newton
  // count, so `solve(l, 1)` compiled and meant something else entirely -- one
  // Newton step became a 1 MPa upstream potential. That is the whole reason the
  // upstream potential is an explicit argument instead of being read off the leaf.
  const phylloptim::closed_form::Solution wet =
      phylloptim::closed_form::solve(l, 0.0, 1);
  ok(std::abs(wet.assim / A_wet - 1.0) < 2e-3,
     "closed form is within 0.2% of the exact solve at h=1 m");
  ok(phylloptim::closed_form::within_guard(l, wet), "h=1 m passes the guard");

  // Error grows steeply as the leaf moves away from the wet end. These bounds
  // record measured behaviour -- they are deliberately loose enough to be stable
  // and tight enough to catch a regression.
  struct Case {
    double h, max_err;
  };
  for (const Case &cs : {Case{3.0, 2e-3}, Case{8.0, 1.5e-2}, Case{12.0, 4e-2}}) {
    setp(l, cs.h, 2.0);
    l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
    const double A_ex = l.assim_colimited_;
    setp(l, cs.h, 2.0);
    const double A_cf = phylloptim::closed_form::solve(l, 0.0, 1).assim;
    ok(std::abs(A_cf / A_ex - 1.0) < cs.max_err,
       "closed-form error is bounded at h=" + std::to_string(cs.h) + " m");
  }

  // The guard is coarse, and saying so is the point. It admits ~2% error in A at
  // h = 12 m (ci/ca ~ 0.55) and only rejects once the error is ~8% (h = 20 m).
  // So it is a filter on gross failure, not an error bound -- and note that the
  // heights it rejects are the dominant canopy trees.
  setp(l, 20.0, 2.0);
  l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
  const double A_tall = l.assim_colimited_;
  setp(l, 20.0, 2.0);
  const phylloptim::closed_form::Solution tall =
      phylloptim::closed_form::solve(l, 0.0, 1);
  ok(!phylloptim::closed_form::within_guard(l, tall), "h=20 m is rejected by the guard");
  ok(std::abs(tall.assim / A_tall - 1.0) > 3e-2,
     "and it is rejected because the error really is large there");

  // The TF24_beta2 = 1/c leaf, where xi is constant and nothing needs solving.
  phylloptim::Leaf exact_leaf(96.0, 2.680147, 3.4, 2.680147, 3.4, 1.0 / 2.680147, 157.44, 0.30 / 0.85, 0.7, 0.99, 1e-3,
                        100, 1e-3, 1000, 7.5);
  exact_leaf.set_supply_multilayer();
  exact_leaf.set_model(phylloptim::Leaf::CostCurve::TF24, true);
  ok(phylloptim::closed_form::beta2_is_exact(exact_leaf),
     "beta2_is_exact recognises TF24_beta2 = 1/stem_c");
  ok(!phylloptim::closed_form::beta2_is_exact(l), "and rejects the default TF24_beta2 = 1.5");
  setp(exact_leaf, 5.0, 1.5);
  exact_leaf.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
  const double A_ref = exact_leaf.assim_colimited_;
  setp(exact_leaf, 5.0, 1.5);
  const phylloptim::closed_form::Solution ex =
      phylloptim::closed_form::solve_exact_beta2(exact_leaf, 0.0);
  // The explicit form recovers its potential from the supply; the check is the
  // round trip, the one statement about it that does not restate the code.
  ok(std::isfinite(ex.psi_leaf), "the explicit form recovers a psi_leaf");
  ok(std::abs(exact_leaf.transpiration(ex.psi_leaf, 0.0) / ex.transpiration -
              1.0) < 1e-6,
     "and it is the potential whose hydraulic supply carries that E");
  ok(std::abs(ex.assim / A_ref - 1.0) < 4e-2,
     "the explicit form is within a few percent of the exact solve");
}

// The aligned defaults (#156): at the reference environment `align_curves()`
// uses, every curve with a level parameter has TF24's conductance at its class
// default. This reads the C++ initialisers directly, where the R side's test
// reads `leaf_traits()`, so the two copies of the four numbers are each checked
// against the property rather than against each other.
void test_aligned_defaults_match_tf24_at_the_reference() {
  printf("the aligned level defaults give TF24's conductance at the reference\n");
  using K = phylloptim::Leaf::CostCurve;
  Drivers d;  // PPFD 900, VPD 2 kPa, ca 40 Pa, 25 C, kmax 3.14e-5: the reference
  auto gc_of = [&](auto curve) {
    phylloptim::Leaf l = fixture::multilayer_leaf();
    l.setup_transpiration(phylloptim::Leaf::ncontrol_default);
    l.setup_root_vulnerability(phylloptim::Leaf::ncontrol_default);
    l.set_supply_single(0.0);
    phylloptim::RootNetwork rn;
    rn.r_R_V_sum = std::vector<double>{1.0e3};
    rn.r_R_H_min = std::vector<double>{0.0};
    rn.r_R_V = std::vector<double>{1.0e3};
    rn.c_r_V = std::vector<double>{0.0};
    rn.c_r_H = std::vector<double>{0.0};
    l.set_physiology(rn, d.PPFD, {1.0}, {1.0}, d.K_s * d.theta / d.h,
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    l.optimise_psi_leaf_single<decltype(curve)::value>();
    return l.stom_cond_CO2_;
  };
  const double tf24 = gc_of(std::integral_constant<K, K::TF24>{});
  ok(tf24 > 0.0, "TF24 transpires at the reference");
  // 1e-6: the defaults carry ten significant figures, and the alignment is
  // re-solved in R against a 400-knot spline, as here.
  near(gc_of(std::integral_constant<K, K::CF77>{}) / tf24, 1.0, 1e-6,
       "CF77 at its default price");
  near(gc_of(std::integral_constant<K, K::JS22>{}) / tf24, 1.0, 1e-6,
       "JS22 at its default gamma");
  near(gc_of(std::integral_constant<K, K::CMax>{}) / tf24, 1.0, 1e-6,
       "CMax at its default slope");
  near(gc_of(std::integral_constant<K, K::LeastCost>{}) / tf24, 1.0, 1e-6,
       "LeastCost at its default beta");
}

// ===========================================================================
// Sperry (2017) ProfitMax
// ---------------------------------------------------------------------------
// THE LOAD-BEARING CLAIM is that the normalised objective and this package's
// older `A - lambda*cost` form are the same function up to a positive scale, so
// they share an argmax when lambda = |A|max/(k_soil - kcrit). Everything the
// Sicangco et al. (2026) replication does rests on it.
// The normalised objective IS a constant-lambda objective, up to a positive
// scale. Asserted as the pointwise identity
//
//   A_max * ProfitMax(psi) + A_max * TC  ==  A(psi) - lambda* (k(psi_s) - k(psi))
//
// rather than by comparing two argmaxes. That is both tighter -- the identity
// holds to rounding, where two searches agree only to their bracket width -- and
// stronger, since a positive affine map preserves the argmax, so the identity
// implies the argmax claim while the converse does not.
//
// The right-hand side is built from `assim_colimited_` and
// `proportion_of_conductivity`, not from the normalisers ProfitMax cached, so the
// two sides do not share their arithmetic.
void test_profitmax_is_a_constant_lambda_objective() {
  printf("ProfitMax is the constant-lambda objective, up to a positive scale\n");
  Drivers d;
  d.PPFD = 1500.0;
  double worst = 0.0;
  int points = 0;

  for (double psi_soil : {0.5, 2.0}) {
    for (double t : {25.0, 40.0}) {
      d.leaf_temp = t;
      phylloptim::Leaf l = make_single_leaf(d, psi_soil);
      l.prepare_profitmax();
      const double A_max = l.profitmax_A_max_, k_span = l.profitmax_k_span_;
      if (!(A_max > 0.0) || !(k_span > 0.0)) continue;
      const double lambda_star = A_max / k_span;
      const double k_soil =
          l.leaf_specific_conductance_max_ *
          l.proportion_of_conductivity(psi_soil);

      for (double frac : {0.2, 0.4, 0.6, 0.8}) {
        const double psi = psi_soil + frac * (l.psi_crit - psi_soil);
        const double lhs_norm = l.profit_psi_leaf_ProfitMax(psi, psi_soil);
        // Read AFTER the call: these describe the point just evaluated.
        const double A = l.assim_colimited_, TC = l.thermal_cost_;
        const double k =
            l.leaf_specific_conductance_max_ * l.proportion_of_conductivity(psi);

        const double lhs = A_max * lhs_norm + A_max * TC;
        const double rhs = A - lambda_star * (k_soil - k);
        const double scale = std::max(std::abs(lhs), 1.0);
        worst = std::max(worst, std::abs(lhs - rhs) / scale);
        ++points;
        near(lhs / scale, rhs / scale, 1e-12,
             "the identity holds at frac=" + std::to_string(frac) +
                 " psi_soil=" + std::to_string(psi_soil));
      }
      ok(std::isfinite(lambda_star) && lambda_star > 0.0,
         "and the reported lambda is finite and positive");
    }
  }
  printf("    %d points | worst relative departure from the identity: %.3e\n",
         points, worst);
}

void test_profitmax_normalisation() {
  printf("ProfitMax: what the normalisation does and does not remove\n");
  Drivers d;
  d.PPFD = 1500.0;

  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  const std::vector<double> curve = l.profitmax_curve(101);
  const std::size_t n = 101;
  ok(curve.size() == 5 * n, "profitmax_curve returns five columns");

  // HC runs from 0 at the soil potential to 1 at psi_crit, monotonically. That is
  // the definition and it is what makes the cost unable to vanish.
  near(curve[2 * n + 0], 0.0, 1e-12, "HC is zero at the soil potential");
  near(curve[2 * n + (n - 1)], 1.0, 1e-9, "HC is one at psi_crit");
  bool hc_monotone = true;
  for (std::size_t i = 1; i < n; ++i) {
    if (!(curve[2 * n + i] >= curve[2 * n + i - 1])) hc_monotone = false;
  }
  ok(hc_monotone, "HC increases monotonically along the supply stream");

  // ⚠️ AND HC DOES NOT DEPEND ON kmax AT ALL. Both the numerator and the
  // denominator carry one factor of leaf_specific_conductance_max_, so it cancels
  // exactly. That is worth pinning because Sicangco's ProfitMaxkmax(T) arm gives
  // kmax a temperature response and the paper describes it as changing the cost:
  // it does not. It changes the SUPPLY, and reaches the cost only through which
  // potentials the leaf can reach and how hot it gets there.
  Drivers d2 = d;
  d2.K_s = d.K_s * 3.0;
  phylloptim::Leaf l2 = make_single_leaf(d2, 0.5);
  const std::vector<double> curve2 = l2.profitmax_curve(101);
  double hc_worst = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    hc_worst = std::max(hc_worst, std::abs(curve2[2 * n + i] - curve[2 * n + i]));
  }
  ok(hc_worst < 1e-12, "HC is invariant to a 3x change in kmax");
  printf("    HC under kmax x3: worst difference %.3e\n", hc_worst);

  // The optimiser lands on the curve's own maximum.
  l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  double best = -1e300;
  std::size_t at = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (curve[4 * n + i] > best) { best = curve[4 * n + i]; at = i; }
  }
  ok(std::abs(l.opt_psi_leaf_ - curve[at]) < 2.0 * (curve[1] - curve[0]),
     "the optimiser lands within a grid step of the curve's maximum");
  ok(l.profit_ >= best - 1e-9, "and at no lower profit than the grid's best");
}

void test_profitmax_thermal_cost() {
  printf("ProfitMax: the thermal cost, and that it is inert when off\n");
  Drivers d;
  d.PPFD = 1500.0;
  d.leaf_temp = 48.0;

  phylloptim::Leaf off = make_single_leaf(d, 0.5);
  off.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  ok(off.thermal_cost_ == 0.0, "gate off: TC is exactly zero");

  phylloptim::Leaf on = make_single_leaf(d, 0.5);
  on.use_thermal_cost_ = true;
  on.T50_ = 50.4;
  on.Tcrit_ = 46.5;
  // The gate reaches jmax_ through the temperature block, so the drivers have to
  // be re-supplied for the cache key to notice. Doing it the way a caller would.
  on = make_single_leaf(d, 0.5);
  on.use_thermal_cost_ = true;
  on.T50_ = 50.4;
  on.Tcrit_ = 46.5;
  phylloptim::RootNetwork rn;
  rn.r_R_V_sum = std::vector<double>{1.0e3};
  rn.r_R_H_min = std::vector<double>{0.0};
  rn.r_R_V = std::vector<double>{1.0e3};
  rn.c_r_V = std::vector<double>{0.0};
  rn.c_r_H = std::vector<double>{0.0};
  on.set_physiology(rn, d.PPFD, {0.5}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd,
                    d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  on.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  ok(on.thermal_cost_ > 0.0 && on.thermal_cost_ < 1.0,
     "gate on at 48 C: TC is in (0,1)");
  ok(on.jmax_ < off.jmax_, "and Jmax is scaled down by (1 - TC)");
  near(on.jmax_, off.jmax_ * (1.0 - on.thermal_cost_at(d.leaf_temp)), 1e-12,
       "by exactly that factor");
  printf("    Tleaf %.1f C, Tcrit %.1f, T50 %.1f: TC %.4f, Jmax %.3f -> %.3f\n",
         d.leaf_temp, on.Tcrit_, on.T50_, on.thermal_cost_, off.jmax_, on.jmax_);

  // ⚠️ THE COST AT Tcrit IS A FIXED 11.9%, AND THAT CONTRADICTS THE PAPER'S PROSE.
  // Sicangco et al. write that "(1 - TC) equals one for temperatures below Tcrit,
  // [so] Eqns 10 and 11 yield the same result under such conditions". The equation
  // they cite does not do that: with r = 2/(T50 - Tcrit), the argument at Tcrit is
  // exactly -2 whatever the two thresholds are, so
  //
  //     TC(Tcrit) = 1/(1 + e^2) = 0.1192...
  //
  // independent of parameterisation. A leaf sitting AT its critical temperature
  // has already lost 11.9% of Jmax in this model, and a leaf 4 K below it still
  // pays 3%. Pinned as an identity because it is a property of the functional
  // form rather than of the values in Table 2.
  near(on.thermal_cost_at(on.Tcrit_), 1.0 / (1.0 + std::exp(2.0)), 1e-14,
       "TC at Tcrit is 1/(1+e^2), not zero");
  {
    phylloptim::Leaf wide = make_single_leaf(d, 0.5);
    wide.use_thermal_cost_ = true;
    wide.T50_ = 55.0;
    wide.Tcrit_ = 43.0;
    near(wide.thermal_cost_at(wide.Tcrit_), on.thermal_cost_at(on.Tcrit_), 1e-14,
         "and is the same 11.9% for a threshold pair three times as wide");
  }
  near(on.thermal_cost_at(on.T50_), 0.5, 1e-12, "and is exactly 0.5 at T50");
}

// The collar route for the linked curves, both gate states. Against a difference
// of the profit the solve REPORTS: the product for the log links, so the
// derivative -- which is of the log -- is compared with a difference of log P.
template <phylloptim::Leaf::CostCurve K>
int check_collar_foc_matches_profit_linked(const Drivers &d, bool eb, bool tc,
                                           double &worst) {
  const struct {
    std::vector<double> psi;
    std::vector<double> depth;
  } soils[] = {{{1.0}, {1.0}}, {{1.0, 2.0, 3.0}, {0.5, 1.0, 1.5}}};
  constexpr bool log_link = phylloptim::Leaf::benefit_link<K>() ==
                            phylloptim::Leaf::BenefitLink::Log;
  int compared = 0;
  for (const auto &s : soils) {
    phylloptim::Leaf l = make_pm_leaf(d, s.psi, s.depth, eb, 2.0, 0.05, tc);
    l.LeastCost_beta = 146.0;
    if (tc) {
      l.T50_ = kLiveT50;
      l.Tcrit_ = kLiveTcrit;
    }
    l.set_model(K, true);
    l.optimise();
    if constexpr (K == phylloptim::Leaf::CostCurve::ProfitMax) {
      if (!(l.profitmax_A_max() > 0.0)) {
        continue;
      }
    }
    const double p0 = l.opt_root_psi_;
    if (!std::isfinite(p0)) {
      continue;
    }
    for (double target : {p0 - 0.05, p0 - 0.1}) {
      const double eps = 1e-5;
      bool f = false;
      const double analytic =
          l.template dprofit_droot_collar_psi_for<K>(target, &f);
      if (!f || !std::isfinite(analytic)) {
        continue;
      }
      double up = l.template evaluate_root_collar_psi_for<K>(target + eps);
      double dn = l.template evaluate_root_collar_psi_for<K>(target - eps);
      if (!std::isfinite(up) || !std::isfinite(dn)) {
        continue;
      }
      if constexpr (log_link) {
        if (!(up > 0.0 && dn > 0.0)) {
          continue;
        }
        up = std::log(up);
        dn = std::log(dn);
      }
      const double fd = (up - dn) / (2 * eps);
      const double err = std::abs(analytic - fd) / std::max(1.0, std::abs(fd));
      worst = std::max(worst, err);
      near(analytic, fd, 1e-4,
           phylloptim::Leaf::curve_name(static_cast<int>(K)) +
               (tc ? "+TC" : "") + " eb=" + std::to_string(int(eb)) +
               " at collar " + std::to_string(target) + ", " +
               std::to_string(s.psi.size()) + " layer(s)");
      ++compared;
    }
  }
  return compared;
}

void test_collar_foc_matches_profit_linked_curves() {
  printf("collar dprofit is the derivative of the reported profit, linked curves, "
         "both gate states\n");
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  d.PPFD = 1500.0;
  d.leaf_temp = 25.0;
  int compared = 0;
  double worst_off = 0.0, worst_on = 0.0;
  for (bool eb : {false, true}) {
    double &worst = eb ? worst_on : worst_off;
    for (bool tc : {false, true}) {
      compared += check_collar_foc_matches_profit_linked<CC::SOX>(d, eb, tc, worst);
      compared += check_collar_foc_matches_profit_linked<CC::JW26>(d, eb, tc, worst);
      compared +=
          check_collar_foc_matches_profit_linked<CC::ProfitMax>(d, eb, tc, worst);
      compared +=
          check_collar_foc_matches_profit_linked<CC::LeastCost>(d, eb, tc, worst);
    }
  }
  printf("    %d points compared | worst |error| gate off %.3e, gate on %.3e\n",
         compared, worst_off, worst_on);
  ok(compared >= 36, "enough points were compared on each route (hazard 16)");
}

// The closed form is confined to the prescribed-temperature path (#116).
//
// `Leaf::optimise_closed` refuses the energy balance, and
// `test_closed_form_is_a_selectable_method` covers that route. This test covers
// the HEADER's own entry points, which a bench, a test or a consumer can call
// directly: each of the four solvers must refuse too, from the same stop, because
// the alternative is a wrong answer that `within_guard` passes. Measured with the
// gate forced through, against a 4000-point argmax of the exact objective,
// assimilation was over by 15% at Tair 25 and 178% at Tair 40 and the guard
// returned true on every row -- it reads the solution and three bounds, not the
// deficit the solution was built from.
void test_closed_form_is_prescribed_temperature_only() {
  printf("closed form: prescribed-temperature path only\n");
  Drivers d;
  d.PPFD = 900.0;
  d.atm_vpd = 2.0;
  d.leaf_temp = 30.0;
  namespace cf = phylloptim::closed_form;

  // 1. Every solver refuses on the energy-balance path, naming the gate. The
  //    stop comes first, so the explicit-form arm is exercised without seating
  //    beta2 = 1/stem_c: what is under test is the refusal, not the arithmetic.
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0, true);
    l.CF77_lambda_ = 1.5e5;
    l.TF24_floor_lambda_o = 1.5e5;
    const double psi_soil = l.supply_psi_soil_scalar();
    ok(throws_with([&] { cf::solve(l, psi_soil); }, "use_energy_balance_"),
       "closed_form::solve refuses the energy-balance path, naming the gate");
    ok(throws_with([&] { cf::solve_exact_beta2(l, psi_soil); },
                   "use_energy_balance_"),
       "and so does solve_exact_beta2");
    ok(throws_with([&] { cf::solve_CF77(l, psi_soil); }, "use_energy_balance_"),
       "and solve_CF77");
    ok(throws_with([&] { cf::solve_TF24_floor(l, psi_soil); },
                   "use_energy_balance_"),
       "and solve_TF24_floor");
    // The production route raises the SAME stop: one function, so the two
    // cannot drift on what they refuse or why.
    l.set_model("TF24", "leaf", "closed");
    ok(throws_with([&] { l.optimise(); }, "up to 178%"),
       "optimise_closed raises the shared stop, measurement included");
  }

  // 2. Off the gate the deficit it divides by IS the air deficit, bit-exactly --
  //    the equality that makes `vpd_leaf_` a no-op substitution there, and whose
  //    absence from the tests is how the header drifted behind the live model.
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0, false);
    ok(l.vpd_leaf_ == l.atm_vpd_, "gate off: vpd_leaf_ IS atm_vpd_, bit-for-bit");
    const cf::Solution sol = cf::solve(l, l.supply_psi_soil_scalar());
    ok(std::isfinite(sol.assim) && std::isfinite(sol.transpiration),
       "and the closed form solves");
    // The closed form's conductance is the live model's expression
    // (`Leaf::stom_cond_CO2`), so at the SAME (E, ci) the two agree exactly.
    const double gc_live = l.atm_kpa_ * sol.transpiration *
                           phylloptim::kg_to_mol_h2o / l.vpd_leaf_ /
                           l.H2O_CO2_stom_diff_ratio_;
    ok(sol.stom_cond_CO2 == gc_live,
       "and its conductance is the live model's expression, bit-for-bit");
  }

  // 3. The guard's scope, pinned so its role cannot be misread: it is the ci/ca
  //    test plus the bracket, and it reads NO driver. Tripling the deficit after
  //    the solve -- the size of the error the energy balance would introduce --
  //    leaves its verdict unchanged, which is the whole reason the refusal above
  //    is a stop in the solvers rather than a case in the guard.
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0, false);
    const double psi_soil = l.supply_psi_soil_scalar();
    const cf::Solution sol = cf::solve(l, psi_soil);
    const bool verdict = cf::within_guard(l, sol);
    ok(verdict == (std::isfinite(sol.psi_leaf) && std::isfinite(sol.ci) &&
                   std::isfinite(sol.assim) && sol.ci / l.ca_ > 0.5 &&
                   sol.psi_leaf > psi_soil && sol.psi_leaf < l.psi_crit),
       "within_guard is exactly the ci/ca test plus the bracket");
    ok(verdict, "and the wet reference row passes it");
    l.vpd_leaf_ *= 3.0;
    l.atm_vpd_ *= 3.0;
    ok(cf::within_guard(l, sol) == verdict,
       "and a deficit wrong by 3x after the fact does not change its verdict");
  }
}

// `CF77_lambda_` is a caller INPUT -- the Cowan-Farquhar marginal value of water -- and
// nothing in the model writes it. ProfitMax reports the marginal cost its own
// normalisation implies in a SEPARATE field, because the two are not the same
// quantity: one is carbon per unit transpiration, the other per unit conductance.
//
// The load-bearing assertion is the third: a ProfitMax solve followed by a
// Cowan-Farquhar one must price water at the CALLER's number. Before the split it
// priced it at ProfitMax's, silently, because both are finite and plausible.
void test_profitmax_reports_an_emergent_lambda() {
  printf("ProfitMax reports an emergent lambda and leaves the input alone\n");
  Drivers d;
  d.PPFD = 1500.0;
  const double prescribed = 1.5e5;

  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  l.CF77_lambda_ = prescribed;
  l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::ProfitMax>();

  ok(l.CF77_lambda_ == prescribed, "the prescribed CF77_lambda_ survives untouched");
  ok(std::isfinite(l.lambda_emergent()) && l.lambda_emergent() > 0.0,
     "and ProfitMax reports the marginal cost its point implies");
  ok(l.lambda_emergent() != prescribed,
     "which is its own number, not the caller's");

  // The hazard the split exists to remove: solving one and then the other now
  // prices water at the caller's value, not at ProfitMax's.
  const double p_after = l.profit_psi_leaf_CF77(2.0, 0.5);
  phylloptim::Leaf fresh = make_single_leaf(d, 0.5);
  fresh.CF77_lambda_ = prescribed;
  const double p_clean = fresh.profit_psi_leaf_CF77(2.0, 0.5);
  ok(p_after == p_clean,
     "a Cowan-Farquhar solve after ProfitMax is bit-identical to a clean one");
}

// ⚠️ THE OBJECTIVE IS NOT UNIMODAL AND ITS MAXIMUM CAN BE AN ENDPOINT. This is
// the test that a bare Brent search fails: at a leaf hot enough that net
// assimilation is negative everywhere, the profit is highest at FULL CLOSURE and
// there is a local maximum out in the interior. Brent steps in from the bounds
// and cannot return an endpoint, so it used to report the local one -- an open
// stoma where the model says the leaf should be shut.
void test_profitmax_finds_a_closed_optimum() {
  printf("ProfitMax finds a boundary optimum, which Brent alone cannot\n");
  Drivers d;
  d.PPFD = 1500.0;
  d.leaf_temp = 50.0;

  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  l.use_thermal_cost_ = true;
  l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::ProfitMax>();

  // Reconstruct the objective on a coarse grid and find its global maximum
  // independently of the solver.
  const std::vector<double> curve = l.profitmax_curve(201);
  const std::size_t n = 201;
  std::size_t best = 0;
  for (std::size_t i = 1; i < n; ++i) {
    if (curve[4 * n + i] > curve[4 * n + best]) best = i;
  }
  printf("    grid argmax at psi = %.4f (index %zu of %zu), solver returned %.4f\n",
         curve[best], best, n, l.opt_psi_leaf_);

  // profitmax_curve re-prepares, so re-solve before reading the operating point.
  l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  const double step = curve[1] - curve[0];
  ok(std::abs(l.opt_psi_leaf_ - curve[best]) < 3.0 * step,
     "the solver lands on the objective's GLOBAL maximum, not a local one");

  // And the profit it reports is at least the grid's best.
  ok(l.profit_ >= curve[4 * n + best] - 1e-9,
     "at no lower profit than the grid's best");
}

// The two PRODUCT objectives go through the same additive machinery as everything
// else, by storing their cost in log form: `C = -log g`, so `h(A) - C` is
// `log(A*g)`. This asserts the consequence as an IDENTITY at interior points rather
// than as a residual at one argmax, because the residual-at-the-argmax version
// passes on any derivative that merely happens to vanish near there:
//
//     (log P)' * P == P',   P = A*g      <=>   (A'/A + g'/g) * A*g == A'g + A*g'
//
// i.e. the framework's first-order condition is the product rule divided through by
// the strictly positive `A*g`, so the two have identical roots and putting the link
// on the benefit alone drops nothing.
//
// ⚠️ NO FINITE DIFFERENCE HERE, AND THAT IS THE POINT. The first version of this
// test differenced `P` and compared. It passed on macOS at 1.08e-05 and FAILED on
// both Linux compilers at 2.51e-04, because differencing through psi_leaf_to_ci's
// root-find has a noise floor of tol/h (hazard 15) and `tol` is set by libm -- so
// the floor is PLATFORM-DEPENDENT and a tolerance measured on one platform does not
// transfer. Every term below is analytic instead: `A'` is a zero-cost JS22
// derivative, `g` and `g'` are the reduction factor and its slope, and the test is
// then pure arithmetic on quantities both sides already computed.
//
// ⚠️ NORMALISED BY THE TERM SCALE, NOT BY THE RESULT. `A'g + A*g'` is a difference
// of two similar magnitudes that cancels to exactly zero at the argmax, so a
// relative error against the result diverges there for a CORRECT implementation.
// Dividing by max(|A'g|, |A*g'|) asks the question that has an answer: are the two
// ways of combining the same terms the same number, to roundoff.
void test_product_link_is_the_product_rule() {
  printf("the product curves' log link is the product rule over A*g\n");

  struct Row { double psi_soil, PPFD, atm_vpd; };
  const Row rows[] = {{0.5, 1500.0, 2.0}, {1.0, 900.0, 2.0},
                      {2.0, 900.0, 3.0}, {1.0, 300.0, 1.0}};
  const int kJS22 = static_cast<int>(phylloptim::Leaf::CostCurve::JS22);
  const int curves[] = {static_cast<int>(phylloptim::Leaf::CostCurve::SOX),
                        static_cast<int>(phylloptim::Leaf::CostCurve::JW26)};

  double worst = 0.0;
  int checked = 0;
  for (const int curve : curves) {
    const bool is_sox = curve == static_cast<int>(phylloptim::Leaf::CostCurve::SOX);
    for (const Row &r : rows) {
      Drivers d;
      d.PPFD = r.PPFD;
      d.atm_vpd = r.atm_vpd;
      phylloptim::Leaf l = make_single_leaf(d, r.psi_soil);
      // A zero unit cost makes the JS22 arm's C' exactly 0, so under the identity
      // link its derivative IS dA/dpsi. Safe as a bare write: JS22_gamma enters
      // one cost expression and no derived state (unlike the curve or capacity
      // traits, which hazard 10 covers).
      l.JS22_gamma = 0.0;
      const double lo = r.psi_soil, hi = l.psi_crit;
      for (int k = 1; k <= 9; ++k) {
        const double psi = lo + (hi - lo) * (0.1 * double(k));
        // P and A first: the derivative calls below re-seat the leaf's state.
        const double P = is_sox ? l.profit_psi_leaf_SOX(psi, r.psi_soil)
                                : l.profit_psi_leaf_JW26(psi, r.psi_soil);
        const double A = l.assim_colimited_;
        const double g = is_sox ? l.sox_reduction(psi) : l.jw26_reduction(psi);
        const double gp =
            is_sox ? l.sox_reduction_deriv(psi) : l.jw26_reduction_deriv(psi);
        // The premise of the log form. Where it fails nothing is claimed: at full
        // closure A = -R_d < 0 and log A does not exist, which is precisely why the
        // OPTIMISER searches the product instead of its logarithm.
        if (!(P > 0.0) || !(g > 0.0)) {
          continue;
        }
        l.set_model(static_cast<phylloptim::Leaf::CostCurve>(curve), false);
        const std::vector<double> dlog = l.dprofit_dpsi_leaf_checked(psi);
        l.set_model(static_cast<phylloptim::Leaf::CostCurve>(kJS22), false);
        const std::vector<double> dA = l.dprofit_dpsi_leaf_checked(psi);
        if (!std::isfinite(dlog[0]) || !std::isfinite(dA[0])) {
          continue;
        }
        const double lhs = dlog[0] * P;          // (log P)' * P
        const double rhs = dA[0] * g + A * gp;   // the product rule
        const double scale = std::max(std::fabs(dA[0] * g), std::fabs(A * gp));
        const double rel = std::fabs(lhs - rhs) / std::max(scale, 1e-30);
        if (rel > worst) {
          worst = rel;
        }
        ++checked;
      }
    }
  }
  // ⚠️ The count is asserted because the two `continue`s above can empty this test
  // silently -- a curve returning NaN everywhere would leave `worst` at 0 and pass
  // (hazard 16).
  printf("    %d interior points | worst departure from (logP)'*P == A'g + Ag': "
         "%.3e\n", checked, worst);
  ok(checked >= 40, "the identity was actually evaluated, on both curves");
  ok(worst < 1e-12, "the log-link derivative IS the product rule divided by A*g");
}

void test_shadow_cost() {
  printf("carbon and shadow_cost split each curve's cost into carbon and price\n");
  // WHAT THIS PINS. `hydraulic_cost_` is the whole cost a difference objective
  // subtracted. Each curve's definition states the part of it that is carbon
  // spent (Curve<K>::carbon_cost); `shadow_cost()` is the rest, a price, and
  // `carbon()` is A less the carbon part, what a consumer growing a plant on
  // this leaf books, never `profit`. The carbon part is TF24's whole cost and
  // TF24_floor's TF24 half; nothing on CF77, JS22 and CMax, whose costs are
  // prices; and the product, ratio and normalised curves have no cost in
  // carbon units to split, so their shadow cost is NA and their carbon is A.
  //
  // ⚠️ BIT EQUALITY against the reported numbers. Where the price reads only
  // the flux, both read the reported objective and E (so a caller that checks
  // `shadow_cost == lambda_o * E` gets it exactly), and carbon equals A less the
  // carbon part to rounding: the cost was subtracted and the price added back.
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  int compared = 0, per_curve = 0;
  for (double ppfd : {300.0, 900.0, 1500.0}) {
    for (double vpd : {1.0, 2.0, 3.0}) {
      for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
        d.PPFD = ppfd;
        d.atm_vpd = vpd;

        for (double lo : {0.0, 1.5e4, 1.5e5}) {
          phylloptim::Leaf l = make_single_leaf(d, psi_soil);
          l.TF24_floor_lambda_o = lo;
          l.set_model(CC::TF24_floor, false);
          l.optimise();
          if (!std::isfinite(l.transpiration_)) continue;
          ok(l.shadow_cost() == lo * l.transpiration_,
             "TF24_floor's shadow cost is lambda_o * E exactly");
          ok(l.carbon() == l.profit_ + l.shadow_cost(),
             "TF24_floor's carbon is profit plus the price, exactly");
          near(l.carbon(),
               l.assim_colimited_ - l.hydraulic_cost_TF24_kernel(l.opt_psi_leaf_), 1e-12,
               "TF24_floor's carbon is A less TF24's cost");
          ok((lo == 0.0) == (l.shadow_cost() == 0.0),
             "the price vanishes exactly when lambda_o does");
          ok(l.carbon() >= l.profit_, "the carbon kept is never below the objective");
          ++compared;
        }

        for (int c = 0; c < phylloptim::Leaf::n_cost_curves; ++c) {
          const CC curve = static_cast<CC>(c);
          const std::string name = phylloptim::Leaf::curve_name(c);
          phylloptim::Leaf l = make_single_leaf(d, psi_soil);
          l.CF77_lambda_ = 1.5e5;
          l.set_model(curve, false);
          l.optimise();
          if (!std::isfinite(l.assim_colimited_)) continue;
          if (curve == CC::TF24) {
            ok(l.shadow_cost() == 0.0, "TF24's shadow cost is exactly zero");
            ok(l.carbon() == l.profit_, "TF24's carbon IS its objective");
          } else if (curve == CC::TF24_floor) {
            continue;  // above, at a price
          } else if (curve == CC::CF77) {
            ok(l.shadow_cost() == l.hydraulic_cost_,
               "CF77's cost is all price: its shadow cost is the whole cost");
            ok(l.carbon() == l.profit_ + l.shadow_cost(),
               "CF77's carbon is profit plus the price, exactly");
            near(l.carbon(), l.assim_colimited_, 1e-12, "CF77's carbon is A");
          } else if (curve == CC::JS22 || curve == CC::CMax) {
            // Read off the reported A and profit, so a crown average holds it.
            ok(l.shadow_cost() == l.assim_colimited_ - l.profit_,
               name + "'s shadow cost is A less profit, exactly");
            near(l.shadow_cost(), l.hydraulic_cost_, 1e-12,
                 name + "'s cost is all price: its shadow cost is the whole cost");
            ok(l.carbon() == l.assim_colimited_, name + "'s carbon is A");
          } else {
            ok(std::isnan(l.shadow_cost()),
               name + " has no cost in carbon to split: its shadow cost is NA");
            ok(l.carbon() == l.assim_colimited_, name + "'s carbon is A");
          }
          ++per_curve;
        }
      }
    }
  }
  printf("  compared %d priced TF24_floor points and %d points across the other "
         "curves\n", compared, per_curve);
  ok(per_curve >= 8 * 30, "every other curve was compared on most of the grid");
}

void test_lambda_emergent_reported() {
  printf("operating_point_values() reports the SEATED curve's lambda\n");
  // R's `lambda` is this slot: the seated curve's own (dC/dpsi)/(dE/dpsi), not
  // TF24's whatever curve ran, which on any other curve answers a question
  // nobody asked.
  //
  // ⚠️ CF77 IS THE CASE THAT PINS IT, because its emergent lambda is KNOWN: the
  // curve prices water at a prescribed constant, so `lambda_emergent` must come
  // back as that constant. Nothing else here has an answer fixed independently of
  // the solve.
  Drivers d;
  const double price = 1.5e5;
  int compared = 0;
  for (double ppfd : {300.0, 900.0, 1500.0}) {
    for (double psi_soil : {0.5, 1.0, 2.0}) {
      d.PPFD = ppfd;

      phylloptim::Leaf c = make_single_leaf(d, psi_soil);
      c.CF77_lambda_ = price;
      c.set_model(phylloptim::Leaf::CostCurve::CF77, false);
      c.optimise();
      if (std::isfinite(c.transpiration_) && c.transpiration_ > 0.0) {
        ok(std::fabs(c.lambda_emergent() - price) <= 1e-6 * price,
           "CF77's emergent lambda is the price it was given");
        ++compared;
      }

      // TF24's is its bound per-curve lambda at the optimum, exactly
      phylloptim::Leaf t = make_single_leaf(d, psi_soil);
      t.set_model(phylloptim::Leaf::CostCurve::TF24, false);
      t.optimise();
      ok(t.lambda_emergent() == t.lambda_TF24(t.opt_psi_leaf_),
         "on TF24 the seated lambda is lambda_TF24 at the optimum exactly");

      // and TF24_floor's is its hydraulic half plus the floor -- the correction a
      // caller used to have to make by hand
      phylloptim::Leaf f = make_single_leaf(d, psi_soil);
      f.TF24_floor_lambda_o = price;
      f.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
      f.optimise();
      near(f.lambda_emergent(), f.lambda_TF24(f.opt_psi_leaf_) + price, 1e-7,
           "TF24_floor's emergent lambda is the hydraulic half plus the floor");

      // the reported vector carries the same number the accessor does, in the
      // slot R names `lambda`, and the residual flux last
      const std::vector<double> ov = f.operating_point_values();
      ok(ov[10] == f.lambda_emergent(),
         "operating_point_values() reports lambda_emergent as `lambda`");
      ok(ov[14] == f.residual_flux_, "and the residual flux after the shadow cost");
      ok(ov.size() == 16 && ov.back() == f.carbon(), "and carbon last");
    }
  }
  printf("  %d CF77 points with a known emergent lambda\n", compared);
}

void test_tf24_floor_reduces() {
  printf("TF24_floor reduces to TF24 and to CF77, bit-for-bit\n");
  // ⚠️ BIT EQUALITY IS THE RIGHT TEST HERE AND A TOLERANCE WOULD BE THE WRONG
  // ONE. The reductions are exact by construction rather than by convergence:
  // each term is the parent's own expression, so zeroing one parameter adds an
  // exact zero to the other curve's exact value -- see
  // `Leaf::hydraulic_cost_TF24_floor`, which also records why fused multiply-add
  // does not break it. A tolerance would pass on a rewrite that lost the
  // property, which is the only thing this test exists to catch.
  Drivers d;
  int compared = 0, cf77_pinned = 0;
  for (double ppfd : {300.0, 900.0, 1500.0}) {
    for (double vpd : {1.0, 2.0, 3.0}) {
      for (double psi_soil : {0.5, 1.0, 2.0, 3.0, 4.0}) {
        d.PPFD = ppfd;
        d.atm_vpd = vpd;

        // --- lambda_o = 0 is TF24, at TF24's own traits --------------------
        // ⚠️ THE ONE-RESTRICTION REDUCTION, and the reason this curve has no
        // hydraulic parameter of its own. Nothing is set on either leaf but the
        // price: if the two disagree, they disagree about `lambda_o` alone.
        {
          phylloptim::Leaf a = make_single_leaf(d, psi_soil);
          a.TF24_floor_lambda_o = 0.0;
          a.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
          a.optimise();

          phylloptim::Leaf b = make_single_leaf(d, psi_soil);
          b.set_model(phylloptim::Leaf::CostCurve::TF24, false);
          b.optimise();

          const std::vector<double> va = reported(a), vb = reported(b);
          for (std::size_t i = 0; i < va.size(); ++i) {
            ok(same_bits(va[i], vb[i]),
               std::string("TF24_floor(lambda_o=0) == TF24 in ") +
                   reported_names[i] + " at psi_soil=" +
                   std::to_string(psi_soil) + " ppfd=" + std::to_string(ppfd) +
                   " vpd=" + std::to_string(vpd));
          }
        }

        // --- TF24_cost_scale = 0 is CF77 -----------------------------------
        // The hydraulic half is scaled by `TF24_cost_scale`, so zeroing that
        // trait is what leaves the price alone. `pow(1-f, beta2)` is finite
        // inside the bracket, so `0 * pow(...)` is an exact zero rather than a
        // NaN -- which is why this reduction is available at all.
        //
        // ⚠️ EXCEPT WHERE CF77 SITS AT psi_crit. CF77's domain ends there
        // (cost_domain_ends_at_psi_crit) and TF24_floor's does not (#183), so
        // where CF77 is pinned at its bound TF24_floor goes on past it, to a
        // profit at least as high. Elsewhere both solve on [psi_soil, psi_crit]
        // and the reduction is bit-for-bit.
        {
          phylloptim::Leaf a = make_single_leaf(d, psi_soil);
          a.TF24_cost_scale = 0.0;
          a.TF24_floor_lambda_o = kFloorLambdaO;
          a.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
          a.optimise();

          phylloptim::Leaf b = make_single_leaf(d, psi_soil);
          b.CF77_lambda_ = kFloorLambdaO;
          b.set_model(phylloptim::Leaf::CostCurve::CF77, false);
          b.optimise();

          const std::string at = " at psi_soil=" + std::to_string(psi_soil) +
                                 " ppfd=" + std::to_string(ppfd) +
                                 " vpd=" + std::to_string(vpd);
          if (b.opt_psi_leaf_ == b.psi_crit) {
            ++cf77_pinned;
            ok(a.opt_psi_leaf_ >= b.opt_psi_leaf_ && a.profit_ >= b.profit_,
               "where CF77 is pinned at psi_crit, TF24_floor(TF24_cost_scale=0) "
               "goes past it to a profit at least as high" + at);
          } else {
            const std::vector<double> va = reported(a), vb = reported(b);
            for (std::size_t i = 0; i < va.size(); ++i) {
              ok(same_bits(va[i], vb[i]),
                 std::string("TF24_floor(TF24_cost_scale=0) == CF77 in ") +
                     reported_names[i] + at);
            }
          }
        }
        ++compared;
      }
    }
  }
  // ⚠️ Hazard 16: an assertion that never runs passes. Both reductions are
  // inside two loops and a refused row would silently skip one.
  printf("    %d driver rows, both reductions, %d fields each; CF77 pinned at "
         "psi_crit on %d\n", compared,
         int(sizeof(reported_names) / sizeof(reported_names[0])), cf77_pinned);
  ok(compared == 45, "every driver row was compared");
  ok(cf77_pinned < compared, "and the CF77 reduction was compared bit-for-bit somewhere");
}

void test_tf24_floor_wet_end_price() {
  printf("TF24_floor's marginal cost of water does not vanish at the wet end\n");
  Drivers d;
  phylloptim::Leaf l = make_single_leaf(d, 1.0);
  l.TF24_floor_lambda_o = kFloorLambdaO;

  // ⚠️ THE LIMIT IS AT ZERO POTENTIAL, NOT AT ZERO DROP, and that is the one
  // behavioural difference from the JS22-based version of this curve. The
  // hydraulic half reads the ABSOLUTE potential, so it vanishes as psi -> 0 --
  // not as psi -> psi_soil, where TF24's own cost is already non-zero.
  double prev_excess = 0.0;
  for (int k = 0; k < 4; ++k) {
    const double psi = 1e-2 / std::pow(2.0, k);
    const double excess = l.lambda_TF24_floor(psi) - l.TF24_floor_lambda_o;
    ok(excess > 0.0, "the excess over lambda_o is positive");
    if (k > 0) {
      // lambda_TF24 ~ psi^(stem_c*beta2 - 1), so halving psi divides the excess
      // by 2^(stem_c*beta2 - 1). Asserting the EXPONENT rather than a ratio of 2
      // is what makes this a statement about the curve rather than about a step.
      const double n = l.stem_c * l.TF24_beta2 - 1.0;
      near(prev_excess / excess, std::pow(2.0, n), 1e-2,
           "halving psi divides the excess by 2^(stem_c*beta2 - 1)");
    }
    prev_excess = excess;
  }
  ok(l.stem_c * l.TF24_beta2 > 1.0,
     "and the wet-end limit exists at all, which needs stem_c*beta2 > 1");

  // ⚠️ THE CONTRAST IS THE POINT. At the same potential every conductance-loss
  // curve here is heading for zero, which is what "water is free when abundant"
  // means as a number rather than as an argument.
  const double psi = 1e-6;
  printf("    at psi = %.0e MPa: TF24_floor %.6g, TF24 %.3e, JS22 %.3e\n", psi,
         l.lambda_TF24_floor(psi), l.lambda_TF24(psi), l.lambda_JS22(psi, 0.0));
  ok(l.lambda_TF24(psi) < 1e-3 * l.TF24_floor_lambda_o,
     "TF24's price at the same potential is negligible beside lambda_o");
  near(l.lambda_TF24_floor(psi) / l.TF24_floor_lambda_o, 1.0, 1e-6,
     "so TF24_floor's is lambda_o to six figures");
}

void test_tf24_floor_emergent_lambda() {
  printf("TF24_floor's emergent lambda is dC/dE at the solved point\n");
  // Against a finite difference of the COST over the TRANSPIRATION, not against
  // a second copy of `lambda_TF24_floor` -- which would only assert that the
  // dispatch table reaches the function it reaches. This tests the definition.
  Drivers d;
  d.PPFD = 1500.0;
  const double h = 1e-3;   // the step `test_every_curve_reports_an_emergent_lambda` uses
  int interior = 0;

  for (double psi_soil : {0.5, 1.0, 2.0}) {
    phylloptim::Leaf l = make_single_leaf(d, psi_soil);
    l.TF24_floor_lambda_o = kFloorLambdaO;
    l.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
    l.optimise();
    const double psi = l.opt_psi_leaf_;
    const double got = l.lambda_emergent();
    if (psi <= psi_soil + h || psi >= l.psi_crit - h) {
      continue;   // a pinned optimum: dC/dE is defined but is not what pinned it
    }
    ++interior;

    // ⚠️ THROUGH `profit_psi_leaf_TF24_floor`, NOT `hydraulic_cost_TF24_floor`. The cost
    // function reads the memoised transpiration but does not SEAT the leaf, so
    // reading `transpiration_` after it returns gives whatever the last solve
    // left there -- the same E at both ends, a zero denominator, and a ratio of
    // exactly 0 that looks like a broken lambda. The profit function seats it.
    l.profit_psi_leaf_TF24_floor(psi + h, psi_soil);
    const double C1 = l.hydraulic_cost_, E1 = l.transpiration_;
    l.profit_psi_leaf_TF24_floor(psi - h, psi_soil);
    const double C0 = l.hydraulic_cost_, E0 = l.transpiration_;
    near(got / ((C1 - C0) / (E1 - E0)), 1.0, 1e-4,
         "TF24_floor's emergent lambda is dC/dE at psi_soil=" +
             std::to_string(psi_soil));

    // And it exceeds lambda_o by the hydraulic term, the decomposition read back
    // off the solved point, to the stem table's slope accuracy: the floor's
    // part reads E' from the table, TF24's from the closed form.
    ok(got > l.TF24_floor_lambda_o,
       "the emergent lambda sits above the wet-end floor at psi_soil=" +
           std::to_string(psi_soil));
    near(got, l.lambda_TF24(psi) + l.TF24_floor_lambda_o, 1e-7,
       "and it IS TF24's own lambda plus the floor at psi_soil=" +
           std::to_string(psi_soil));
  }
  // Hazard 16 again: print what was actually exercised.
  printf("    %d of 3 rows interior and compared\n", interior);
  ok(interior >= 2, "the identity was evaluated on at least two rows");
}

void test_tf24_floor_kmax_vpd_invariance() {
  printf("TF24_floor sits between TF24 and CF77 under kmax/D rescaling\n");
  // ⚠️ A PREDICTION, NOT A REGRESSION CHECK. Scaling `kmax` and the vapour
  // deficit together leaves the stomatal conductance -- and so ci, A and dA/dpsi
  // -- untouched at every potential, because `gc ~ kmax*G(psi)/D`. So a cost
  // whose MARGINAL value depends on psi alone is invariant, and one that prices
  // the FLUX is not: `dC/dpsi` is TF24's `|f'|` form for the first part and
  // `lambda_o*kmax*f(psi)` for the second. TF24_floor carries both, so it must move
  // by less than CF77 and by more than TF24, and by more as lambda_o rises.
  const double factor = 4.0;
  Drivers lo, hi;
  lo.PPFD = hi.PPFD = 1500.0;
  // Both deficits below saturation at 25 C (3.17 kPa).
  lo.atm_vpd = 0.75;
  hi.atm_vpd = lo.atm_vpd * factor;
  hi.K_s = lo.K_s * factor;
  const double psi_soil = 1.0;

  auto solved = [&](const Drivers &d, phylloptim::Leaf::CostCurve curve,
                    double lambda_o) {
    phylloptim::Leaf l = make_single_leaf(d, psi_soil);
    l.CF77_lambda_ = kFloorLambdaO;
    l.TF24_floor_lambda_o = lambda_o;
    l.set_model(curve, false);
    l.optimise();
    return l.opt_psi_leaf_;
  };
  auto shift = [&](phylloptim::Leaf::CostCurve curve, double lambda_o) {
    const double a = solved(lo, curve, lambda_o);
    const double b = solved(hi, curve, lambda_o);
    return std::abs(b - a) / a;
  };

  using Curve = phylloptim::Leaf::CostCurve;
  const double tf24 = shift(Curve::TF24, 0.0);
  const double cf77 = shift(Curve::CF77, 0.0);
  printf("    kmax and D both x%.0f:  TF24 %.3e  CF77 %.3e\n", factor, tf24,
         cf77);
  ok(tf24 < 1e-12, "TF24's optimum is invariant");
  ok(cf77 > 0.1, "CF77's is not, because it prices the flux");

  // The interpolation, and that it is monotone in lambda_o.
  double previous = 0.0;
  for (double lambda_o : {1.5e3, 1.5e4, 1.5e5}) {
    const double s = shift(Curve::TF24_floor, lambda_o);
    printf("    TF24_floor at lambda_o = %.1e: %.3e\n", lambda_o, s);
    ok(s > tf24, "TF24_floor moves more than TF24 at lambda_o=" +
                     std::to_string(lambda_o));
    ok(s < cf77, "and less than CF77 at lambda_o=" + std::to_string(lambda_o));
    ok(s > previous, "and by more as lambda_o rises");
    previous = s;
  }
}

void test_tf24_floor_defaults_to_tf24() {
  printf("TF24_floor's price defaults to 0, where it is TF24; a cleared one is refused\n");
  Drivers d;
  {
    // The default is plant's, so a leaf that seats TF24_floor and sets nothing
    // solves TF24, bit for bit (test_tf24_floor_reduces says why exactly).
    phylloptim::Leaf a = make_single_leaf(d, 1.0);
    ok(a.TF24_floor_lambda_o == 0.0, "a fresh leaf's lambda_o is 0");
    a.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
    a.optimise();
    phylloptim::Leaf b = make_single_leaf(d, 1.0);
    b.set_model(phylloptim::Leaf::CostCurve::TF24, false);
    b.optimise();
    const std::vector<double> va = reported(a), vb = reported(b);
    for (std::size_t i = 0; i < va.size(); ++i) {
      ok(same_bits(va[i], vb[i]),
         std::string("the default TF24_floor is TF24 in ") + reported_names[i]);
    }
  }
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.TF24_floor_lambda_o = phylloptim::util::na_value;
    l.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
    ok(throws_with([&] { l.optimise(); },
                   "needs TF24_floor_lambda_o set"),
       "a cleared lambda_o is refused rather than optimised as a NaN objective");
    // It has to say where a price can be SET FROM, since the field is not
    // reachable through leaf_traits().
    ok(throws_with([&] { l.optimise(); }, "leaf_solve()"),
       "and names a route that can supply one");
  }
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.TF24_floor_lambda_o = -1.0;
    l.set_model(phylloptim::Leaf::CostCurve::TF24_floor, false);
    ok(throws_with([&] { l.optimise(); }, "non-negative"),
       "a negative price -- which would PAY the leaf to transpire -- is refused");
  }
}

void test_tf24_floor_closed_form() {
  printf("TF24_floor's closed form: the reduction, the coverage, and the tail\n");
  Drivers d;
  d.PPFD = 1500.0;

  // --- at TF24_cost_scale = 0 this curve is CF77, and the two closed arms are
  // compared. ⚠️ NOT BIT-FOR-BIT, AND THE DIFFERENCE IS DELIBERATE.
  // `solve_CF77` takes no Newton step at all -- with lambda constant there is no
  // circularity to break -- and stops at the spline inversion of the supply. This
  // arm still takes its step, and with a zeroed hydraulic term there is nothing in
  // lambda for it to correct, so what it does instead is drive the
  // supply-minus-demand residual to zero from a start clamped strictly inside the
  // bracket. So the comparison worth making is against the EXACT answer, and the
  // finding is that the extra step does not hurt.
  int both_served = 0, different_decision = 0;
  double worst_gap = 0.0;
  for (double psi_soil : {0.5, 1.0, 2.0}) {
    phylloptim::Leaf exact = make_single_leaf(d, psi_soil);
    exact.TF24_cost_scale = 0.0;
    exact.TF24_floor_lambda_o = 1.5e5;
    exact.set_model("TF24_floor", "leaf", "exact");
    exact.optimise();

    phylloptim::Leaf a = make_single_leaf(d, psi_soil);
    a.TF24_cost_scale = 0.0;
    a.TF24_floor_lambda_o = 1.5e5;
    a.set_model("TF24_floor", "leaf", "closed");
    a.optimise();

    phylloptim::Leaf b = make_single_leaf(d, psi_soil);
    b.CF77_lambda_ = 1.5e5;
    b.set_model("CF77", "leaf", "closed");
    b.optimise();

    // The EXACT solves must agree bit-for-bit -- that is the reduction, on the
    // method both curves share, and it is asserted without a tolerance.
    phylloptim::Leaf cf_exact = make_single_leaf(d, psi_soil);
    cf_exact.CF77_lambda_ = 1.5e5;
    cf_exact.set_model("CF77", "leaf", "exact");
    cf_exact.optimise();
    ok(same_bits(exact.opt_psi_leaf_, cf_exact.opt_psi_leaf_),
       "exact TF24_floor(TF24_cost_scale=0) IS exact CF77 at psi_soil=" +
           std::to_string(psi_soil));

    if (a.last_solve_fell_back_ != b.last_solve_fell_back_) {
      ++different_decision;
      continue;
    }
    if (a.last_solve_fell_back_) {
      continue;
    }
    ++both_served;
    const double ea = std::abs(a.opt_psi_leaf_ - exact.opt_psi_leaf_);
    const double eb = std::abs(b.opt_psi_leaf_ - exact.opt_psi_leaf_);
    ok(ea <= eb * 1.05,
       "the TF24_floor arm is no further from the exact optimum than CF77's at "
       "psi_soil=" + std::to_string(psi_soil));
    worst_gap = std::max(worst_gap,
                         std::abs(a.opt_psi_leaf_ - b.opt_psi_leaf_) /
                             b.opt_psi_leaf_);
  }
  printf("    scale=0: %d rows served by both (worst |dpsi*| between arms %.3e), "
         "%d where only one was served\n", both_served, worst_gap,
         different_decision);
  ok(both_served >= 2, "the scale = 0 comparison was actually exercised");

  // --- coverage and the worst served error at real parameters -----------------
  // ⚠️ THIS BOUND IS A MEASUREMENT, NOT A SPECIFICATION. What limits the accuracy
  // is the USO collapse (see closed_form.hpp), and on this curve the leading order
  // additionally discards the whole TF24 term -- so the tail is worse than CF77's
  // and the shared ci/ca guard does not see it. The assertion exists to notice
  // DRIFT; read the printed numbers, not the pass.
  int served = 0, fell_back = 0;
  double worst = 0.0;
  for (double scale : {0.5, 7.5}) {
    for (double psi_soil : {0.5, 1.0, 1.5, 2.0}) {
      for (double vpd : {1.0, 2.0}) {
        Drivers dd = d;
        dd.atm_vpd = vpd;
        phylloptim::Leaf exact = make_single_leaf(dd, psi_soil);
        exact.TF24_cost_scale = scale;
        exact.TF24_floor_lambda_o = 1.5e5;
        exact.set_model("TF24_floor", "leaf", "exact");
        exact.optimise();

        phylloptim::Leaf closed = make_single_leaf(dd, psi_soil);
        closed.TF24_cost_scale = scale;
        closed.TF24_floor_lambda_o = 1.5e5;
        closed.set_model("TF24_floor", "leaf", "closed");
        closed.optimise();

        if (closed.last_solve_fell_back_) {
          ++fell_back;
          // ⚠️ A FALLBACK IS THE WHOLE EXACT SOLVE, so it must be bit-identical
          // rather than close. This is what makes the fallback fraction a cost
          // rather than an accuracy question.
          ok(same_bits(closed.opt_psi_leaf_, exact.opt_psi_leaf_) &&
                 same_bits(closed.profit_, exact.profit_),
             "a fallback IS the exact solve, bit-for-bit");
          continue;
        }
        ++served;
        ok(closed.opt_psi_leaf_ > psi_soil && closed.opt_psi_leaf_ < closed.psi_crit,
           "a served row is strictly inside the bracket");
        worst = std::max(worst,
                         std::abs(closed.opt_psi_leaf_ - exact.opt_psi_leaf_) /
                             exact.opt_psi_leaf_);
      }
    }
  }
  printf("    %d served, %d fell back (phi %.2f) | worst served |dpsi*|/psi* %.3e\n",
         served, fell_back, double(fell_back) / double(served + fell_back), worst);
  ok(served >= 2, "the closed form served at least two rows");
  ok(worst < 0.6, "and no served row is worse than the measured tail");

  // --- and where it is refused -------------------------------------------------
  // At lambda_o = 0 the curve is TF24, whose wet-end lambda is zero rather than a
  // constant: `xi = sqrt(Q/lambda)` diverges and there is nothing to expand about.
  // Refused, naming which of the two curves it has become -- and TF24 has its own
  // start, so the message can point somewhere.
  phylloptim::Leaf zero = make_single_leaf(d, 1.0);
  zero.TF24_floor_lambda_o = 0.0;
  zero.set_model("TF24_floor", "leaf", "closed");
  ok(throws_with([&] { zero.optimise(); }, "needs TF24_floor_lambda_o > 0"),
     "the closed form refuses lambda_o = 0");
  ok(throws_with([&] { zero.optimise(); }, "Seat TF24 instead"),
     "and points at the curve that does have a start");
}

// ===========================================================================
// Least-cost / cost-per-unit-carbon (Prentice et al. 2014)
// ---------------------------------------------------------------------------
// THIS IS THE ONLY CURVE HERE WITH A PUBLISHED ANALYTICAL SOLUTION, so it is the
// only one whose solver can be checked against something outside this package.
// Prentice et al. minimise `(c_E*E + c_V*Vcmax)/A` and solve the first-order
// condition in `ci` directly -- where the `A` the log link puts into lambda
// cancels against the `A` the ratio divides by. Writing
//
//     xi = sqrt(beta*(Km + Gamma*)/1.6),     beta = c_V/c_E,
//
// their condition is
//
//     (ci - Gamma*) / (ca - ci)  =  xi / sqrt(D)                        (FOC)
//
// and the form they report is (FOC) with `Gamma*` dropped from the numerator:
//
//     chi = ci/ca = xi / (xi + sqrt(D)).                          (published)
//
// ⚠️ THE TWO ARE NOT THE SAME EQUATION AND THE DIFFERENCE IS NOT SMALL. Solving
// (FOC) for chi gives `chi = (Gamma*/ca + r)/(1 + r)` with `r = xi/sqrt(D)`, so
// the published form understates chi by `(Gamma*/ca)/(1 + r)` -- some 4% at
// D = 1 kPa and 7% at D = 3 kPa on this package's defaults, growing with deficit
// because `r` falls as D rises. This test checks against (FOC), which is the
// model; it then checks the gap to the published form against that expression,
// which is what establishes that the gap is the approximation and not a defect.
//
// ⚠️ THREE THINGS HAVE TO BE ARRANGED BEFORE THE COMPARISON MEANS ANYTHING, and
// each is a real difference between this package and the published derivation.
//
//   1. THE ASSIMILATION MUST BE THE SAME FUNCTION. Prentice et al. use the
//      Rubisco-limited rate; this package co-limits it against electron
//      transport through a hyperbolic minimum and subtracts `R_d_`. Zeroing
//      `R_d_` removes the second difference, and driving `electron_transport_`
//      up removes the first -- the hyperbolic minimum tends to the Rubisco
//      branch exactly as the other limit goes to infinity. At this package's
//      DEFAULT co-limited settings the agreement is a few percent rather than
//      exact, which the last block below measures rather than hides.
//   2. THE OPTIMUM MUST BE INTERIOR. There is no vulnerability curve in
//      least-cost theory; here `E` is delivered through the supply relation, so
//      a solve sitting on the wet bound or on `psi_crit` is answering a
//      different question. `K_s` is raised so that it does not, and interiority
//      is asserted rather than assumed.
//   3. `xi` MUST CARRY THIS PACKAGE'S DIFFUSION RATIO, NOT THE PUBLISHED 1.6.
//      The 1.6 in `xi` is the H2O:CO2 stomatal diffusion ratio, and this package
//      uses 1.67 (`H2O_CO2_stom_diff_ratio`, a convention rather than a property
//      of the leaf). Testing against `xi` built on 1.6 while the solver diffuses
//      at 1.67 leaves a 0.6% error that looks like a solver defect and is not.
//      Both settings are run below, and the check holds under each.
//
// ⚠️ THE VISCOSITY TERM IS ABSENT. Their `xi` carries `eta*`, water viscosity
// relative to 25 C, which this package does not model, so every case here is run
// at 25 C where `eta*` is 1 by definition. Away from 25 C this package's beta is
// their `beta/eta*`.
void test_leastcost_matches_prentice_closed_form() {
  printf("LeastCost against Prentice et al. (2014)'s analytical solution\n");
  const double beta = 146.0;   // the published unit-cost ratio
  const double vpds[3] = {1.0, 2.0, 3.0};

  // `K_s` raised for reason 2. Checked across 10-1000 while this test was
  // written: the recovered chi is identical to six figures over that whole
  // range, so the supply relation is genuinely out of the comparison rather
  // than merely loose.
  auto solved = [&](double vpd_kpa, double b, double ratio, bool rubisco_only) {
    Drivers d;
    d.K_s = 200.0;
    d.PPFD = 3000.0;
    d.leaf_temp = 25.0;
    d.ca = 40.0;
    d.atm_vpd = vpd_kpa;
    phylloptim::Leaf l = make_single_leaf(d, 0.1);
    if (rubisco_only) {
      // Reason 1, both halves. `R_d_` and `electron_transport_` rather than
      // `R_d_25` and `jmax_25`, because `set_physiology` has already run inside
      // `make_single_leaf` and these are the derived fields the kernel reads.
      l.R_d_ = 0.0;
      l.electron_transport_ = 1.0e9;
    }
    l.H2O_CO2_stom_diff_ratio_ = ratio;
    l.LeastCost_beta = b;
    l.set_model(phylloptim::Leaf::CostCurve::LeastCost, false);
    l.optimise();
    return l;
  };
  // `r = xi/sqrt(D)`, the whole content of (FOC). Km and Gamma* are read off the
  // solved leaf rather than recomputed here, so a change to the Arrhenius curves
  // moves the prediction with the model instead of breaking the test.
  auto r_of = [&](const phylloptim::Leaf &l, double vpd_kpa, double b,
                  double ratio) {
    const double gamma_star = l.gamma_ * l.umol_per_mol_to_Pa_;         // Pa
    return std::sqrt(b * (l.km_ + gamma_star) / ratio) /
           std::sqrt(vpd_kpa * 1000.0);
  };

  // --- (FOC), which is the model, under both diffusion conventions -----------
  for (double ratio : {1.67, 1.6}) {
    for (double vpd : vpds) {
      phylloptim::Leaf l = solved(vpd, beta, ratio, true);
      ok(std::isfinite(l.profit_), "the curve solves at the published beta");
      ok(l.opt_psi_leaf_ > 0.1 + 1e-6 && l.opt_psi_leaf_ < l.psi_crit - 1e-6,
         "the optimum is interior, so the supply relation is not binding");
      const double gamma_star = l.gamma_ * l.umol_per_mol_to_Pa_;
      const double r = r_of(l, vpd, beta, ratio);
      const double chi = l.ci_ / l.ca_;
      const double foc = (gamma_star / l.ca_ + r) / (1.0 + r);
      printf("    ratio %.2f, D = %.1f kPa: chi = %.8f, (FOC) = %.8f, "
             "rel = %.1e\n", ratio, vpd, chi, foc, (chi - foc) / foc);
      ok(std::fabs(chi - foc) / foc < 1e-6,
         "the solved chi is Prentice et al.'s first-order condition to solver "
         "precision");
    }
  }

  // --- and the gap to the PUBLISHED form is the Gamma* it drops --------------
  //
  // ⚠️ THIS IS WHAT MAKES THE BLOCK ABOVE EVIDENCE RATHER THAN A DEFINITION. A
  // check against (FOC) alone could pass on a solver that had been fitted to
  // (FOC); predicting the published form's ERROR, in size and in how it moves
  // with deficit, cannot be arranged after the fact.
  for (double vpd : vpds) {
    phylloptim::Leaf l = solved(vpd, beta, 1.67, true);
    const double gamma_star = l.gamma_ * l.umol_per_mol_to_Pa_;
    const double r = r_of(l, vpd, beta, 1.67);
    const double chi = l.ci_ / l.ca_;
    const double published = r / (1.0 + r);
    const double predicted_gap = (gamma_star / l.ca_) / (1.0 + r);
    printf("    D = %.1f kPa: published = %.6f, gap = %.6f, "
           "Gamma*/ca/(1+r) = %.6f\n", vpd, published, chi - published,
           predicted_gap);
    ok(std::fabs((chi - published) - predicted_gap) < 1e-6,
       "the gap to the published form is exactly the Gamma* it neglects");
  }

  // --- that the solve READS beta, and reads it where (FOC) puts it -----------
  //
  // ⚠️ EVERY CHECK ABOVE RAN AT ONE BETA and would pass on a leaf that ignored
  // it. `xi` goes as sqrt(beta), and (FOC) says the quantity proportional to it
  // is `(ci - Gamma*)/(ca - ci)` -- NOT `chi/(1-chi)`, which is the published
  // form's invariant and is out by 4% here for the same reason the levels are.
  // This is what fails if `leastcost_flux` drops the capacity term, or weighs
  // `vcmax_25` where the condition wants `vcmax_`.
  {
    double q[2];
    const double betas[2] = {146.0, 292.0};
    for (int i = 0; i < 2; ++i) {
      phylloptim::Leaf l = solved(2.0, betas[i], 1.67, true);
      const double gamma_star = l.gamma_ * l.umol_per_mol_to_Pa_;
      q[i] = (l.ci_ - gamma_star) / (l.ca_ - l.ci_);
    }
    const double got = q[1] / q[0];
    printf("    doubling beta scales (ci-Gamma*)/(ca-ci) by %.8f "
           "(sqrt(2) = %.8f)\n", got, std::sqrt(2.0));
    ok(std::fabs(got - std::sqrt(2.0)) / std::sqrt(2.0) < 1e-6,
       "doubling beta scales (ci-Gamma*)/(ca-ci) by sqrt(2), so the solve reads "
       "the unit-cost ratio where the published condition puts it");
  }

  // --- what the agreement is worth at this package's OWN settings ------------
  //
  // Reported rather than asserted tightly, because it is a property of the
  // co-limitation this package models and Prentice et al. do not. The bound is
  // loose on purpose: it is here to catch a change of KIND, not to pin a number
  // that has no published counterpart.
  for (double vpd : vpds) {
    phylloptim::Leaf l = solved(vpd, beta, 1.67, false);
    const double gamma_star = l.gamma_ * l.umol_per_mol_to_Pa_;
    const double r = r_of(l, vpd, beta, 1.67);
    const double chi = l.ci_ / l.ca_;
    const double foc = (gamma_star / l.ca_ + r) / (1.0 + r);
    printf("    co-limited, D = %.1f kPa: chi = %.6f, (FOC) = %.6f, "
           "rel = %+.2f%%\n", vpd, chi, foc, 100.0 * (chi - foc) / foc);
    ok(std::fabs(chi - foc) / foc < 0.05,
       "at this package's co-limited default the solved chi is within 5% of "
       "(FOC), the residual being the co-limitation Prentice et al. omit");
  }
}

// ⚠️ THE GUARD IS STRICTER THAN `TF24_floor`'s, AND NOT OUT OF SYMMETRY.
// `TF24_floor` accepts an explicit zero because zero is meaningful there -- the
// curve is `TF24`. Zero here is not a limiting case of least-cost: it deletes
// the capacity term, leaving `A/E`, which is water-use efficiency and rises
// without bound as transpiration falls, so the optimiser returns the wet bound
// for every driver set. That looks like a solve and is not one.
void test_leastcost_refuses_an_unset_or_degenerate_beta() {
  printf("LeastCost refuses an unset, zero or negative beta, and says why\n");
  Drivers d;
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.LeastCost_beta = std::numeric_limits<double>::quiet_NaN();
    l.set_model(phylloptim::Leaf::CostCurve::LeastCost, false);
    ok(throws_with([&] { l.optimise(); }, "needs LeastCost_beta set"),
       "a cleared beta is refused rather than optimised as a NaN objective");
    // The message has to carry the units. A beta derived against this package's
    // kg H2O and umol CO2 rather than mol and mol is out by 5.6e7, and would
    // otherwise solve to something entirely plausible.
    ok(throws_with([&] { l.optimise(); }, "DIMENSIONLESS"),
       "and says what units the published number is in");
    ok(throws_with([&] { l.optimise(); }, "leaf_solve()"),
       "and names a route that can supply one");
  }
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.LeastCost_beta = 0.0;
    l.set_model(phylloptim::Leaf::CostCurve::LeastCost, false);
    ok(throws_with([&] { l.optimise(); }, "strictly positive"),
       "zero is refused: it deletes the capacity term rather than limiting it");
  }
  {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.LeastCost_beta = -1.0;
    l.set_model(phylloptim::Leaf::CostCurve::LeastCost, false);
    ok(throws_with([&] { l.optimise(); }, "strictly positive"),
       "a negative unit-cost ratio is refused");
  }
}

// The price of water a curve reports is the one quantity every member here is
// compared on, so it has to be the price this objective actually pays. Checked
// by finite difference against dA/dE at the solved optimum, which is the check
// `test_lambda_emergent_reported` applies to the others.
void test_leastcost_lambda_is_dA_dE() {
  printf("LeastCost's reported lambda is dA/dE at its own optimum\n");
  Drivers d;
  d.PPFD = 1500.0;
  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  l.LeastCost_beta = 146.0;
  l.set_model(phylloptim::Leaf::CostCurve::LeastCost, false);
  l.optimise();
  ok(std::isfinite(l.lambda_emergent_), "a lambda is reported");

  const double psi = l.opt_psi_leaf_;
  const double eps = 1e-6;
  auto A_of = [&](double p) {
    l.set_leaf_states_rates_from_psi_leaf(p, 0.5);
    return l.assim_colimited_;
  };
  auto E_of = [&](double p) { return l.transpiration(p, 0.5); };
  const double dA = (A_of(psi + eps) - A_of(psi - eps)) / (2.0 * eps);
  const double dE = (E_of(psi + eps) - E_of(psi - eps)) / (2.0 * eps);
  const double fd = dA / dE;
  printf("    lambda = %.6g, dA/dE = %.6g, rel = %+.3e\n",
         l.lambda_emergent_, fd, (l.lambda_emergent_ - fd) / fd);
  ok(std::fabs(l.lambda_emergent_ - fd) / std::fabs(fd) < 1e-4,
     "the reported price equals dA/dE by finite difference");

  // With the energy balance on, the capacity term moves with E through Tleaf,
  // and at an interior optimum the reported price is still dA/dE along the
  // stream. Beta 2 puts these optima in the interior (146 pins them).
  for (double tair : {25.0, 30.0}) {
    Drivers de;
    de.PPFD = 1500.0;
    de.leaf_temp = tair;
    phylloptim::Leaf m = make_linked_eb_leaf(de, 1.0, false);
    m.LeastCost_beta = 2.0;
    m.set_model(phylloptim::Leaf::CostCurve::LeastCost, false);
    m.optimise();
    const double p = m.opt_psi_leaf_;
    ok(p > 1.0 + 1e-6 && p < m.psi_crit - 1e-6,
       "the energy-balance optimum is interior (hazard 16)");
    const double lam = m.lambda_emergent_;
    auto Am = [&](double q) {
      m.set_leaf_states_rates_from_psi_leaf(q, 1.0);
      return m.assim_colimited_;
    };
    const double h = 1e-5;
    const double dAm = (Am(p + h) - Am(p - h)) / (2.0 * h);
    const double dEm =
        (m.transpiration(p + h, 1.0) - m.transpiration(p - h, 1.0)) / (2.0 * h);
    const double fde = dAm / dEm;
    printf("    energy balance on, Tair %.0f: lambda = %.6g, dA/dE = %.6g, "
           "rel = %+.3e\n",
           tair, lam, fde, (lam - fde) / fde);
    ok(std::fabs(lam - fde) / std::fabs(fde) < 1e-4,
       "with the energy balance on, the reported price equals dA/dE");
  }
}

void test_every_curve_returns_its_own_maximum() {
  printf("every cost curve's optimum beats a scan of its own objective\n");

  // ⚠️ THE HOT ROWS ARE NOT PADDING. A boundary optimum is where a scan oracle
  // earns its keep: an optimiser that steps in from the bounds reports an open
  // stoma where the objective says shut, and only a comparison of VALUES over a
  // closed interval catches it (hazard 11). Without rows hot enough to make
  // assimilation negative across the supply stream, this test would exercise the
  // interior branch alone -- which is the half that was never in doubt.
  struct Row { double psi_soil, PPFD, atm_vpd, leaf_temp; };
  const Row rows[] = {
      {0.5, 1500.0, 2.0, 25.0},
      {1.0, 900.0, 2.0, 25.0},
      {2.0, 900.0, 4.0, 30.0},
      {3.0, 1500.0, 2.0, 25.0},
      {1.0, 300.0, 1.0, 25.0},
      {0.5, 1500.0, 2.0, 50.0},   // shut
      {1.0, 1500.0, 3.0, 45.0},   // shut
      {2.0, 1200.0, 4.0, 48.0},   // shut
  };
  const int n_curves = phylloptim::Leaf::n_cost_curves;
  const int N = 20000;

  int total_interior = 0, total_pinned = 0, total_rows = 0;
  int min_rows_per_curve = std::numeric_limits<int>::max();
  double worst_shortfall = 0.0;
  std::string worst_where;
  int dark_shut = 0;
  bool dark_at_wet = true;
  for (int curve = 0; curve < n_curves; ++curve) {
    double curve_worst = 0.0;
    int interior = 0, pinned = 0, curve_rows = 0;
    for (const Row &r : rows) {
      Drivers d;
      d.PPFD = r.PPFD;
      d.atm_vpd = r.atm_vpd;
      d.leaf_temp = r.leaf_temp;

      phylloptim::Leaf l = make_single_leaf(d, r.psi_soil);
      // The three curves with no default: prescribed, not derived, and unset is
      // NA. Set unconditionally, because a field a curve does not read costs
      // nothing -- and because an unset one does not degrade gracefully. This
      // loop runs to `n_cost_curves`, so a curve added to the enum arrives here
      // whether or not anyone remembers to; `LeastCost` did, and aborted the
      // whole binary on its refusal until its beta was set here.
      l.CF77_lambda_ = 1.5e5;
      l.TF24_floor_lambda_o = kFloorLambdaO;
      l.LeastCost_beta = 146.0;
      l.set_model(static_cast<phylloptim::Leaf::CostCurve>(curve), false);
      l.optimise();
      const double got = l.profit_;
      if (!std::isfinite(got)) {
        continue;   // a refused row; psi_leaf_optima.tsv covers the refusals
      }
      // Fresh leaf, so the solve's own state cannot influence the oracle.
      phylloptim::Leaf m = make_single_leaf(d, r.psi_soil);

      // A log-link objective with A <= 0 throughout has no maximum (a
      // product's supremum is at psi_crit, a ratio's at the most open leaf),
      // and the leaf is shut by rule rather than by its objective
      // (objective_has_log_link): a scan beats it there by construction. What it must do instead is sit at the wet end. Decided
      // on the oracle, at psi_crit, not on the solve: a shut leaf at g_min = 0
      // has A = 0 exactly, so testing the solve's A would excuse a light leaf
      // the solve shut by mistake.
      const bool log_link =
          curve == static_cast<int>(phylloptim::Leaf::CostCurve::SOX) ||
          curve == static_cast<int>(phylloptim::Leaf::CostCurve::JW26) ||
          curve == static_cast<int>(phylloptim::Leaf::CostCurve::LeastCost);
      if (log_link) {
        m.set_leaf_states_rates_from_psi_leaf(m.psi_crit, r.psi_soil);
        if (!(m.assim_colimited_ > 0.0)) {
          ++dark_shut;
          dark_at_wet = dark_at_wet && l.opt_psi_leaf_ == l.leaf_wet_bound(r.psi_soil);
          continue;
        }
      }
      m.CF77_lambda_ = 1.5e5;
      m.TF24_floor_lambda_o = kFloorLambdaO;
      m.LeastCost_beta = 146.0;
      // ⚠️ ProfitMax's objective is UNDEFINED on a leaf that has not scanned for
      // |A|max, so a fresh oracle leaf returns NaN at every point and the curve
      // silently contributes nothing -- which is how this test first passed while
      // testing six curves and reporting seven. Seed it with the SAME normaliser
      // the solver used: a different |A|max is a different objective, so comparing
      // against one would not be comparing against this curve at all.
      if (curve == static_cast<int>(phylloptim::Leaf::CostCurve::ProfitMax)) {
        m.prepare_profitmax_at(l.profitmax_A_max_);
      }
      m.set_model(static_cast<phylloptim::Leaf::CostCurve>(curve), false);
      const double lo = r.psi_soil, hi = m.psi_crit;
      double best = -std::numeric_limits<double>::infinity(), best_psi = lo;
      for (int i = 0; i <= N; ++i) {
        const double q = lo + (hi - lo) * (double(i) / double(N));
        const double v = m.evaluate_psi_leaf_at(q);
        if (std::isfinite(v) && v > best) { best = v; best_psi = q; }
      }
      if (!std::isfinite(best)) {
        continue;
      }
      ++total_rows;
      ++curve_rows;
      // Relative slack, because the eight objectives differ in units and by orders
      // of magnitude -- ProfitMax is dimensionless and O(1), TF24 is carbon.
      const double shortfall =
          (best - got) / std::max(std::fabs(best), 1.0);
      if (shortfall > curve_worst) {
        curve_worst = shortfall;
        if (shortfall > worst_shortfall) {
          worst_shortfall = shortfall;
          worst_where = phylloptim::Leaf::curve_name(curve) + " psi_soil=" +
                        std::to_string(r.psi_soil);
        }
      }
      const double edge = (hi - lo) * 1e-3;
      if (best_psi <= lo + edge || best_psi >= hi - edge) { ++pinned; } else { ++interior; }
    }
    printf("    %-10s %d rows | worst shortfall %.3e | %d interior, %d pinned\n",
           phylloptim::Leaf::curve_name(curve).c_str(), curve_rows, curve_worst,
           interior, pinned);
    if (curve_rows < min_rows_per_curve) { min_rows_per_curve = curve_rows; }
    total_interior += interior;
    total_pinned += pinned;
  }

  // The premise, asserted so the block cannot pass by covering one regime or by
  // quietly evaluating nothing (hazard 16).
  printf("    %d rows | %d interior, %d pinned | worst %.3e at %s\n", total_rows,
         total_interior, total_pinned, worst_shortfall, worst_where.c_str());
  // ⚠️ PER CURVE, not in total. A total is satisfied by the other six while one
  // contributes nothing, which is exactly what happened here before the ProfitMax
  // seeding above: 30 rows from six curves passed a `>= 30` check.
  ok(min_rows_per_curve > 0, "EVERY curve contributed solved rows, not just six");
  ok(total_interior > 0 && total_pinned > 0,
     "and the rows span both the interior and the pinned regime");
  ok(worst_shortfall < 1e-6,
     "no curve's solver is beaten by a 20001-point scan of its own objective");
  printf("    %d SOX/JW26/LeastCost rows with A <= 0 throughout, shut by rule\n", dark_shut);
  ok(dark_shut > 0 && dark_at_wet,
     "a log-link objective's dark rows sit at the wet end, not at psi_crit");
}

// ProfitMax with no objective (in the dark, or with the soil past psi_crit so
// there is no stream to normalise over) takes the leaf route's degenerate exit
// on the collar route too: zero outputs at the wet bound, classified, never a
// thrown solve or a profit from unusable normalisers. A replay hands the same
// zeros back, and the surface refuses the point by name.
void test_profitmax_collar_route_has_a_degenerate_exit() {
  printf("ProfitMax's collar route takes the degenerate exit\n");
  namespace pl = phylloptim;
  using CC = pl::Leaf::CostCurve;
  using Kind = pl::Leaf::OperatingPointKind;
  int rows = 0, shut = 0;
  bool zeros = true, classified = true, replayed = true, refused = true, open = true;
  bool solved = true;
  for (int layers : {1, 3}) for (double psi : {0.5, 6.0, 17.0})
  for (double ppfd : {0.0, 900.0}) {
    Drivers d;
    d.PPFD = ppfd;
    std::vector<double> ps, dp;
    for (int i = 0; i < layers; ++i) { ps.push_back(psi + 0.25 * i); dp.push_back(i + 1.0); }
    pl::Leaf l = make_leaf(d, ps, dp);
    l.set_model(CC::ProfitMax, true);
    try {
      l.optimise();
    } catch (const std::exception&) {
      solved = false;
      continue;
    }
    const std::string at = " (" + std::to_string(layers) + " layers, psi_soil " +
                           std::to_string(psi) + ", PPFD " + std::to_string(int(ppfd)) + ")";
    if (ppfd > 0.0 && psi == 0.5) {
      open = open && l.profitmax_has_objective() &&
             l.operating_point_kind() == Kind::Interior && l.profit_ > 0.0;
      continue;
    }
    ++rows;
    const Kind kind = l.operating_point_kind();
    if (kind == Kind::HydraulicShutdown) ++shut;
    classified = classified && !l.profitmax_has_objective() &&
                 (kind == Kind::BoundarySoil || kind == Kind::HydraulicShutdown);
    zeros = zeros && l.profit_ == 0.0 && l.carbon_gain_ == 0.0 &&
            l.hydraulic_cost_norm_ == 0.0 && std::isnan(l.lambda_emergent_);
    {
      pl::Leaf r = make_leaf(d, ps, dp);
      r.set_model(CC::ProfitMax, true);
      try {
        replayed = replayed && r.replay_operating_point(l.opt_root_psi_, kind) == 0.0 &&
                   r.profit_ == 0.0 && r.opt_psi_leaf_ == l.opt_psi_leaf_ &&
                   r.operating_point_kind() == kind;
      } catch (const std::exception&) {
        replayed = false;
      }
    }
    if (kind == Kind::BoundarySoil) {
      const auto supply = l.roots_.held_supply();
      const auto draw = l.supply_draw_at<double>(l.opt_root_psi_, supply);
      bool threw = false;
      try {
        l.normalisers_at<CC::ProfitMax, double>(draw, l.passive_pars());
      } catch (const std::exception& e) {
        threw = std::string(e.what()).find("no objective") != std::string::npos;
      }
      refused = refused && threw;
      if (!threw) printf("    not refused%s\n", at.c_str());
    }
  }
  printf("    %d rows without an objective, %d of them shut on water\n", rows, shut);
  ok(solved, "no solve throws");
  ok(rows == 10 && shut > 0 && classified,
     "each is classified, at the wet bound or shut on water");
  ok(zeros, "with the leaf route's zero outputs and lambda NA");
  ok(replayed, "a replay hands back the same zeros");
  ok(refused, "the surface refuses it by name");
  ok(open, "and a lit leaf in moist soil still has its objective");
}

// ===========================================================================
// The residual (cuticular) conductance, #182
// ---------------------------------------------------------------------------
// ProfitMax's |A|max with the energy balance off is |A| at psi_crit, read in one
// evaluation: at a fixed leaf temperature |A| rises with conductance at either
// sign of A. Held against a 400-point walk of the stream, with and without a
// cuticle and the damage layer, on dark and bright, cool and hot leaves, so a
// normaliser read anywhere but the peak fails here.
void test_profitmax_abs_assim_max_is_at_psi_crit() {
  printf("ProfitMax: |A|max is |A| at psi_crit with the energy balance off\n");
  int cases = 0, negative = 0;
  double worst_excess = 0.0;
  bool exact_at_crit = true;
  for (double g_min : {0.0, 5e-4, 5e-3})
    for (bool damage : {false, true})
      for (double psi_soil : {0.5, 2.0, 4.0})
        for (double ppfd : {1.0, 200.0, 1500.0})
          for (double temp : {25.0, 50.0}) {
            Drivers d;
            d.PPFD = ppfd;
            d.leaf_temp = temp;
            phylloptim::Leaf l = make_single_leaf(d, psi_soil);
            l.g_min = g_min;
            if (damage) {
              l.use_thermal_damage_ = true;
              set_atls_inputs(l);
            }
            // Re-seated, since set_physiology reads g_min and the damage gate.
            l = [&] {
              phylloptim::Leaf m = l;
              phylloptim::RootNetwork rn = fixture::series_resistance(1.0e3);
              m.set_physiology(rn, d.PPFD, {psi_soil}, {1.0}, d.K_s * d.theta / d.h,
                               d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa,
                               d.atm_kpa);
              return m;
            }();
            l.prepare_profitmax();
            const double a_max = l.profitmax_A_max_;
            const double ref = l.profitmax_reference_psi();
            const double wet = l.leaf_wet_bound(ref);
            if (!(a_max > 0.0) || !(l.psi_crit > wet)) {
              continue;
            }
            ++cases;
            const int n = 400;
            for (int i = 1; i <= n; ++i) {
              const double psi = wet + (l.psi_crit - wet) * double(i) / n;
              l.set_leaf_states_rates_from_psi_leaf(psi, ref);
              worst_excess = std::max(
                  worst_excess, (std::abs(l.assim_colimited_) - a_max) / a_max);
            }
            // The walk's last point is psi_crit itself.
            if (std::abs(l.assim_colimited_) != a_max) exact_at_crit = false;
            if (l.assim_colimited_ < 0.0) ++negative;
          }
  printf("    %d cases (%d with A < 0), worst |A| above |A|max %.3g of it\n", cases,
         negative, worst_excess);
  ok(cases == 108 && negative >= 20, "the walk reaches leaves of either sign of A");
  ok(worst_excess <= 0.0, "no point of the stream has a larger |A| than |A|max");
  ok(exact_at_crit, "|A|max is |A| at psi_crit, bit for bit");

  // A dark leaf with no respiration has A = 0 along the whole stream: the one
  // evaluation finds no |A|max, and the search it falls back to finds none
  // either, so the leaf has no objective rather than a stale normaliser.
  Drivers dark;
  dark.PPFD = 0.0;
  phylloptim::Leaf z = make_single_leaf(dark, 1.0);
  z.profitmax_A_max_ = 1.0;  // poisoned: the preparation must write it
  z.R_d_25 = 0.0;
  z.set_physiology(fixture::series_resistance(1.0e3), dark.PPFD, {1.0}, {1.0},
                   dark.K_s * dark.theta / dark.h, dark.atm_vpd, dark.ca,
                   dark.leaf_temp, dark.atm_o2_kpa, dark.atm_kpa);
  z.prepare_profitmax();
  ok(z.profitmax_A_max_ == 0.0 && !z.profitmax_has_objective(),
     "a leaf with A = 0 throughout has no |A|max and no objective");
}

// SOX, JW26 and LeastCost on the collar route, where A <= 0 at every opening:
// the leaf is shut at the wet bound (objective_has_log_link) rather than
// reported solver-refused or opened wide, on one and three layers, with and
// without a residual conductance or the energy balance. The same leaf in light
// is open (LeastCost mostly at psi_crit), so a shut that fired on light would
// show here.
void test_log_link_curves_shut_dark_collar_leaves() {
  printf("SOX, JW26 and LeastCost shut a collar-route leaf that cannot fix carbon\n");
  using CC = phylloptim::Leaf::CostCurve;
  using Kind = phylloptim::Leaf::OperatingPointKind;
  int dark = 0, light = 0;
  bool dark_shut = true, dark_finite = true, light_open = true;
  for (CC curve : {CC::SOX, CC::JW26, CC::LeastCost})
  for (int layers : {1, 3}) for (double g : {0.0, residual::kGmin}) for (bool eb : {false, true})
  for (double psi : {0.5, 1.5}) for (double ppfd : {0.0, 900.0}) {
    Drivers d;
    d.PPFD = ppfd;
    d.leaf_temp = 25.0;
    std::vector<double> ps, dp;
    for (int i = 0; i < layers; ++i) { ps.push_back(psi + 0.25 * i); dp.push_back(i + 1.0); }
    phylloptim::Leaf l = residual::make(d, ps, dp, g, eb);
    l.set_model(curve, true);
    l.optimise();
    if (ppfd == 0.0) {
      ++dark;
      dark_shut = dark_shut && l.operating_point_kind() == Kind::BoundarySoil &&
                  !(l.assim_colimited_ > 0.0);
      // A Log link has no slope at A <= 0, so lambda is NA, written not stale.
      dark_finite = dark_finite && std::isfinite(l.profit_) &&
                    std::isnan(l.lambda_emergent_);
    } else {
      ++light;
      const Kind open = curve == CC::LeastCost ? Kind::BoundaryDry : Kind::Interior;
      light_open = light_open &&
                   (l.operating_point_kind() == Kind::Interior ||
                    l.operating_point_kind() == open) &&
                   l.assim_colimited_ > 0.0;
    }
  }
  printf("    %d dark rows, %d light rows\n", dark, light);
  ok(dark == 48 && dark_shut, "every dark row is shut at the wet bound (boundary-soil)");
  ok(dark_finite, "with finite profit, and lambda NA");
  ok(light == 48 && light_open, "and the same leaves in light are open");
}

}  // namespace

void run_curves() {
  test_collar_lambda_emergent_is_dA_dE();
  test_cowan_farquhar_equates_dA_dE_to_lambda();
  test_cowan_farquhar_reproduces_the_TF_optimum();
  test_cowan_farquhar_closed_state();
  test_cowan_farquhar_refuses_an_unset_lambda();
  test_aligned_defaults_match_tf24_at_the_reference();
  test_closed_form();
  test_profitmax_reports_an_emergent_lambda();
  test_every_curve_reports_an_emergent_lambda();
  test_profitmax_emergent_lambda_matches_a_finite_difference();
  test_profitmax_writes_hydraulic_cost_on_every_exit();
  test_profitmax_is_a_constant_lambda_objective();
  test_profitmax_normalisation();
  test_profitmax_thermal_cost();
  test_collar_foc_matches_profit_linked_curves();
  test_closed_form_is_prescribed_temperature_only();
  test_profitmax_finds_a_closed_optimum();
  test_product_link_is_the_product_rule();
  test_shadow_cost();
  test_lambda_emergent_reported();
  test_tf24_floor_reduces();
  test_tf24_floor_wet_end_price();
  test_tf24_floor_emergent_lambda();
  test_tf24_floor_kmax_vpd_invariance();
  test_tf24_floor_defaults_to_tf24();
  test_tf24_floor_closed_form();
  test_leastcost_matches_prentice_closed_form();
  test_leastcost_refuses_an_unset_or_degenerate_beta();
  test_leastcost_lambda_is_dA_dE();
  test_every_curve_returns_its_own_maximum();
  test_log_link_curves_shut_dark_collar_leaves();
  test_profitmax_collar_route_has_a_degenerate_exit();
  test_profitmax_abs_assim_max_is_at_psi_crit();
}
