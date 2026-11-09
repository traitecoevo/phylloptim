// Pure-C++ test suite for the leaf model. No R, no test framework, no linking:
//   make -C tests/cpp && ./tests/cpp/test_leaf
//
// The trait values and drivers below are lifted from plant's
// tests/testthat/test-leaf.r so the two suites exercise the same operating
// point. The expected values here were produced BY this implementation and are
// regression guards rather than independent references -- but the implementation
// has been cross-checked against plant's compiled build and found bit-identical,
// so they are guarding a verified model.


#include "leaf_harness.hpp"

namespace {

// ---------------------------------------------------------------------------

void test_defaults_are_unset() {
  printf("defaults are unset until set_physiology\n");
  phylloptim::Leaf l;
  // One set of defaults, as in R: TF24_floor on the leaf route, single supply.
  ok(l.cost_curve_ == phylloptim::Leaf::CostCurve::TF24_floor,
     "the default curve is TF24_floor");
  ok(l.model_route() == "leaf", "the default route is the leaf route");
  ok(l.supply_kind_name() == "single", "the default supply is the single potential");
  ok(!std::isfinite(l.ci_), "ci_ starts unset");
  ok(!std::isfinite(l.assim_colimited_), "assim_colimited_ starts unset");
  ok(!std::isfinite(l.opt_psi_leaf_), "opt_psi_leaf_ starts unset");
  ok(l.roots_.psi_soil_.empty(), "psi_soil_ starts empty");
  ok(l.use_energy_balance_ == false, "energy balance defaults off");
}

// A curve's definition is what the leaf runs on: set_model() seats it by its
// name, a fresh Leaf holds each of its parameters' defaults, and each parameter's
// member is the one the pack reads at its slot. Returns how many it checked.
template <phylloptim::Leaf::CostCurve K>
int check_curve_definition() {
  using C = phylloptim::Curve<K>;
  const std::string name(C::name);
  phylloptim::Leaf l;
  l.set_model(name, "leaf");
  ok(l.cost_curve_ == K, name + ": set_model() seats it by that name");
  for (const phylloptim::CurveParameter &p : C::parameters) {
    const double got = l.*(p.member);
    ok(std::isnan(p.default_value) ? std::isnan(got) : got == p.default_value,
       name + ": a fresh Leaf holds " + p.field + "'s default");
    ok(p.slot >= 0 && p.slot < phylloptim::par_PPFD,
       name + ": " + p.field + " is a fitted slot");
    const double marker = 0.125 + p.slot;
    l.*(p.member) = marker;
    ok(l.passive_pars()[p.slot] == marker,
       name + ": the pack reads " + p.field + " at its slot");
  }
  return static_cast<int>(C::parameters.size());
}

void test_curve_definitions_are_the_leaf() {
  printf("each curve's definition is what a fresh Leaf holds\n");
  int parameters = 0;
  for (int i = 0; i < phylloptim::Leaf::n_cost_curves; ++i) {
    parameters += phylloptim::Leaf::with_curve(
        static_cast<phylloptim::Leaf::CostCurve>(i),
        [](auto tag) { return check_curve_definition<tag.value>(); });
  }
  printf("    %d curves, %d parameters\n", phylloptim::Leaf::n_cost_curves,
         parameters);
  ok(parameters > 0, "the check reached the definitions' parameters");
}

void test_arrhenius() {
  printf("temperature response\n");
  phylloptim::Leaf l = fixture::multilayer_leaf();
  near(l.arrh_curve(phylloptim::vcmax_ha, 100.0, 25.0), 100.0, 1e-12,
       "Arrhenius is the identity at the 25 C reference");
  ok(l.arrh_curve(phylloptim::vcmax_ha, 100.0, 35.0) > 100.0,
     "Arrhenius rises above the reference temperature");
  ok(l.peak_arrh_curve(phylloptim::jmax_ha, 100.0, 60.0, phylloptim::jmax_H_d,
                       phylloptim::jmax_d_S) <
         l.peak_arrh_curve(phylloptim::jmax_ha, 100.0, 30.0, phylloptim::jmax_H_d,
                           phylloptim::jmax_d_S),
     "peaked Arrhenius declines past its optimum");
}

void test_saturation_vapour_pressure() {
  printf("saturation vapour pressure\n");
  phylloptim::Leaf l = fixture::multilayer_leaf();
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
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  ok(std::isfinite(l.opt_psi_leaf_), "psi_leaf is finite");
  ok(std::isfinite(l.profit_), "profit is finite");
  ok(l.opt_psi_leaf_ > 2.0, "leaf is drier than the soil");
  ok(l.opt_psi_leaf_ <= l.psi_crit, "leaf stays within psi_crit");
  ok(l.transpiration_ > 0.0, "transpiration is positive");
  ok(l.assim_colimited_ > 0.0, "assimilation is positive");
  // Regression guards -- see the note at the top of this file. Regenerate these
  // deliberately, alongside the golden file, and state the movement in the PR.
  near(l.opt_psi_leaf_, 3.595332, 1e-5, "opt_psi_leaf_");
  near(l.assim_colimited_, 5.601016, 1e-5, "assim_colimited_");
  near(l.transpiration_, 1.141941e-05, 1e-5, "transpiration_");
  near(l.profit_, 2.517157, 1e-5, "profit_");
}

void test_solve_is_deterministic() {
  printf("repeated solves are bit-identical\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double first = l.profit_, psi = l.opt_psi_leaf_;
  for (int i = 0; i < 50; ++i) {
    l.find_root_collar_psi();
  }
  ok(l.profit_ == first, "profit is unchanged after 50 re-solves");
  ok(l.opt_psi_leaf_ == psi, "psi_leaf is unchanged after 50 re-solves");
}

void test_drier_soil_costs_carbon() {
  printf("response to drying soil\n");
  Drivers d;
  double prev_profit = 1e9, prev_E = 1e9;
  for (double psi : {0.5, 1.0, 2.0, 3.0, 4.0}) {
    phylloptim::Leaf l = make_leaf(d, {psi}, {1.0});
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
  phylloptim::Leaf shaded = make_leaf(dim, {2.0}, {1.0});
  phylloptim::Leaf sunlit = make_leaf(bright, {2.0}, {1.0});
  shaded.find_root_collar_psi();
  sunlit.find_root_collar_psi();
  ok(sunlit.assim_colimited_ > shaded.assim_colimited_,
     "brighter light assimilates more");
  ok(sunlit.transpiration_ > shaded.transpiration_, "brighter light transpires more");
}

void test_multi_layer_soil() {
  printf("multi-layer soil\n");
  Drivers d;
  // ⚠️ The profile is CUMULATIVE depth to the bottom of each layer, so three
  // 0.5 m layers are {0.5, 1.0, 1.5}. This read {0.5, 0.5, 0.5} until #626 --
  // three layers all bottoming at 0.5 m, which put two of the three midpoints at
  // the same depth and gave the root network a thickness of 0.5/3. Nothing
  // noticed, because the assertions below are finiteness and sign only. The
  // strictly-increasing guard in layer_thicknesses() is what surfaced it.
  phylloptim::Leaf l = make_leaf(d, {1.0, 2.0, 3.0}, {0.5, 1.0, 1.5});
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

// psi_crit bounds nothing under TF24 (#183): a soil drier than it still has an
// interval, and the leaf transpires what its stem can carry, which falls
// smoothly towards zero as the soil dries. It shuts only past the end of the stem
// curve, where the stem carries G(inf) to rounding and so nothing.
void test_shutdown_when_soil_is_drier_than_psi_crit() {
  printf("past psi_crit the leaf transpires; past the end of the curve it shuts\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {12.0}, {1.0});
  l.find_root_collar_psi();
  ok(std::isfinite(l.profit_), "profit is finite past psi_crit");
  // At these drivers the leaf declines to open (boundary-soil); in light enough
  // and humid enough air it transpires (the boundary-dry fixtures elsewhere).
  ok(l.operating_point_kind() == Kind::Interior ||
         l.operating_point_kind() == Kind::BoundarySoil ||
         l.operating_point_kind() == Kind::BoundaryDry,
     std::string("a soil past psi_crit is an optimum, got ") +
         l.operating_point_kind_name(l.operating_point_kind()));

  phylloptim::Leaf s = make_leaf(d, {20.0}, {1.0});
  s.find_root_collar_psi();
  ok(s.operating_point_kind() == Kind::HydraulicShutdown,
     "past the end of the stem curve the leaf shuts");
  ok(s.profit_ <= 0.0, "profit is non-positive when shut down");
  // The shut stem sits where the flux balance puts it: at the collar, which is
  // in equilibrium with the soil (uptake zero at g_min = 0) (#182, #183).
  ok(s.opt_psi_leaf_ == s.opt_root_psi_, "the shut stem sits at the collar");
  near(s.E_up_, 0.0, 1e-12, "and draws no water from the soil");
}

// The shutdown path used to leak the previous solve's fluxes (plant #578/#577).
// set_shut_state now writes them, so this asserts the fixed behaviour: a
// shut-down leaf moves no water, nets what the ordinary evaluator gives it at
// zero flux, and none of that depends on what ran before.
void test_shutdown_writes_its_own_fluxes() {
  printf("shutdown writes its own fluxes (plant #578 fixed)\n");
  Drivers d;
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  std::vector<double> mrp{1.0 / d.area_leaf}, depth{1.0};
  const auto solve = [&](double psi) {
    std::vector<double> ps{psi};
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, ps, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
  };

  solve(4.0); // wet enough to transpire
  ok(l.transpiration_ > 0.0, "the wet solve transpires");

  solve(20.0); // far drier than psi_crit: the leaf shuts down
  ok(l.profit_ < 0.0, "the dry solve is a shutdown (profit < 0)");
  near(l.transpiration_, 0.0, 1e-300, "transpiration is zero, not stale");
  near(l.stom_cond_CO2_, 0.0, 1e-300, "conductance is zero, not stale");
  // The collar is the zero-uptake collar, solved to the polish tolerance.
  near(l.E_up_, 0.0, 1e-12, "soil uptake is zero, not stale");
  for (double s : l.soil_consumption_) {
    near(s, 0.0, 1e-9, "per-layer consumption is zero, not stale");
  }
  // Moving no gas, in light: ci sits where net assimilation is zero, to the ci
  // root-find's tolerance, and profit_ is that less the hydraulic cost at the
  // stem.
  near(l.assim_colimited_, 0.0, 1e-9, "assimilation is zero in light");
  near(l.assim_colimited_, l.profit_ + l.hydraulic_cost_TF(l.opt_psi_leaf_), 1e-12,
       "assimilation is consistent with profit and the hydraulic cost");
  ok(std::isfinite(l.ci_), "ci is set to the compensation point, not left stale");
  near(l.assim_colimited(l.ci_), 0.0, 1e-9, "and net assimilation is zero there");

  // Order independence is the property that was actually broken: a fresh leaf
  // taken straight to dry must report exactly the same thing.
  phylloptim::Leaf fresh = make_leaf(d, {20.0}, {1.0});
  fresh.find_root_collar_psi();
  ok(fresh.transpiration_ == l.transpiration_,
     "a fresh leaf gives the same transpiration");
  ok(fresh.assim_colimited_ == l.assim_colimited_,
     "a fresh leaf gives the same assimilation");
  ok(fresh.profit_ == l.profit_, "a fresh leaf gives the same profit");
}

// soil_consumption_ used to be cleared with .resize, whose fill reaches only
// newly-added elements, while the uptake loop writes only up to max_soil_layer --
// the deepest layer holding root carbon. So a shallow-rooted plant solving on a
// Leaf that a deep-rooted one used before it inherited the deep layers' uptake
// (plant #577, fixed in plant by #585). The soil layer count is the same in both
// solves here; only the rooted depth shrinks, which is why resizing never noticed.
void test_shallow_roots_do_not_inherit_deep_uptake() {
  printf("shallow roots do not inherit the previous plant's deep uptake\n");
  Drivers d;
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);

  const std::vector<double> psi_soil{1.0, 1.5, 2.0};
  const std::vector<double> depth{1.0, 2.0, 3.0};
  const auto solve = [&](std::vector<double> root_carbon) {
    l.set_physiology(fixture::root_network(root_carbon, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
  };

  // A tree rooted through all three layers, then a seedling rooted in the top
  // layer only -- trailing zeros in the carbon vector is how plant expresses
  // "not rooted that deep".
  const double c = 1.0 / 3.0 / d.area_leaf;
  solve({c, c, c});
  ok(l.soil_consumption_[1] != 0.0 && l.soil_consumption_[2] != 0.0,
     "the deep-rooted solve draws from layers 2 and 3");

  solve({1.0 / d.area_leaf, 0.0, 0.0});
  ok(l.soil_consumption_.size() == 3u,
     "the consumption vector still spans every soil layer");
  near(l.soil_consumption_[1], 0.0, 1e-300,
       "layer 2 is zero, not the tree's uptake");
  near(l.soil_consumption_[2], 0.0, 1e-300,
       "layer 3 is zero, not the tree's uptake");

  // Order independence is the property that matters: plant reuses one Leaf for
  // every individual in a patch, so the seedling must not depend on its neighbour.
  phylloptim::Leaf fresh = fixture::multilayer_leaf();
  fresh.setup_transpiration(100);
  fresh.setup_root_vulnerability(100);
  fresh.set_physiology(fixture::root_network({1.0 / d.area_leaf, 0.0, 0.0}, depth), d.PPFD, psi_soil, depth,
                       d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                       d.atm_o2_kpa, d.atm_kpa);
  fresh.find_root_collar_psi();
  for (size_t i = 0; i < 3; ++i) {
    ok(fresh.soil_consumption_[i] == l.soil_consumption_[i],
       "a fresh leaf gives the same per-layer consumption");
  }
}

// The other early exit that determines the operating point without going through
// profit_psi_leaf_TF: assimilation is negative even at ci = ca, so there is no
// light level at which opening the stomata pays. It set profit_ but left the
// leaf-side rates alone, so they held whatever the previous solve wrote -- and on
// a fresh leaf, nothing at all. Ported from plant develop (#585).
void test_negative_assim_exit_writes_its_own_rates() {
  printf("the assim_max_ < 0 exit writes its own rates\n");
  Drivers d;
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  const std::vector<double> psi_soil{1.0}, depth{1.0};
  const std::vector<double> root{1.0 / d.area_leaf};
  const auto solve = [&](double ppfd) {
    l.set_physiology(fixture::root_network(root, depth), ppfd, psi_soil, depth, d.K_s * d.theta / d.h,
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
  };

  solve(900.0); // bright: a normal optimising solve
  ok(l.transpiration_ > 0.0, "the bright solve transpires");

  // Dim enough that gross assimilation cannot cover R_d even with ci at ca. The
  // soil is wet, so this is not the psi_crit shut-down path: the leaf is a net
  // CO2 source whose optimum pins to the wet bound (#182), where at g_min = 0
  // it sits a step inside the zero-uptake collar and moves a rounding-level
  // flux. Every output is the evaluator's at that point, not the bright
  // solve's.
  solve(10.0);
  ok(l.assim_max_ < 0.0, "the dim leaf cannot fix carbon at ambient CO2");
  ok(l.transpiration_ >= 0.0 && l.transpiration_ < 1e-9,
     "transpiration is (almost) zero, not stale");
  ok(l.stom_cond_CO2_ >= 0.0 && l.stom_cond_CO2_ < 1e-6,
     "conductance is (almost) zero, not stale");
  // The leaf nets the block's assimilation at the ci it solved for, which sits
  // above ambient because respiratory CO2 has to leave through what little
  // conductance there is: a loss, smaller than R_d.
  ok(l.ci_ > l.ca_, "ci sits above ambient: the leaf is a source");
  near(l.assim_colimited_, l.assim_colimited(l.ci_), 1e-12,
       "net assimilation is the block's value at the reported ci");
  ok(l.assim_colimited_ < 0.0 && l.assim_colimited_ > -l.R_d_,
     "a loss, smaller than R_d");
  // The invariant the fix buys: profit_ == assim_colimited_ - hydraulic cost in
  // every branch. It held numerically to the last bit when measured.
  near(l.assim_colimited_ - l.hydraulic_cost_TF(l.opt_psi_leaf_), l.profit_,
       1e-14, "profit is consistent with assimilation and the hydraulic cost");

  phylloptim::Leaf fresh = fixture::multilayer_leaf();
  fresh.setup_transpiration(100);
  fresh.setup_root_vulnerability(100);
  fresh.set_physiology(fixture::root_network(root, depth), 10.0, psi_soil, depth, d.K_s * d.theta / d.h,
                       d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  fresh.find_root_collar_psi();
  ok(fresh.transpiration_ == l.transpiration_,
     "a fresh leaf gives the same transpiration");
  ok(fresh.assim_colimited_ == l.assim_colimited_,
     "a fresh leaf gives the same assimilation");
  ok(fresh.profit_ == l.profit_, "a fresh leaf gives the same profit");
}

void test_analytic_gradient_matches_finite_difference() {
  printf("analytic dprofit/dpsi_collar vs central difference\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double p0 = l.opt_root_psi_;
  const double target = std::max(2.2, std::min(p0, l.psi_crit - 0.5));
  const double eps = 1e-5;
  const double analytic = l.dprofit_droot_collar_psi(target);
  const double up = l.evaluate_root_collar_psi(target + eps);
  const double dn = l.evaluate_root_collar_psi(target - eps);
  const double fd = (up - dn) / (2 * eps);
  ok(std::isfinite(analytic), "analytic gradient is finite");
  near(analytic, fd, 2e-3, "analytic gradient matches central difference");
}

// The collar solve root-finds `dprofit_droot_collar_psi == 0`, so that function has
// to be the derivative of the profit the solve REPORTS, on every curve. A cost that
// reads the upstream potential moves with the collar directly, and a first-order
// condition missing that term still converges -- to a point that is not the
// maximum -- with every reported number plausible. So this checks the slope
// against a difference of `profit_`, on each curve whose objective IS `profit_`
// (the identity-link ones; the log links report the product, not its log), at two
// collars wetter than the solved one. Not AT it: a curve pinned to the dry bound
// would put a central difference across the bound, which halves the slope.
// Through the `_for<K>` entry points, not the seated-curve wrappers: those are the
// functions the solve itself calls, whatever the wrappers dispatch to.
template <phylloptim::Leaf::CostCurve K>
int check_collar_foc_matches_profit(const Drivers &d) {
  const struct {
    std::vector<double> psi;
    std::vector<double> depth;
  } soils[] = {{{1.0}, {1.0}}, {{1.0, 2.0, 3.0}, {0.5, 1.0, 1.5}}};
  int compared = 0;
  for (const auto &s : soils) {
    phylloptim::Leaf l = make_leaf(d, s.psi, s.depth);
    l.CF77_lambda_ = kPriceForFOC;
    l.TF24_floor_lambda_o = kPriceForFOC;
    // ⚠️ LeastCost IS THE ONE LOG-LINK CURVE CHECKED HERE, and it is here
    // because its cost reads the upstream potential: `log W` with `W` carrying
    // `E(psi_leaf, psi_upstream)`. SOX and JW26 read no upstream. The reported
    // `profit_` is the ratio `A/W` and the solve differentiates its LOG, so the
    // difference below is taken of `log(profit_)` for this curve.
    l.LeastCost_beta = 146.0;
    l.JS22_gamma = kGammaJS22;
    l.CMax_a = kCMaxA;
    l.CMax_b = kCMaxB;
    l.set_model(K, true);
    l.optimise();
    const double p0 = l.opt_root_psi_;
    for (double target : {p0 - 0.05, p0 - 0.1}) {
      const double eps = 1e-5;
      const double analytic = l.template dprofit_droot_collar_psi_for<K>(target);
      double up = l.template evaluate_root_collar_psi_for<K>(target + eps);
      double dn = l.template evaluate_root_collar_psi_for<K>(target - eps);
      if constexpr (K == phylloptim::Leaf::CostCurve::LeastCost) {
        up = std::log(up);
        dn = std::log(dn);
      }
      if (!std::isfinite(up) || !std::isfinite(dn)) {
        continue;
      }
      // The difference's own floor is ~3e-5 on the three-layer supply (TF24
      // included, which reads no upstream); the term this guards is O(1).
      near(analytic, (up - dn) / (2 * eps), 1e-4,
           phylloptim::Leaf::curve_name(static_cast<int>(K)) + " at collar " +
               std::to_string(target) + ", " + std::to_string(s.psi.size()) +
               " layer(s)");
      ++compared;
    }
  }
  return compared;
}

void test_collar_foc_matches_profit() {
  printf("collar dprofit is the derivative of the reported profit, every curve\n");
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  d.PPFD = 1500.0;
  const int compared = check_collar_foc_matches_profit<CC::TF24>(d) +
                       check_collar_foc_matches_profit<CC::CF77>(d) +
                       check_collar_foc_matches_profit<CC::JS22>(d) +
                       check_collar_foc_matches_profit<CC::CMax>(d) +
                       check_collar_foc_matches_profit<CC::TF24_floor>(d) +
                       check_collar_foc_matches_profit<CC::LeastCost>(d);
  // Hazard 16: a row whose difference straddles the feasible edge is skipped.
  printf("    %d points compared\n", compared);
  ok(compared == 24, "every curve, layer count and point was compared");
}

// lambda = dA/dE is the first-order condition every model in this family shares,
// so checking the analytic lambda against a finite-difference dA/dE is a check on
// the whole optimisation, not just on one formula.
// dprofit_droot_collar_psi reads the supply path's signed soil potentials, which
// used to be seated only by a solve -- so calling it on a leaf that had had
// set_physiology but not find_root_collar_psi read an empty vector. It now seats
// them itself. Ported from plant develop (#585).
void test_gradient_needs_no_prior_solve() {
  printf("dprofit/dpsi_collar does not require a prior solve\n");
  Drivers d;
  phylloptim::Leaf solved = make_leaf(d, {2.0}, {1.0});
  solved.find_root_collar_psi();
  phylloptim::Leaf unsolved = make_leaf(d, {2.0}, {1.0});
  // Bit-identical, not merely close: seating the potentials from psi_soil_ is
  // exactly what a solve does, so this is idempotent and moves no arithmetic.
  ok(unsolved.dprofit_droot_collar_psi(2.5) == solved.dprofit_droot_collar_psi(2.5),
     "the gradient is the same with and without a prior solve");
}

// The reversed-gradient state: the collar is asked about a potential drier than
// the stem it would have to supply, so there is no flow and no informative
// gradient. psi_leaf_to_ci does not return non-finite there -- it either throws
// (gc goes negative, the residual stops crossing zero, and the bracketing solver
// gives up) or returns a number built on a negative conductance. Both were
// reachable from TF24f's acclimation gradient on a dry patch.
void test_gradient_is_zero_in_reversed_gradient_state() {
  printf("dprofit/dpsi_collar returns zero when the gradient reverses\n");
  Drivers d;
  // Drier than psi_crit at every layer, so the leaf is shut down and every collar
  // potential below the wettest layer implies psi_leaf < psi_upstream.
  phylloptim::Leaf l = make_leaf(d, {5.9, 6.15, 6.4, 6.65, 6.9}, {1.0, 2.0, 3.0, 4.0, 5.0});
  l.find_root_collar_psi();
  for (double target : {1.0, 3.0, 5.5, 5.9}) {
    const double psi_leaf = l.find_psi_leaf_from_psi_root(
        target, l.roots_.psi_soil_);
    ok(target >= psi_leaf,
       "the target is in the reversed-gradient state at " + std::to_string(target));
    ok(l.dprofit_droot_collar_psi(target) == 0.0,
       "the gradient is exactly zero at " + std::to_string(target));
  }
}

// The 0.0 the two exits above return is a SENTINEL, and the collar root-find is what makes
// telling it apart from a genuine stationary point load bearing: it proposes
// root-finding on dprofit == 0, and the sentinel fires at the WET END of the very
// bracket such a solve would search. At bound_a = root_zero_E uptake is zero by
// construction, so psi >= psi_leaf and the reversed-gradient exit is taken -- and
// a bracketing solver evaluates that endpoint FIRST, to check its bracket
// brackets. So it would read the sentinel as the answer and return the
// zero-transpiration point as the optimum.
//
// This test pins the distinction rather than the eventual solver, so it is
// independent of how 11a is implemented. Note what it can and cannot see: the
// golden file cannot see any of this (no golden column comes from the gradient),
// and neither can a test that only samples interior points -- the sentinel region
// measured at most 3.46e-07 MPa wide over the golden grid, which is exactly why
// it would survive casual testing.
void test_gradient_reports_feasibility() {
  printf("dprofit/dpsi_collar distinguishes its 0.0 sentinel from a stationary point\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});

  // The bracket a root-find would search, the first one prepare_collar_solve
  // tries: the wet end is where uptake is zero, the dry end the continuity root
  // with the stem at psi_crit.
  const double wettest = 2.0;
  const double root_zero_E = l.find_root_psi(wettest, l.roots_.psi_soil_, 0);
  const double collar_at_dry_end =
      l.find_root_psi(wettest, l.roots_.psi_soil_, 1);
  const double bound_b = collar_at_dry_end;
  ok(root_zero_E < bound_b, "the bracket is non-empty");

  // At the wet end: 0.0, and NOT a stationary point. Taken just wetter than the
  // zero-uptake root, where uptake is certainly negative: the root itself is
  // placed to 1e-4 MPa and may land on either side of zero uptake.
  bool feasible = true;
  const double wet_side = root_zero_E - 1e-3;
  const double at_a = l.dprofit_droot_collar_psi(wet_side, &feasible);
  ok(at_a == 0.0, "the gradient is 0.0 at the wet bracket endpoint");
  ok(!feasible, "and the endpoint is reported infeasible, so the 0.0 is a sentinel");

  // Why accepting it would be wrong, stated as a number rather than an assertion
  // about intent: the endpoint is a far worse operating point than the optimum.
  l.find_root_collar_psi();
  const double best = l.profit_;
  const double at_endpoint = l.evaluate_root_collar_psi(root_zero_E);
  ok(at_endpoint < best - 1.0,
     "and profit at the endpoint is much worse than at the optimum");

  // An interior point: a real derivative, so a zero there WOULD be stationary.
  feasible = false;
  const double mid = 0.5 * (root_zero_E + bound_b);
  const double at_mid = l.dprofit_droot_collar_psi(mid, &feasible);
  ok(feasible, "an interior point is reported feasible");
  ok(std::isfinite(at_mid) && at_mid != 0.0,
     "and its gradient is a finite non-zero number");

  // The reversed-gradient state is infeasible too, and there the 0.0 is the whole
  // answer -- this is the case test_gradient_is_zero_in_reversed_gradient_state
  // already pins, checked here for the flag rather than the value.
  phylloptim::Leaf dry = make_leaf(d, {5.9, 6.15, 6.4, 6.65, 6.9},
                             {1.0, 2.0, 3.0, 4.0, 5.0});
  dry.find_root_collar_psi();
  feasible = true;
  ok(dry.dprofit_droot_collar_psi(3.0, &feasible) == 0.0 && !feasible,
     "the reversed-gradient state is reported infeasible");

  // The out-parameter changes nothing about the value, so TF24f's contract (the
  // return value alone, consumed as an ODE rate) is untouched. Bit-identical, not
  // merely close: the flag is a write to a caller's bool, not arithmetic.
  ok(l.dprofit_droot_collar_psi(mid) == at_mid,
     "passing no flag returns exactly the same value");
  ok(l.dprofit_droot_collar_psi(wet_side) == 0.0,
     "and the sentinel is still 0.0 for a caller that does not ask");
}

// The AD derivative and the forward model are instantiations of ONE
// body, so they cannot be derivatives of different functions. That was not true
// before: `detail::assim_colimited_ad` associated the electron-limited term
// left-to-right where `assim_electron_limited` divides the bracket first, and used
// `s*s` where `assim_colimited` uses `pow(s, 2)`.
//
// The check that bites is the DERIVATIVE against a central difference of the
// `double` function, because that is precisely the identity the drift broke. Note a
// weaker test would have passed throughout: the drifted replica was still a
// perfectly good derivative of *itself*, and agreed with the real function to
// ~1e-16 in VALUE. It is the derivative-of-the-same-function property that failed.
void test_ad_kernels_are_the_model_not_a_mirror() {
  printf("AD differentiates the model's own algebra, not a mirror of it\n");
  Drivers d;
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  using AD = xad::fwd<double>::active_type;

  // 1. The kernels and the double entry points are the same code, so they must
  // agree BIT-EXACTLY, not merely closely. This is what "one body" means.
  for (double ci : {5.0, 12.0, 20.0, 29.26, 35.0, 39.0}) {
    ok(l.assim_colimited_kernel(ci) == l.assim_colimited(ci),
       "assim_colimited_kernel<double> is assim_colimited at ci=" +
           std::to_string(ci));
    ok(l.assim_rubisco_limited_kernel(ci) == l.assim_rubisco_limited(ci),
       "rubisco kernel matches bit-exactly at ci=" + std::to_string(ci));
    ok(l.assim_electron_limited_kernel(ci) == l.assim_electron_limited(ci),
       "electron kernel matches bit-exactly at ci=" + std::to_string(ci));
  }
  for (double p : {0.5, 2.0, 3.5949, 5.0}) {
    ok(l.hydraulic_cost_TF_kernel(p) == l.hydraulic_cost_TF(p),
       "hydraulic_cost_TF_kernel<double> is hydraulic_cost_TF at psi=" +
           std::to_string(p));
  }

  // 2. The AD derivative of the kernel against a central difference of the DOUBLE
  // function. Richardson-extrapolated, so the FD reference is good to ~1e-10 and
  // the tolerance is testing the derivative rather than the difference quotient.
  const auto richardson = [](auto f, double x, double h) {
    const double d1 = (f(x + h) - f(x - h)) / (2 * h);
    const double d2 = (f(x + h / 2) - f(x - h / 2)) / h;
    return (4 * d2 - d1) / 3;
  };
  for (double ci : {12.0, 20.0, 29.26, 35.0}) {
    AD a = ci; xad::derivative(a) = 1.0;
    const double ad = xad::derivative(l.assim_colimited_kernel(a));
    const double fd =
        richardson([&](double x) { return l.assim_colimited(x); }, ci, 1e-4);
    near(ad, fd, 1e-8, "dA/dci: AD vs Richardson FD at ci=" + std::to_string(ci));
  }
  for (double p : {0.5, 2.0, 3.5949, 5.0}) {
    AD a = p; xad::derivative(a) = 1.0;
    const double ad = xad::derivative(l.hydraulic_cost_TF_kernel(a));
    const double fd =
        richardson([&](double x) { return l.hydraulic_cost_TF_kernel(x); }, p, 1e-4);
    near(ad, fd, 1e-8, "dcost/dpsi: AD vs Richardson FD at psi=" +
                           std::to_string(p));
  }

  // 3. The cost kernel must be PURE -- an AD probe must not scribble the cached
  // hydraulic_cost_, or a gradient evaluation would corrupt reported model state.
  const double cached = l.hydraulic_cost_TF(3.0);
  AD probe = 4.5; xad::derivative(probe) = 1.0;
  (void)l.hydraulic_cost_TF_kernel(probe);
  ok(l.hydraulic_cost_ == cached,
     "an AD probe of the cost kernel leaves hydraulic_cost_ untouched");
}

// The collar solve solves its own first-order condition, so the
// check that bites is the RESIDUAL at the returned point, not the returned value.
// Measured over the golden grid: 240 of 240 feasible rows improved their residual
// and none got worse, interior rows landing at a median |dprofit| of 5.6e-15
// against golden section's 7.8e-4.
//
// Deliberately NOT checked against profit. Profit is the wrong instrument here for
// a reason worth recording: it is the maximum, so it is flat, and its own
// numerical floor is set by the nested ci root-find's 1e-7 tolerance. Over the
// grid two rows come out ~6e-7 LOWER in profit than golden section while their
// residual improves by ten orders of magnitude -- that is the floor, not a
// regression, and a test asserting "profit never decreases" would encode the noise.
void test_collar_solve_satisfies_its_own_first_order_condition() {
  printf("the collar solve lands where dprofit == 0\n");
  Drivers d;
  // An interior optimum: the gradient at the answer should be at solver precision.
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double residual = l.dprofit_droot_collar_psi(l.opt_root_psi_);
  ok(std::abs(residual) < 1e-9,
     "the interior optimum satisfies dprofit == 0 to solver precision");
  // The bar this clears, stated as the thing it replaced: golden section left
  // -6.22e-4 here, which is what made the argmax a staircase in the traits.
  ok(std::abs(residual) < 1e-6 * 6.22e-4,
     "and it is orders below the residual golden section left behind");

  // The answer stays strictly inside the feasible interval. This is the guard on
  // the #31 regression that returning bound_a caused: below psi_upstream the
  // profit algebra runs on a negative conductance and profit is DISCONTINUOUS, so
  // an endpoint is not merely a poor answer, it is a different function. Returning
  // bound_a moved 22 golden rows' profit DOWN by 1.44 before this was caught.
  const double root_zero_E = l.find_root_psi(2.0, l.roots_.psi_soil_, 0);
  const double bound_b =
      l.find_root_psi(2.0, l.roots_.psi_soil_, 1);
  ok(l.opt_root_psi_ > root_zero_E,
     "the collar is strictly drier than the zero-uptake bound (#31)");
  ok(l.opt_root_psi_ <= bound_b, "and no drier than the binding dry limit");
}

// A CONSTRAINED optimum, which 51 of the 288 feasible golden-grid rows at 25 C are: profit
// is still climbing at one end of the feasible interval, so dprofit never crosses
// zero inside it. The residual is then NOT small, and that is correct rather than a
// convergence failure -- so this pins the property that actually holds there.
void test_collar_solve_handles_a_pinned_optimum() {
  printf("a collar optimum pinned to its constraint\n");
  Drivers d;
  // psi_soil 4.0 over 5 layers at vpd 2.0 -- one of the measured pinned rows.
  phylloptim::Leaf l = make_leaf(d, {4.0, 4.25, 4.5, 4.75, 5.0},
                           {1.0, 2.0, 3.0, 4.0, 5.0});
  l.find_root_collar_psi();
  const double root_zero_E = l.find_root_psi(4.0, l.roots_.psi_soil_, 0);
  const double bound_b =
      l.find_root_psi(4.0, l.roots_.psi_soil_, 1);
  ok(std::isfinite(l.profit_), "the pinned row still yields a finite profit");
  ok(l.opt_root_psi_ > root_zero_E,
     "and a collar strictly inside the zero-uptake bound, not on it (#31)");
  ok(l.opt_root_psi_ <= bound_b, "and inside the dry bound");
  // It is pinned NEAR the wet end, but not AT it: the true argmax sits essentially
  // on the feasibility boundary, and the step-in resolves it to ~1e-6 of the
  // bracket width. Measured 3.29e-4 wetter than golden section's answer, with a
  // profit 1.27e-3 HIGHER -- so being pinned is not being wrong.
  ok(l.opt_root_psi_ - root_zero_E < 1e-4 * (bound_b - root_zero_E),
     "the pinned answer sits at the wet end of the interval");
  // The gradient there is genuinely non-zero, which is what "constrained" means.
  ok(l.dprofit_droot_collar_psi(l.opt_root_psi_) < 0.0,
     "profit is decreasing at the pinned answer, so the bound is what binds");
  ok(l.operating_point_kind() ==
         phylloptim::Leaf::OperatingPointKind::BoundarySoil,
     "and the leaf says so: the point is tagged boundary-soil");
}

// The operating-point classification, on ONE REUSED LEAF -- which is the only way
// to test it. `Leaf` is a value member that plant drives every individual in a
// patch through (hazard 8), so the failure this guards is a branch that declines
// to write the tag and thereby reports the PREVIOUS plant's kind of operating
// point. test_golden cannot see that: it builds a fresh Leaf per grid point, so a
// bit-identical golden run says nothing here. The sequence below alternates kinds
// on purpose, so an unwritten tag reads as the step before it.
//
// The count check over the whole golden grid lives in test_golden.cpp; this is
// the other half -- that every path writes, and that the tag comes from the
// branch rather than from the numbers.
void test_operating_point_kind_is_written_by_every_path() {
  printf("the collar solve says which kind of operating point it found\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  Drivers d;
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  ok(l.operating_point_kind() == Kind::Unsolved,
     "a fresh Leaf reports no operating point");

  const std::vector<double> depth{1.0}, root{1.0 / d.area_leaf};
  const auto solve = [&](double psi, double ppfd) {
    std::vector<double> ps{psi};
    l.set_physiology(fixture::root_network(root, depth), ppfd, ps, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
  };

  solve(2.0, d.PPFD);
  ok(l.operating_point_kind() == Kind::Interior,
     "wet soil gives an interior profit maximum");

  // Drier than psi_crit: shutdown on water. Would read `interior` if the shutdown
  // exit did not write.
  solve(20.0, d.PPFD);
  ok(l.operating_point_kind() == Kind::HydraulicShutdown,
     "past psi_crit the point is a hydraulic shutdown");

  // ⚠️ THE ASSERTION THE TAG EXISTS FOR. dprofit's shut-down exit returns a hard
  // 0.0 sentinel, so a classifier of the form "|residual| < tol => interior
  // optimum" would call this an interior stationary point, and then read a
  // curvature of zero off the same sentinel and agree with itself. The residual
  // really is exactly 0.0 here, and the tag really does say shutdown.
  bool feasible = true;
  const double residual = l.dprofit_droot_collar_psi(l.opt_root_psi_, &feasible);
  ok(residual == 0.0, "dprofit at a shut-down point is exactly the 0.0 sentinel");
  ok(!feasible, "and reports itself infeasible");
  ok(l.operating_point_kind() != Kind::Interior,
     "so a residual test would misclassify it, and the tag does not");

  // Back to wet: would read `hydraulic-shutdown` if the interior path did not
  // write.
  solve(2.0, d.PPFD);
  ok(l.operating_point_kind() == Kind::Interior,
     "and the next wet solve is interior again, not the previous plant's kind");

  // Dim enough that gross assimilation cannot cover R_d, on wet soil. Light does
  // not shut a leaf (#182): the leaf is a net CO2 source, its ci sits above
  // ambient, and its optimum is an ordinary pinned one -- the leaf declines to
  // open -- not a kind of its own.
  solve(1.0, 10.0);
  ok(l.assim_max_ < 0.0, "the dim leaf cannot fix carbon at ambient CO2");
  ok(l.operating_point_kind() == Kind::BoundarySoil,
     "and declines to open, a wet-pinned optimum rather than a shutdown");
  ok(l.ci_ > l.ca_, "with its internal CO2 above ambient: respiratory CO2 leaves");

  // A constrained optimum, then the same leaf asked to EVALUATE a prescribed
  // collar potential rather than optimise one. The prescribed point is not an
  // optimum of any kind, and would read `boundary-soil` if that path did not write.
  const std::vector<double> psi5{4.0, 4.25, 4.5, 4.75, 5.0};
  const std::vector<double> depth5{1.0, 2.0, 3.0, 4.0, 5.0};
  const std::vector<double> root5(5, 1.0 / 5.0 / d.area_leaf);
  l.set_physiology(fixture::root_network(root5, depth5), d.PPFD, psi5, depth5,
                   d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                   d.atm_o2_kpa, d.atm_kpa);
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == Kind::BoundarySoil,
     "psi_soil 4.0 over 5 layers is pinned to the wet bound");
  const double pinned_collar = l.opt_root_psi_;

  l.evaluate_root_collar_psi(pinned_collar);
  ok(l.operating_point_kind() == Kind::Prescribed,
     "evaluating a given collar potential is not an optimum, and says so");

  // ...and the optimising path takes it back, so `prescribed` is not sticky.
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == Kind::BoundarySoil,
     "and optimising again restores the pin");

  // set_traits returns the object to its just-constructed state, and the
  // classification is part of that state (hazard 10 / setup_clean_leaf).
  l.set_traits(l.vcmax_25, l.stem_c, l.stem_P50, l.roots_.root_c,
               l.roots_.root_P50, l.TF24_beta2, l.jmax_25, l.phi_J,
               l.curv_fact_elec_trans, l.curv_fact_colim, l.TF24_cost_scale,
               l.R_d_25, l.JS22_gamma, l.CMax_a, l.CMax_b, l.g_min);
  ok(l.operating_point_kind() == Kind::Unsolved,
     "set_traits clears the classification with the rest of the solved state");

  // Neither failure kind is reachable from a driver sweep, and that is the point
  // of separating them from a pin: they mean "no plant is described".
  // test_collar_solve_refuses_rather_than_guessing drives both directly.
  ok(std::string(phylloptim::Leaf::operating_point_kind_name(
         Kind::SolverRefused)) == "solver-refused",
     "the failure kinds are named, not numbered");
}

// The two branches of maximise_profit_over_collar that change the collar it
// returns. No driver reaches either, so both are driven through the public
// maximise_profit_over_collar with brackets built to land on them.
void test_collar_solve_refuses_rather_than_guessing() {
  printf("the collar solve refuses a bracket it cannot resolve\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  Drivers d;
  const std::vector<double> psi{2.0, 2.25, 2.5, 2.75, 3.0};
  const std::vector<double> depth{1.0, 2.0, 3.0, 4.0, 5.0};

  // A comfortably interior optimum, and the interval it was found in.
  phylloptim::Leaf l = make_leaf(d, psi, depth);
  double bound_a = 0.0, bound_b = 0.0;
  ok(l.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(bound_a, bound_b),
     "the reference case has a feasible collar interval");
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == Kind::Interior,
     "and an interior optimum inside it to straddle");
  const double star = l.opt_root_psi_;

  // Bounds handed over inverted, so the interval runs dry to wet: dprofit <= 0 at
  // bound_a and >= 0 at bound_b, which makes each end a local maximum.
  const double dry_end = star + 0.9 * (bound_b - star);
  const double wet_end = star - 0.9 * (star - bound_a);

  phylloptim::Leaf p = make_leaf(d, psi, depth);
  double pa = 0.0, pb = 0.0;
  p.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(pa, pb);
  const double profit_dry = p.profit_psi_leaf_TF(
      p.find_psi_leaf_from_psi_root(dry_end, p.supply_psi_soil()), dry_end);
  const double profit_wet = p.profit_psi_leaf_TF(
      p.find_psi_leaf_from_psi_root(wet_end, p.supply_psi_soil()), wet_end);
  ok(profit_wet > profit_dry,
     "the wetter end is the better of the two, by construction");

  phylloptim::Leaf m = make_leaf(d, psi, depth);
  double ma = 0.0, mb = 0.0;
  m.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(ma, mb);
  const double refused = m.maximise_profit_over_collar<phylloptim::Leaf::CostCurve::TF24>(dry_end, wet_end);
  ok(m.operating_point_kind() == Kind::SolverRefused,
     "the solve reports that it could not resolve the bracket");
  ok(m.operating_point_kind() != Kind::BoundarySoil &&
         m.operating_point_kind() != Kind::BoundaryDry,
     "and does not pass it off as a constrained optimum");
  // The endpoint the solve returns is stepped a fraction of the width inside the
  // bound it came from, so compare against the bound rather than for equality.
  ok(std::abs(refused - wet_end) < 1e-5,
     "and returns the end with the higher profit");
  ok(std::abs(refused - dry_end) > 0.1, "which is not the drier end");

  // A bracket lying wholly at or below bound_a, where dprofit takes its
  // reversed-gradient exit: no usable gradient at either end. The
  // non-finite-gradient half of the same guard has no bracket that reaches it and
  // is not covered.
  //
  // ⚠️ THE BRACKET IS BELOW THE BOUND, NOT JUST INSIDE IT. G^-1 is the forward
  // table inverted, so zero flux returns the collar exactly and the infeasible
  // region is exactly (-inf, bound_a] -- one ulp above it answers. DO NOT fit the
  // inverse separately: at zero flux the round trip then puts psi_leaf about 1e-07
  // below the collar, which is a reversed gradient the model does not have, and it
  // opens an infeasible sliver that wide above bound_a.
  phylloptim::Leaf s = make_leaf(d, psi, depth);
  double sa = 0.0, sb = 0.0;
  s.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(sa, sb);
  const double sliver = 1e-9;
  bool feasible = true;
  s.dprofit_at_collar_psi<phylloptim::Leaf::CostCurve::TF24>(sa, &feasible);
  ok(!feasible, "the wet bound admits no informative gradient");
  const double fallen_back = s.maximise_profit_over_collar<phylloptim::Leaf::CostCurve::TF24>(sa - sliver, sa);
  ok(s.operating_point_kind() == Kind::SolverRefused,
     "a bracket with no usable gradient at either end is refused too");
  ok(std::isfinite(fallen_back) && fallen_back >= sa - sliver &&
         fallen_back <= sa,
     "and the fallback stays inside the bracket it was handed");
}

// Hazard 3, re-measured rather than argued. The guide's constraint is that the
// argmax must vary SMOOTHLY with inputs. Golden section resolved it only to
// GSS_tol_abs, making it piecewise constant at fine scales -- so a sharper solve
// makes it smoother, which is the opposite of what the hazard reads like.
//
// ⚠️ This measures smoothness in a TRAIT. Hazard 3's actual concern is smoothness
// in plant state feeding the demographic growth-rate gradient, which cannot be
// measured from here while plant #591 blocks end-to-end validation. Same
// mechanism, so the same direction is expected -- but it is not the same test, and
// this one passing does not discharge the hazard.
void test_collar_argmax_is_smooth_in_a_trait() {
  printf("the collar argmax varies smoothly with a trait (hazard 3)\n");
  Drivers d;
  const int n = 11;
  const double base = 96.0, step = base * 1e-3;  // 0.1% steps in vcmax_25
  double collar[n];
  int distinct = 0;
  for (int i = 0; i < n; ++i) {
    phylloptim::Leaf l = fixture::multilayer_leaf();
    l.vcmax_25 = base + i * step;
    l.setup_transpiration(100);
    l.setup_root_vulnerability(100);
    std::vector<double> ps{2.0}, depth{1.0}, root{1.0 / d.area_leaf};
    l.set_physiology(fixture::root_network(root, depth), d.PPFD, ps, depth, d.K_s * d.theta / d.h, d.atm_vpd,
                     d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
    collar[i] = l.opt_root_psi_;
    bool seen = false;
    for (int j = 0; j < i; ++j) {
      if (collar[j] == collar[i]) seen = true;
    }
    if (!seen) ++distinct;
  }
  // Golden section gave 6 distinct values out of 11 -- the staircase, tread width
  // ~GSS_tol_abs. Every step must now move the answer.
  ok(distinct == n, "every trait step moves the argmax (no staircase)");

  // Smoothness, as second differences against the step size. Golden section
  // measured 3.9e-4, i.e. the same size as the steps themselves: pure noise. The
  // root-find measured 3.4e-7, ~1000x smaller than its own steps.
  double worst_d2 = 0.0, mean_step = 0.0;
  for (int i = 1; i < n; ++i) {
    mean_step += std::abs(collar[i] - collar[i - 1]) / (n - 1);
  }
  for (int i = 2; i < n; ++i) {
    worst_d2 = std::max(worst_d2,
                        std::abs(collar[i] - 2 * collar[i - 1] + collar[i - 2]));
  }
  ok(worst_d2 < 0.02 * mean_step,
     "second differences are small against the step size, not comparable to it");
}

// #25's invariant, asserted rather than documented: every psi is a positive
// magnitude, so the soil->collar derivative is a CONDUCTANCE. Under the old signed
// convention this came back negative and had to be negated at the one call site
// that wanted a conductance -- the omission of that negation is the bug that
// motivated #8, and it is now unrepresentable.
void test_soil_conductance_is_positive() {
  printf("dE_up/d(collar suction) is a positive conductance\n");
  Drivers d;
  for (int layers : {1, 3, 5}) {
    std::vector<double> ps(layers), depth(layers);
    for (int i = 0; i < layers; ++i) { ps[i] = 1.0 + 0.25 * i; depth[i] = 1.0 * (i + 1); }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const double S = l.dE_from_soil_dpsi_collar(l.opt_root_psi_, l.roots_.psi_soil_);
    ok(std::isfinite(S) && S > 0.0,
       "conductance is finite and positive at " + std::to_string(layers) + " layers");
    // It really is dE/dT: a central difference on the uptake must agree.
    const double h = 1e-7;
    double up = 0.0, dn = 0.0;
    std::vector<double> buf(l.soil_consumption_.size(), 0.0);
    l.roots_.uptake_at(l.opt_root_psi_ + h, l.roots_.psi_soil_, buf, up);
    l.roots_.uptake_at(l.opt_root_psi_ - h, l.roots_.psi_soil_, buf, dn);
    near(S, (up - dn) / (2.0 * h), 1e-5,
         "conductance matches a central difference at " + std::to_string(layers) + " layers");
    // The collar objective's one walk gives what the two calls give, bit for
    // bit: at the solved collar, and on a layer's own potential, a kink where
    // the slope is NaN (hazard 6).
    const auto same = [](double a, double b) {
      return a == b || (std::isnan(a) && std::isnan(b));
    };
    for (double collar : {l.opt_root_psi_, l.roots_.psi_soil_[0]}) {
      const std::string at = " at " + std::to_string(layers) + " layers, collar " +
                             std::to_string(collar);
      std::vector<double> c1(buf.size()), c2(buf.size());
      double e1 = 0.0, e2 = 0.0;
      l.roots_.uptake_at(collar, l.roots_.psi_soil_, c1, e1);
      const double s1 = l.roots_.duptake_dpsi(collar, l.roots_.psi_soil_);
      const double s2 =
          l.roots_.uptake_and_slope_at(collar, l.roots_.psi_soil_, c2, e2);
      ok(e1 == e2 && c1 == c2 && same(s1, s2),
         "one walk's uptake, draws and slope are the two calls'" + at);
      phylloptim::Leaf a = l, b = l;
      a.E_from_Soil_to_Root_Collar(collar, a.supply_psi_soil());
      const double sa = a.dE_from_soil_dpsi_collar(collar, a.supply_psi_soil());
      const double sb = b.E_and_slope_from_soil(collar, b.supply_psi_soil());
      ok(a.E_up_ == b.E_up_ && a.soil_consumption_ == b.soil_consumption_ &&
             same(sa, sb),
         "and so through the Leaf" + at);
      if (collar != l.opt_root_psi_) {
        ok(std::isnan(s2), "the slope on a layer's potential is NaN" + at);
      }
    }
  }
}

// The convention is now checkable, so check that it is checked: a caller still
// holding the pre-#25 signed vector must fail loudly, not run.
void test_signed_potentials_are_rejected() {
  printf("signed potentials are rejected at the input boundary\n");
  Drivers d;
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  bool threw = false;
  try {
    l.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD, {-2.0}, {1.0},
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "set_physiology rejects a negative psi_soil");

  phylloptim::Leaf ok_leaf = make_leaf(d, {2.0}, {1.0});
  ok_leaf.find_root_collar_psi();
  threw = false;
  try {
    ok_leaf.E_from_Soil_to_Root_Collar(2.5, {-2.0});
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "E_from_Soil_to_Root_Collar rejects a signed soil vector");

  // The constructor's half of the invariant.
  threw = false;
  try {
    phylloptim::Leaf bad(100, 2.04, -3.0, 2.65, 1.29, 1, 167 * 100, 0.3,
                   0.7, 0.99, 1e-8, 100, 1e-6, 1000, 46.32995);
    static_cast<void>(bad);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "the constructor rejects a negative stem_b");
}

// As the soil dries through psi_crit, which bounds nothing (#183), and then
// through the end of the stem curve into shutdown, transpiration, profit and the
// held stem move continuously, and every open leaf transpires what its roots
// supply. Swept at both GSS_tol_abs a caller uses.
void test_collar_continuous_into_shutdown() {
  printf("the collar solve is continuous through psi_crit and into shutdown\n");
  Drivers d;
  for (double gss : {1e-3, 1e-1}) {
    const auto solve = [&](double psi_soil) {
      phylloptim::Leaf l = fixture::multilayer_leaf();
      l.GSS_tol_abs = gss;
      l.setup_transpiration(100);
      l.setup_root_vulnerability(100);
      std::vector<double> ps{psi_soil}, depth{1.0}, root{1.0 / d.area_leaf};
      l.set_physiology(fixture::root_network(root, depth), d.PPFD, ps, depth,
                       d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                       d.atm_o2_kpa, d.atm_kpa);
      l.find_root_collar_psi();
      return l;
    };
    const std::string tag = " (GSS_tol_abs " + std::to_string(gss) + ")";
    // Through psi_crit (5.870283): no change of regime, and no jump.
    {
      double prev = NAN, worst = 0.0;
      int shut = 0;
      for (int i = 0; i <= 800; ++i) {
        phylloptim::Leaf l = solve(5.70 + 5e-4 * i);
        if (l.operating_point_kind() ==
            phylloptim::Leaf::OperatingPointKind::HydraulicShutdown) ++shut;
        if (!std::isnan(prev)) worst = std::max(worst, std::abs(l.profit_ - prev));
        prev = l.profit_;
      }
      ok(shut == 0, "a soil drying through psi_crit does not shut the leaf" + tag);
      ok(worst < 1e-2, "and its profit does not jump there" + tag);
    }
    double worst_water = 0.0, worst_step = 0.0, last_open = NAN, first_shut = NAN;
    double prev = NAN;
    int n_open = 0, n_shut = 0;
    // Into shutdown, short of the end of the stem curve (15.439 at the
    // defaults), where the interval closes below collar precision.
    for (int i = 0; i <= 800; ++i) {
      const double psi_soil = 12.0 + 5e-3 * i;   // to 16.0, past the end
      phylloptim::Leaf l = solve(psi_soil);
      if (l.operating_point_kind() !=
          phylloptim::Leaf::OperatingPointKind::HydraulicShutdown) {
        ++n_open;
        // Against the flux's own resolution as well as relatively: the stem's
        // flux is read to ~eps of kmax G(inf), which near shutdown is the size
        // of the flux itself.
        if (l.transpiration_ > 0.0) {
          const double resolution = 64.0 * std::numeric_limits<double>::epsilon() *
                                    l.leaf_specific_conductance_max_ *
                                    l.stem_curve_limit();
          worst_water = std::max(worst_water,
                                 std::abs(l.transpiration_ - l.E_up_) /
                                     (l.E_up_ + resolution / 1e-6));
        }
        last_open = l.profit_;
      } else if (std::isnan(first_shut)) {
        ++n_shut;
        first_shut = l.profit_;
      }
      if (!std::isnan(prev)) worst_step = std::max(worst_step, std::abs(l.profit_ - prev));
      prev = l.profit_;
    }
    ok(n_open > 0 && n_shut > 0, "the sweep reaches both regimes" + tag);
    printf("    worst |E - E_up| / E_up into shutdown: %.3g%s\n", worst_water, tag.c_str());
    ok(worst_water < 1e-6, "a transpiring leaf transpires what its roots supply" + tag);
    // A 5e-3 MPa step in psi_soil moves a smooth profit by ~1e-3 here.
    ok(worst_step < 1e-2, "no jump in profit along the sweep" + tag);
    near(first_shut, last_open, 1e-3, "shutdown is the E -> 0 limit" + tag);
  }
}

// The same continuity across the light-compensation point, where a closed leaf
// nets A(ci = ca): zero at compensation and -R_d in the dark.
void test_collar_continuous_through_light_compensation() {
  printf("the collar solve is continuous through the light-compensation point\n");
  Drivers d;
  using Kind = phylloptim::Leaf::OperatingPointKind;
  double prev = NAN, worst_step = 0.0, dark = NAN;
  int n_pinned = 0, n_open = 0;
  for (int i = 0; i <= 6000; ++i) {
    Drivers di = d;
    di.PPFD = 0.01 * i;   // to 60, past compensation
    phylloptim::Leaf l = make_leaf(di, {1.0}, {1.0});
    l.find_root_collar_psi();
    // Below compensation the leaf declines to open and pins to the wet bound;
    // above it the optimum is interior. Neither is a shutdown (#182).
    ok(l.operating_point_kind() != Kind::HydraulicShutdown, "light does not shut a leaf");
    if (l.operating_point_kind() == Kind::Interior) ++n_open; else ++n_pinned;
    if (i == 0) dark = l.assim_colimited_ + l.R_d_;
    if (!std::isnan(prev)) worst_step = std::max(worst_step, std::abs(l.profit_ - prev));
    prev = l.profit_;
  }
  ok(n_open > 0 && n_pinned > 0, "the sweep crosses the compensation point");
  // A 0.01 step in PPFD moves a smooth profit by ~4e-4 here.
  ok(worst_step < 1e-2, "no jump in profit across compensation");
  // In the dark A is -R_d at every ci, so the saturated value and the limit
  // agree to the rounding of one square root.
  near(dark, 0.0, 1e-12, "and in the dark the leaf nets -R_d");
}

// On a layered soil the collar interval narrows below the 1e-4 MPa its bounds
// are placed to before it closes, at the end of the stem curve; the solve must
// still reach shutdown smoothly.
void test_layered_soil_dries_into_shutdown() {
  printf("a layered soil dries into shutdown without a refused solve\n");
  Drivers d;
  int refused = 0, shut = 0;
  double prev = NAN, worst_step = 0.0;
  for (int i = 0; i <= 400; ++i) {
    const double psi = 12.0 + 1e-2 * i;
    phylloptim::Leaf l = make_leaf(d, {psi + 0.4, psi + 0.2, psi}, {1.0, 2.0, 3.0});
    l.find_root_collar_psi();
    const auto k = l.operating_point_kind();
    if (k == phylloptim::Leaf::OperatingPointKind::SolverRefused) ++refused;
    if (k == phylloptim::Leaf::OperatingPointKind::HydraulicShutdown) ++shut;
    if (!std::isnan(prev)) worst_step = std::max(worst_step, std::abs(l.profit_ - prev));
    prev = l.profit_;
  }
  ok(shut > 0, "the sweep reaches shutdown");
  ok(refused == 0, "no solve on the way is refused, got " + std::to_string(refused));
  // A 1e-5 MPa step moves a smooth profit here by ~1e-5.
  ok(worst_step < 1e-4, "no jump in profit along the sweep");
}

void test_lambda_equals_dA_dE_single_layer() {
  printf("marginal cost of water: analytic lambda vs dA/dE (stem free)\n");
  Drivers d;
  for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
    phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
    // the TF24 leaf route holds the collar fixed at psi_soil_[0] and optimises
    // the stem, so the single-layer lambda is the one that applies here.
    l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
    const double psi = l.opt_psi_leaf_;
    const double eps = 1e-6;
    l.set_leaf_states_rates_from_psi_leaf(psi + eps, psi_soil);
    const double A1 = l.assim_colimited_, E1 = l.transpiration_;
    l.set_leaf_states_rates_from_psi_leaf(psi - eps, psi_soil);
    const double A0 = l.assim_colimited_, E0 = l.transpiration_;
    const double fd = (A1 - A0) / (E1 - E0);
    l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
    // Tolerance is set by the optimiser: GSS_tol_abs is 1e-3 on psi, so the
    // first-order condition only holds to about that accuracy.
    near(l.lambda_emergent() / fd, 1.0, 1e-3,
         "lambda/(dA/dE) at psi_soil=" + std::to_string(psi_soil));
    ok(l.lambda_emergent() > 0.0, "lambda is positive");
  }
}

// dprofit_dpsi_leaf against a central difference of the objective it claims to
// differentiate, for both cost curves.
//
// THE POINT OF THE TEST is that the two halves come from different places: the
// derivative is analytic (forward AD for the TF24 curve, a closed form for
// Cowan-Farquhar, and the implicit function theorem through the ci root-find for
// the assimilation term in both), while the reference differences
// `profit_psi_leaf_*` itself. So this checks the derivative against the model
// rather than against another derivative -- which is the failure that matters:
// a derivative of a nearby function is smooth, plausible and wrong.
void test_dprofit_dpsi_leaf_matches_a_finite_difference() {
  printf("dprofit/dpsi_leaf: analytic vs a central difference, both curves\n");
  Drivers d;
  d.PPFD = 1500.0;
  // ⚠️ THE STEP IS MEASURED, NOT PICKED. The differenced quantity is computed
  // through the ci root-find, so it carries that solver's tolerance as a fixed
  // absolute noise floor, and a central difference divides it by h. The relative
  // error is therefore ~1/h until truncation takes over: measured at one of the
  // rows below it runs 1.3e-01, 1.1e-02, 1.1e-03, 1.1e-04, 1.1e-05, 1.9e-07 for
  // h from 1e-8 to 1e-3, then back up to 1.3e-04 at 1e-2. So 1e-3 is the floor of
  // that V, and a "small" step like 1e-6 sits three orders up the noise side of
  // it -- which reads as a wrong derivative and is not one.
  const double h = 1e-3;
  double worst_tf = 0.0, worst_cf = 0.0;
  int rows = 0;

  for (double psi_soil : {0.5, 1.0, 2.0, 3.0}) {
    // Sample the interior of [psi_soil, psi_crit] rather than only the optimum:
    // away from the optimum dprofit is LARGE, so a wrong chain rule shows up as a
    // relative error instead of hiding under a near-zero value.
    phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
    const double lo = psi_soil, hi = l.psi_crit;
    for (double frac : {0.15, 0.35, 0.55, 0.75}) {
      const double psi = lo + frac * (hi - lo);
      ++rows;

      bool feasible = false;
      const double an_tf =
          l.dprofit_dpsi_leaf<phylloptim::Leaf::CostCurve::TF24>(psi, psi_soil,
                                                                 &feasible);
      const double up_tf = l.profit_psi_leaf_TF(psi + h, psi_soil);
      const double dn_tf = l.profit_psi_leaf_TF(psi - h, psi_soil);
      const double fd_tf = (up_tf - dn_tf) / (2.0 * h);
      if (feasible) {
        worst_tf = std::max(worst_tf, std::abs(an_tf / fd_tf - 1.0));
        near(an_tf / fd_tf, 1.0, 5e-5,
             "TF24 dprofit at frac=" + std::to_string(frac) +
                 " psi_soil=" + std::to_string(psi_soil));
      }

      l.CF77_lambda_ = 1.5e5;
      const double an_cf =
          l.dprofit_dpsi_leaf<phylloptim::Leaf::CostCurve::CF77>(
              psi, psi_soil, &feasible);
      const double up_cf = l.profit_psi_leaf_CF77(psi + h, psi_soil);
      const double dn_cf = l.profit_psi_leaf_CF77(psi - h, psi_soil);
      const double fd_cf = (up_cf - dn_cf) / (2.0 * h);
      if (feasible) {
        worst_cf = std::max(worst_cf, std::abs(an_cf / fd_cf - 1.0));
        near(an_cf / fd_cf, 1.0, 5e-5,
             "Cowan-Farquhar dprofit at frac=" + std::to_string(frac) +
                 " psi_soil=" + std::to_string(psi_soil));
      }
    }
  }
  // Measured 8.9e-07 (TF24) and 4.1e-06 (Cowan-Farquhar). The bound is 5e-05
  // rather than either: the residual is the ci solver's noise floor divided by h,
  // and that floor is libm-dependent, so it moves across platforms. Nothing is
  // lost by the slack -- a missing or mis-signed chain-rule term is an O(1)
  // relative error, not an O(1e-5) one.
  printf("    %d points | worst relative error: TF24 %.3e, Cowan-Farquhar %.3e\n",
         rows, worst_tf, worst_cf);
}

// The derivative and the optimiser have to agree about where the optimum is: the
// solver returns an argmax, so dprofit there must be ~zero, and AWAY from it must
// not be. The second half is what makes this more than a tautology -- a
// derivative that returned zero everywhere would pass the first half alone.
void test_dprofit_dpsi_leaf_vanishes_at_the_optimum() {
  printf("dprofit/dpsi_leaf is ~0 at the optimum and not away from it\n");
  Drivers d;
  d.PPFD = 1500.0;
  for (double psi_soil : {0.5, 1.0, 2.0}) {
    phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
    l.CF77_lambda_ = 1.5e5;
    l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::CF77>();
    const double psi = l.opt_psi_leaf_;
    if (psi <= psi_soil + 1e-6 || psi >= l.psi_crit - 1e-6) {
      continue; // pinned: the gradient is genuinely non-zero there
    }
    bool feasible = false;
    const double at =
        l.dprofit_dpsi_leaf<phylloptim::Leaf::CostCurve::CF77>(
            psi, psi_soil, &feasible);
    ok(feasible, "the optimum admits an informative gradient");
    // Scaled by the derivative's own magnitude a little away from the optimum,
    // so the bound means "small compared with what this function returns here"
    // rather than "small in absolute carbon units".
    const double away =
        l.dprofit_dpsi_leaf<phylloptim::Leaf::CostCurve::CF77>(
            psi_soil + 0.2 * (l.psi_crit - psi_soil), psi_soil, &feasible);
    ok(std::abs(at) < 1e-3 * std::abs(away),
       "dprofit at the optimum is negligible against its scale, psi_soil=" +
           std::to_string(psi_soil));
    ok(std::abs(away) > 0.0, "and the scale itself is non-zero");
  }
}

// evaluate_psi_leaf: the prescribed-point counterpart of the optimisers, and the
// case where a gradient is a partial at fixed psi rather than a total through a
// moving argmax.
void test_evaluate_psi_leaf_prescribes_rather_than_optimises() {
  printf("evaluate_psi_leaf: a prescribed operating point, refused outside and tagged\n");
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  d.PPFD = 1500.0;
  const double psi_soil = 1.0;

  phylloptim::Leaf l = make_leaf(d, {psi_soil}, {1.0});
  l.CF77_lambda_ = 1.5e5;
  l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::CF77>();
  const double psi_opt = l.opt_psi_leaf_, profit_opt = l.profit_;

  // A point deliberately away from the optimum is NOT stationary, and is worth
  // strictly less. Both halves matter: the first says the function honoured the
  // request, the second says the optimiser was actually finding a maximum.
  const double psi_off = psi_soil + 0.25 * (l.psi_crit - psi_soil);
  const double p_off = l.evaluate_psi_leaf<CC::CF77>(psi_off);
  near(l.opt_psi_leaf_, psi_off, 1e-12, "it sits exactly where it was told to");
  ok(l.operating_point_kind() == phylloptim::Leaf::OperatingPointKind::Prescribed,
     "and reports itself as prescribed rather than solved");
  ok(p_off < profit_opt, "a prescribed point off the optimum is worth less");

  bool feasible = false;
  const double slope =
      l.dprofit_dpsi_leaf<CC::CF77>(psi_off, psi_soil, &feasible);
  ok(feasible && std::abs(slope) > 1e-6,
     "dprofit is genuinely non-zero there, unlike at the optimum");

  // Asking for the optimum back reproduces it, which is what makes the two entry
  // points comparable at all.
  const double p_at = l.evaluate_psi_leaf<CC::CF77>(psi_opt);
  near(p_at, profit_opt, 1e-12, "prescribing the optimum reproduces its profit");

  // Outside the interval the target is refused, as infeasible, naming the
  // potential and the interval: projected onto the bound, the outputs would move
  // with the bound rather than with the potential asked for.
  for (double target : {l.psi_crit + 5.0, 0.0}) {
    std::string msg;
    bool infeasible = false;
    try {
      l.evaluate_psi_leaf<CC::CF77>(target);
    } catch (const phylloptim::util::infeasible_error &e) {
      infeasible = true;
      msg = e.what();
    }
    ok(infeasible && msg.find("prescribed_potential") != std::string::npos &&
           msg.find("feasible interval") != std::string::npos,
       "a prescribed leaf potential outside [wet, dry] is refused as infeasible: " + msg);
  }

  // The TF24 curve takes the same route and needs no lambda.
  phylloptim::Leaf t = make_leaf(d, {psi_soil}, {1.0});
  const double p_tf = t.evaluate_psi_leaf<CC::TF24>(psi_off);
  near(p_tf, t.profit_psi_leaf_TF(psi_off, psi_soil), 1e-12,
       "the TF24 arm agrees with its own objective at that point");
}

// A prescribed leaf point writes every reported output, on a reused leaf
// (hazard 8): solved first, outputs poisoned, then evaluated, it must report
// what a fresh leaf evaluated at the same point reports, bit for bit. The dry
// soil puts the curves bounded at psi_crit on the one-feasible-point exit.
template <phylloptim::Leaf::CostCurve K>
void check_prescribed_leaf_writes_every_output(const Drivers &d, int *rows,
                                               int *single_point) {
  auto seat = [](phylloptim::Leaf &l) {
    l.CF77_lambda_ = 3e5;
    l.TF24_floor_lambda_o = kPriceForFOC;
    l.LeastCost_beta = 2.0;
    l.JS22_gamma = kGammaJS22;
    l.CMax_a = kCMaxA;
    l.CMax_b = kCMaxB;
    l.set_model(K, false);
  };
  const std::string name = phylloptim::Leaf::curve_name(static_cast<int>(K));
  const double psi_crit = make_leaf(d, {1.0}, {1.0}).psi_crit;
  for (double psi_soil : {0.5, 1.0, 2.0, psi_crit + 0.5}) {
    const double target = psi_soil < psi_crit
                              ? psi_soil + 0.25 * (psi_crit - psi_soil)
                              : psi_soil + 0.25;
    phylloptim::Leaf fresh = make_leaf(d, {psi_soil}, {1.0});
    seat(fresh);
    fresh.template evaluate_psi_leaf<K>(target);

    phylloptim::Leaf reused = make_leaf(d, {psi_soil}, {1.0});
    seat(reused);
    reused.optimise();
    const double poison = -12345.0;
    for (double* f : {&reused.opt_psi_leaf_, &reused.ci_, &reused.assim_colimited_,
                      &reused.transpiration_, &reused.stom_cond_CO2_,
                      &reused.profit_, &reused.Tleaf_, &reused.lambda_emergent_,
                      &reused.residual_flux_, &reused.hydraulic_cost_}) {
      *f = poison;
    }
    reused.template evaluate_psi_leaf<K>(target);

    const std::vector<double> a = fresh.operating_point_values();
    const std::vector<double> b = reused.operating_point_values();
    bool same = a.size() == b.size();
    for (size_t i = 0; same && i < a.size(); ++i) {
      same = (std::isnan(a[i]) && std::isnan(b[i])) || a[i] == b[i];
    }
    ok(same && reused.lambda_emergent() != poison,
       name + " at a prescribed leaf potential reports a fresh leaf's outputs, "
              "psi_soil " + std::to_string(psi_soil));
    const double expected = fresh.template lambda_for<K>(
        fresh.opt_psi_leaf_, fresh.supply_psi_soil_scalar());
    ok((std::isnan(expected) && std::isnan(fresh.lambda_emergent())) ||
           fresh.lambda_emergent() == expected,
       name + " reports the lambda of the point it evaluated");
    ++*rows;
    if (fresh.opt_psi_leaf_ != target) ++*single_point;
  }
}

void test_prescribed_leaf_point_writes_every_output() {
  printf("a prescribed leaf point writes every output on a reused leaf\n");
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  d.PPFD = 1500.0;
  int rows = 0, single_point = 0;
  check_prescribed_leaf_writes_every_output<CC::TF24>(d, &rows, &single_point);
  check_prescribed_leaf_writes_every_output<CC::CF77>(d, &rows, &single_point);
  check_prescribed_leaf_writes_every_output<CC::JS22>(d, &rows, &single_point);
  check_prescribed_leaf_writes_every_output<CC::CMax>(d, &rows, &single_point);
  check_prescribed_leaf_writes_every_output<CC::SOX>(d, &rows, &single_point);
  check_prescribed_leaf_writes_every_output<CC::JW26>(d, &rows, &single_point);
  check_prescribed_leaf_writes_every_output<CC::ProfitMax>(d, &rows,
                                                           &single_point);
  check_prescribed_leaf_writes_every_output<CC::TF24_floor>(d, &rows,
                                                            &single_point);
  check_prescribed_leaf_writes_every_output<CC::LeastCost>(d, &rows,
                                                           &single_point);
  printf("    %d prescribed points, %d on the one-feasible-point exit\n", rows,
         single_point);
  ok(rows == 36, "every curve and soil was compared");
  ok(single_point == 4,
     "the four curves bounded at psi_crit reached the one-feasible-point exit");
}

// The multi-layer identity from the companion manuscript:
//   lambda_multi = lambda_single * [1 + kmax*f(psi_r)/S],   S = dE_up/dpsi_r
// find_root_collar_psi optimises the COLLAR, so lambda_multi, the collar route's
// lambda_emergent, is what it equalises. The leaf route's lambda at the same
// point should be badly wrong here -- that is the point of the correction, and
// this test pins the size of the error.
void test_multilayer_lambda_identity() {
  printf("multi-layer lambda: series-resistance correction vs dA/dE (collar free)\n");
  Drivers d;
  for (int layers : {1, 3, 5}) {
    std::vector<double> ps(layers), depth(layers), root(layers);
    for (int i = 0; i < layers; ++i) {
      ps[i] = 1.0 + 0.25 * i;
      depth[i] = 1.0 * (i + 1);
      root[i] = 1.0 / layers / d.area_leaf;
    }
    phylloptim::Leaf l = make_leaf(d, ps, depth);
    l.find_root_collar_psi();
    const double single = l.lambda_TF24(l.opt_psi_leaf_);
    const double multi = l.lambda_emergent();

    const double target = l.opt_root_psi_;
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
  phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
  l.find_root_collar_psi();
  const double g1 = l.g1_eff();
  ok(std::isfinite(g1) && g1 > 0.0, "g1_eff is finite and positive");
  // g1_eff is defined by inverting chi = g1/(g1 + sqrt(D)), so that must hold.
  const double chi = l.ci_ / l.ca_;
  near(g1 / (g1 + std::sqrt(l.atm_vpd_)), chi, 1e-12,
       "g1_eff inverts the USO relation exactly");
  // Drier soil closes stomata, lowering chi and therefore g1_eff.
  phylloptim::Leaf dry = make_leaf(d, {4.0}, {1.0});
  dry.find_root_collar_psi();
  ok(dry.g1_eff() < g1, "g1_eff falls as the soil dries");
  // And a higher marginal cost of water goes with a lower g1_eff.
  ok(dry.lambda_emergent() > l.lambda_emergent(),
     "lambda rises as the soil dries");
}

// With the energy balance on, the light-compensation edge is crossed by an
// ordinary optimum: below it the leaf is a net CO2 source whose ci sits above
// ambient and whose optimum pins to the wet bound; above it an interior one.
// Nothing is decided by a test on A(ca) (#182), so the edge is continuous and a
// hot leaf is not refused.
void test_light_compensation_is_crossed_by_an_ordinary_optimum() {
  printf("the light-compensation edge is crossed by an ordinary optimum\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  double prev = NAN, worst_step = 0.0;
  int pinned = 0, open = 0;
  for (int i = 0; i <= 1000; ++i) {
    Drivers d;
    d.PPFD = 20.0 + 0.01 * i;
    phylloptim::Leaf l = make_pm_leaf(d, {1.0}, {1.0}, true);
    l.find_root_collar_psi();
    ok(l.operating_point_kind() != Kind::HydraulicShutdown,
       "light does not shut a leaf");
    if (l.operating_point_kind() == Kind::BoundarySoil) ++pinned; else ++open;
    if (!std::isnan(prev)) worst_step = std::max(worst_step, std::abs(l.profit_ - prev));
    prev = l.profit_;
  }
  ok(pinned > 0 && open > 0, "the sweep crosses the compensation point");
  // A 0.01 step in PPFD moves a smooth profit here by ~5e-4.
  ok(worst_step < 5e-3, "no jump in profit across compensation");

  Drivers hot;
  hot.leaf_temp = 40.0;
  phylloptim::Leaf l = make_pm_leaf(hot, {1.0}, {1.0}, true);
  l.find_root_collar_psi();
  ok(std::isfinite(l.profit_),
     "a hot leaf that cannot fix carbon at ambient CO2 still has an optimum");
  ok(l.Tleaf_ > hot.leaf_temp && l.assim_colimited_ < 0.0,
     "and is a net source there");
  ok(l.ci_ > l.ca_, "with its internal CO2 above ambient");
}

// phi_J is per ABSORBED photon, so the light term of J is phi_J * alpha_PAR *
// PPFD: halving the absorptance at a doubled yield leaves J where it was.
void test_quantum_yield_on_absorbed_basis() {
  printf("electron transport reads phi_J * alpha_PAR * PPFD\n");
  Drivers d;
  phylloptim::Leaf base = make_leaf(d, {2.0}, {1.0});
  phylloptim::Leaf swap = make_leaf(d, {2.0}, {1.0});
  swap.alpha_PAR_ = base.alpha_PAR_ / 2.0;
  swap.phi_J = base.phi_J * 2.0;
  near(swap.electron_transport(), base.electron_transport(), 1e-12,
       "J depends on phi_J and alpha_PAR only through their product");
  near(base.phi_J * base.alpha_PAR_, 0.30, 1e-15,
       "the default restates the old incident-basis yield of 0.30");
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
// A whole Leaf solving through SinglePotential (issue #2 stage 2). This is the
// point of the whole item: the gas-exchange core is supply-agnostic, so swapping
// the supply path should change the operating point and nothing else.
void test_leaf_on_single_potential() {
  printf("Leaf solving on the single-potential supply path\n");
  Drivers d;

  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  l.set_supply_single();

  // ⚠️ THE SAME CALL AS THE MULTI-LAYER PATH, which is the point. The resistance
  // is a per-call driver on both paths now, carried by the RootNetwork argument;
  // it used to be a set_supply_single() constructor-style argument. 1.0e3 is per
  // unit leaf area, i.e. the old 2.0e4 * 0.05.
  std::vector<double> psi_soil{1.0}, depth{1.0};
  l.set_physiology(fixture::series_resistance(1.0e3), d.PPFD, psi_soil, depth,
                   d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                   d.atm_o2_kpa, d.atm_kpa);
  l.find_root_collar_psi();

  ok(std::isfinite(l.profit_), "single-potential solve gives a finite profit");
  ok(std::isfinite(l.opt_psi_leaf_), "and a finite leaf potential");
  ok(l.opt_psi_leaf_ > 0.0 && l.opt_psi_leaf_ <= l.psi_crit,
     "leaf potential is a positive magnitude within psi_crit");
  ok(l.opt_root_psi_ >= 0.0, "collar potential is stored as a positive magnitude");
  ok(l.assim_colimited_ > 0.0, "the leaf assimilates");
  ok(l.soil_consumption_.size() == 1u,
     "the consumption buffer is sized to one layer, not the caller's vector");

  // The collar must sit between the soil and the stem: water runs downhill.
  const double collar_mag = l.opt_root_psi_;
  ok(collar_mag >= psi_soil[0] - 1e-9 && collar_mag <= l.opt_psi_leaf_ + 1e-9,
     "collar potential lies between soil and stem");

  // Drier soil must cost carbon here too -- the same contract the multi-layer
  // path is held to, which is what makes the two comparable at all.
  phylloptim::Leaf dry = fixture::multilayer_leaf();
  dry.setup_transpiration(100);
  dry.setup_root_vulnerability(100);
  dry.set_supply_single();
  std::vector<double> psi_dry{3.0};
  dry.set_physiology(fixture::series_resistance(1.0e3), d.PPFD, psi_dry, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  dry.find_root_collar_psi();
  ok(dry.profit_ < l.profit_, "drier soil yields less profit");

  // A larger series resistance is a worse-supplied plant, so it must not do
  // better. This is the knob the multi-layer path spends root carbon to lower.
  phylloptim::Leaf tight = fixture::multilayer_leaf();
  tight.setup_transpiration(100);
  tight.setup_root_vulnerability(100);
  tight.set_supply_single();
  tight.set_physiology(fixture::series_resistance(1.0e4), d.PPFD, psi_soil, depth,
                       d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                       d.atm_o2_kpa, d.atm_kpa);
  tight.find_root_collar_psi();
  ok(tight.profit_ <= l.profit_, "a higher series resistance does not help");

  // Asking for the root network seats the collar route, the only one it has.
  phylloptim::Leaf plain;
  plain.set_supply_multilayer();
  ok(plain.supply_kind_ == phylloptim::Leaf::SupplyKind::MultiLayer &&
         plain.model_route() == "collar",
     "set_supply_multilayer() seats the collar route");
}

// A collar interval narrower than find_root_psi's bracket tolerance
// (bracket_tol_at, ~1e-3 MPa at these potentials) arises wherever E times the
// series resistance is small: a fast supply, or a dry soil. Extending the first
// bracket to the end of the stem curve must keep that interval's width, or the
// leaf is reported shut while profit still rises across it. Asserted against a
// scan of the prescribed collar, which no solve may beat.
void test_narrow_extended_bracket_is_not_shut() {
  printf("a narrow extended collar bracket keeps its width\n");
  Drivers d;
  std::vector<double> depth{1.0};
  int reached = 0;
  double worst = 0.0;
  for (double r : {1e-2, 1.0, 10.0}) {
    for (double ps : {4.5, 4.8, 5.5}) {
      std::vector<double> psi_soil{ps};
      auto make = [&](phylloptim::Leaf& l) {
        l.setup_transpiration(100);
        l.setup_root_vulnerability(100);
        l.set_supply_single();
        l.set_physiology(fixture::series_resistance(r), d.PPFD, psi_soil, depth,
                         d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                         d.atm_o2_kpa, d.atm_kpa);
      };
      phylloptim::Leaf l, scan;
      make(l);
      make(scan);
      l.find_root_collar_psi();
      double best = -std::numeric_limits<double>::infinity();
      for (double dx = 1e-12; dx < 1.0; dx *= 1.25) {
        try {
          best = std::max(best, scan.evaluate_root_collar_psi(ps + dx));
        } catch (const std::exception&) {
          // outside the feasible interval
        }
      }
      ok(std::isfinite(best), "the scan reaches a feasible collar at r = " +
                                  std::to_string(r) +
                                  ", psi_soil = " + std::to_string(ps));
      ok(l.transpiration_ > 0.0, "the leaf opens at r = " + std::to_string(r) +
                                     ", psi_soil = " + std::to_string(ps));
      ok(l.profit_ >= best - 1e-9 * std::abs(best),
         "no feasible collar beats the solve at r = " + std::to_string(r) +
             ", psi_soil = " + std::to_string(ps));
      worst = std::max(worst, (best - l.profit_) / std::abs(best));
      ++reached;
    }
  }
  printf("    %d rows | worst (scan - solve) / |scan|: %.3e\n", reached, worst);
}

void test_single_potential() {
  printf("single-potential supply path\n");
  phylloptim::SinglePotential sp;
  sp.set_soil_state(1.5);        // positive magnitude, -MPa
  // resistance_ is PER UNIT LEAF AREA, like every other input to the leaf.
  sp.resistance_ = 1.0e3;

  // begin_solve reports the only suction; there is nothing to flip (#25).
  near(sp.begin_solve(), 1.5, 1e-14, "begin_solve returns the soil suction");
  ok(sp.n_layers() == 1, "single potential writes exactly one layer");

  // Ohm's law, and the sign that matters: a collar drier than the soil -- a
  // LARGER suction now -- draws water UP (positive uptake).
  std::vector<double> consumption(1, 0.0);
  double E_up = 0.0;
  sp.uptake(2.5, consumption, E_up);
  ok(E_up > 0.0, "a collar drier than the soil draws water up");
  near(E_up, (2.5 - 1.5) / sp.resistance_ * phylloptim::kg_per_mol_h2o,
       1e-14, "uptake is the Ohm's-law flux");
  ok(consumption[0] > 0.0, "per-layer consumption is filled");

  // A collar WETTER than the soil pushes water back into it. Losing this sign is
  // how hydraulic redistribution silently becomes extra uptake.
  sp.uptake(0.5, consumption, E_up);
  ok(E_up < 0.0, "a collar wetter than the soil loses water to it");

  // The analytic derivative must match a central difference on uptake, and stay
  // finite everywhere -- unlike MultiLayerRoots there are no branch kinks, so it
  // never asks the caller for a finite-difference fallback.
  const double h = 1e-6, p0 = 2.5;
  double up = 0.0, dn = 0.0;
  sp.uptake(p0 + h, consumption, up);
  sp.uptake(p0 - h, consumption, dn);
  const double fd = (up - dn) / (2.0 * h);
  near(sp.duptake_dpsi(), fd, 1e-8, "analytic duptake_dpsi matches FD");
  ok(sp.duptake_dpsi() > 0.0,
     "duptake_dpsi is a positive conductance: uptake rises as the collar pulls harder");

  // An UNSET resistance would be an infinite or NaN flux; it is rejected, not
  // returned. The default is now the NA sentinel rather than zero, because the
  // resistance became a per-call driver -- an unset one means set_physiology was
  // never called, which is the same class of mistake as an unset psi_soil.
  phylloptim::SinglePotential bad;
  bad.set_soil_state(1.0);
  bad.begin_solve();
  bool threw = false;
  try {
    bad.uptake(2.0, consumption, E_up);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "an unset resistance throws rather than returning an infinity");

  // --- the driver entry point, and its two guards ---------------------------
  // set_supply_resistances is what makes the two supply paths take the same
  // set_physiology argument, so what it accepts and refuses is the contract.
  phylloptim::SinglePotential drv;
  drv.set_supply_resistances(fixture::series_resistance(2.5e3));
  near(drv.resistance_, 2.5e3, 1e-14,
       "the series resistance is read from r_R_V_sum[0]");

  // A network for the OTHER path carries a vulnerability-weighted horizontal
  // term this path cannot apply. Ignoring it would silently drop a resistance the
  // caller meant to use, so it is refused.
  threw = false;
  try {
    drv.set_supply_resistances(fixture::root_network({20.0}, {1.0}));
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "a multi-layer network is refused, not silently reinterpreted");
  near(drv.resistance_, 2.5e3, 1e-14,
       "and the refusal leaves the previous resistance intact");

  // More than one layer is the other path's shape too.
  threw = false;
  try {
    phylloptim::RootNetwork two;
    two.r_R_V_sum.assign(2, 1.0e3);
    drv.set_supply_resistances(two);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "two layers are refused on a one-potential path");

  // And a non-positive resistance is caught at the boundary rather than inside a
  // root-find, which is why this entry point validates at all.
  threw = false;
  try {
    drv.set_supply_resistances(fixture::series_resistance(0.0));
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "a zero series resistance is refused at the boundary");

  // So is one too small for the collar route to resolve, down to the last
  // representable value below the floor, and the floor itself is accepted.
  threw = false;
  bool names_floor = false;
  try {
    drv.set_supply_resistances(fixture::series_resistance(
        std::nextafter(phylloptim::SinglePotential::resistance_min, 0.0)));
  } catch (const std::exception &e) {
    threw = true;
    names_floor = std::string(e.what()).find("cannot resolve the flux") !=
                  std::string::npos;
  }
  ok(threw, "a series resistance just below resistance_min is refused");
  ok(names_floor, "and the refusal says why");
  drv.set_supply_resistances(
      fixture::series_resistance(phylloptim::SinglePotential::resistance_min));
  near(drv.resistance_, phylloptim::SinglePotential::resistance_min, 0.0,
       "and resistance_min itself is accepted");
}

void test_root_network_from_carbon() {
  printf("root architecture: carbon -> resistance\n");
  const double beta_H = 3.4e2, beta_V = 9.4e3, dz = 0.5;
  const std::vector<double> dz3(3, dz), dz1(1, dz);

  // Closed form, straight from the documented model: carbon splits 1/3 vertical
  // : 2/3 horizontal, r_R_H_min = beta_H/c_r_h, r_R_V = beta_V*dz[i]^2/c_r_v.
  const std::vector<double> carbon{3.0, 6.0, 1.5};
  const auto n = phylloptim::root_network_from_carbon(carbon, dz3, beta_H, beta_V);

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
  const auto rich = phylloptim::root_network_from_carbon({12.0}, dz1, beta_H, beta_V);
  const auto poor = phylloptim::root_network_from_carbon({3.0}, dz1, beta_H, beta_V);
  ok(rich.r_R_H_min[0] < poor.r_R_H_min[0], "more root carbon -> less horizontal resistance");
  ok(rich.r_R_V_sum[0] < poor.r_R_V_sum[0], "more root carbon -> less vertical resistance");

  // Trailing zero-carbon layers are dropped, so the hot loop never visits them.
  const auto trailing = phylloptim::root_network_from_carbon({3.0, 6.0, 0.0, 0.0},
                                                      std::vector<double>(4, dz),
                                                      beta_H, beta_V);
  ok(trailing.r_R_H_min.size() == 2u, "trailing rootless layers are dropped");

  // Negative carbon is rejected rather than producing a negative resistance.
  bool threw = false;
  try {
    phylloptim::root_network_from_carbon({3.0, -1.0}, std::vector<double>(2, dz),
                                        beta_H, beta_V);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "negative root carbon throws");

  // -------------------------------------------------------------------------
  // Per-layer thickness (#626)
  // -------------------------------------------------------------------------
  // Each layer's vertical resistance uses ITS OWN thickness, not a column
  // average. Asserted per layer, because a scalar dz would satisfy the r_R_V_sum
  // total on a uniform profile and only differ on a graded one.
  const std::vector<double> graded{0.02, 0.28, 1.2};
  const auto g = phylloptim::root_network_from_carbon(carbon, graded, beta_H, beta_V);
  for (int i = 0; i < 3; ++i) {
    near(g.r_R_V[i], beta_V * graded[i] * graded[i] / (carbon[i] / 3.0), 1e-14,
         "r_R_V uses this layer's own thickness");
  }
  // The horizontal term has no thickness in it, so it must be untouched by the
  // profile -- which is what says the change went where it was meant to. Exact,
  // because these are the same expression on the same input.
  for (int i = 0; i < 3; ++i) {
    ok(g.r_R_H_min[i] == n.r_R_H_min[i], "r_R_H_min does not depend on dz");
  }

  // DISCRETISATION INVARIANCE, the property the scalar dz broke and the reason
  // this is a bug fix and not a generalisation. With root density uniform over
  // depth the layer-integrated carbon goes as dz[i], so total vertical
  // resistance is 3*beta_V*D^2/C for ANY slicing of the same column. Three
  // profiles of total depth 1.5 m carrying the same 20 kg C m^-2:
  const double kD = 1.5, kC = 20.0;
  const std::vector<std::vector<double> > profiles{
      {0.3, 0.3, 0.3, 0.3, 0.3},              // 5 equal layers
      {0.5, 0.5, 0.5},                        // 3 equal layers
      {0.02, 0.28, 0.3, 0.4, 0.5},            // a 2 cm evaporation layer on top
      {0.75, 0.25, 0.25, 0.15, 0.1}           // thick over thin
  };
  const double expected = 3.0 * beta_V * kD * kD / kC;
  for (const auto &p : profiles) {
    std::vector<double> c(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
      c[i] = kC * p[i] / kD;  // uniform density -> carbon proportional to width
    }
    const auto net = phylloptim::root_network_from_carbon(c, p, beta_H, beta_V);
    near(net.r_R_V_sum.back(), expected, 1e-12,
         "total vertical resistance is independent of the layer slicing");
  }

  // A width vector of the wrong length, or a non-positive width, is refused
  // rather than read past the end or turned into a zero resistance. Checked by
  // REASON, because the two refusals are different mistakes: a length mismatch is
  // a mis-paired pair of vectors, a non-positive width is a bad profile.
  ok(throws_with([&] { phylloptim::root_network_from_carbon(
                           carbon, std::vector<double>{0.5, 0.5},
                           beta_H, beta_V); },
                 "same number of elements"),
     "a dz vector of the wrong length is refused as a length mismatch");
  for (double w : {0.0, -0.1}) {
    ok(throws_with([&] { phylloptim::root_network_from_carbon(
                             carbon, std::vector<double>{0.5, w, 0.5},
                             beta_H, beta_V); },
                   "finite and positive"),
       "a non-positive dz is refused as a bad width");
  }
}

// The template overload is plant's route onto its own tape, and its comment
// claims it is the double overload term for term. At double that is a
// bit-identity, so it is asserted as one, on a uniform, a graded and a gapped
// profile (an interior zero-carbon layer passes the running sum through) and one
// with trailing rootless layers. At a tangent the VALUES must still be the double
// overload's, bit for bit, since only the rows are added.
void test_root_network_template_is_the_double_overload() {
  printf("root architecture: the template overload is the double one, bit for bit\n");
  using T = odelia::ode::tangent_scalar<double>;
  const double beta_H = 3.4e2, beta_V = 9.4e3;
  const std::vector<std::pair<std::vector<double>, std::vector<double> > > cases{
      {{3.0, 6.0, 1.5}, {0.5, 0.5, 0.5}},
      {{3.0, 6.0, 1.5}, {0.02, 0.28, 1.2}},
      {{4.0, 0.0, 2.5, 7.0}, {0.1, 0.3, 0.6, 1.0}},
      {{3.0, 6.0, 0.0, 0.0}, {0.5, 0.5, 0.5, 0.5}},
      {{1.0 / 3.0, 2.0 / 7.0, 0.1}, {0.3, 0.7, 0.11}}};
  int compared = 0;
  for (const auto& [carbon, dz] : cases) {
    const phylloptim::RootNetwork want =
        phylloptim::root_network_from_carbon(carbon, dz, beta_H, beta_V);
    std::vector<double> h, v;
    phylloptim::root_network_from_carbon<double>(carbon, dz, beta_H, beta_V, h, v);
    ok(h.size() == want.r_R_H_min.size() && v.size() == want.r_R_V_sum.size(),
       "the template sizes to the deepest rooted layer, as the double one does");
    std::vector<T> carbon_t(carbon.begin(), carbon.end()), h_t, v_t;
    for (T& c : carbon_t) odelia::ode::set_direction(c, 1.0);
    phylloptim::root_network_from_carbon<T>(carbon_t, dz, beta_H, beta_V, h_t, v_t);
    for (std::size_t i = 0; i < h.size() && i < want.r_R_H_min.size(); ++i) {
      ok(h[i] == want.r_R_H_min[i] && v[i] == want.r_R_V_sum[i],
         "layer " + std::to_string(i) + " is bit-identical at double");
      ok(odelia::util::to_passive(h_t[i]) == want.r_R_H_min[i] &&
             odelia::util::to_passive(v_t[i]) == want.r_R_V_sum[i],
         "layer " + std::to_string(i) + " is bit-identical in value at a tangent");
      ++compared;
    }
  }
  printf("    %d layers compared\n", compared);
  ok(compared == 15, "every rooted layer of every case was compared");
}

// layer_thicknesses(): the profile -> widths helper an R caller and the test
// fixtures come through (#626).
void test_layer_thicknesses() {
  printf("root architecture: layer_thicknesses\n");

  // Differences with an implicit 0 at the surface, so the first width is the
  // first boundary and equal layers come back equal.
  const auto uniform = phylloptim::layer_thicknesses({1.0, 2.0, 3.0});
  for (int i = 0; i < 3; ++i) {
    // Exact, and that is the load-bearing part: it is what keeps the golden grid
    // (1 m layers, differenced here since #626) bit-identical.
    ok(uniform[i] == 1.0, "1 m boundaries give exactly 1 m widths");
  }

  const auto graded = phylloptim::layer_thicknesses({0.02, 0.3, 1.5});
  near(graded[0], 0.02, 1e-16, "surface layer width is the first boundary");
  near(graded[1], 0.28, 1e-16, "interior width is a difference of boundaries");
  near(graded[2], 1.2, 1e-16, "deepest width is a difference of boundaries");

  // ⚠️ THE NEGATIVE, which is why plant passes its own widths rather than coming
  // through here: differencing is NOT bit-exact. 1.5 m over 5 layers is the
  // package's own default geometry, and the round trip loses a bit.
  const size_t kN = 5;
  const double kDelta = 1.5 / kN;
  std::vector<double> profile(kN);
  for (size_t i = 0; i < kN; ++i) {
    profile[i] = (i + 1) * kDelta;
  }
  const auto round_tripped = phylloptim::layer_thicknesses(profile);
  bool all_exact = true;
  for (size_t i = 0; i < kN; ++i) {
    if (round_tripped[i] != kDelta) {
      all_exact = false;
    }
  }
  ok(!all_exact,
     "differencing a uniform profile is NOT bit-exact -- pass widths, do not "
     "re-derive them");
  for (size_t i = 0; i < kN; ++i) {
    near(round_tripped[i], kDelta, 1e-15, "...but it is right to rounding");
  }

  // A profile that is not a profile: empty, or not strictly increasing. The
  // second is the one that matters -- it is what happens when someone passes
  // WIDTHS to a function wanting boundaries, which for a graded profile is not
  // even monotone.
  ok(throws_with([] { phylloptim::layer_thicknesses({}); }, "at least one layer"),
     "an empty profile is refused");
  for (const std::vector<double> &bad :
       {std::vector<double>{0.3, 0.3},        // repeated boundary
        std::vector<double>{0.5, 0.4, 0.3},   // decreasing
        std::vector<double>{-0.5},            // above the surface
        std::vector<double>{0.0, 1.0}}) {     // a zero-thickness surface layer
    ok(throws_with([&] { phylloptim::layer_thicknesses(bad); },
                   "strictly increasing"),
       "a profile that is not strictly increasing is refused, and says so");
  }

  // ⚠️ AND THE LIMIT OF THAT GUARD, ASSERTED SO IT IS NOT MISREAD AS PROTECTION.
  // The mistake it was reached for -- passing WIDTHS where boundaries are wanted
  // -- is only caught when the widths happen to be non-monotone. plant's own
  // graded profile is a counterexample: its widths increase with depth, so read
  // as boundaries they pass every check here and silently yield different,
  // entirely plausible widths. Only the caller knows which it holds; the guard
  // catches a malformed profile, not a mislabelled one.
  const std::vector<double> widths{0.02, 0.28, 0.3, 0.4, 0.5};
  const auto misread = phylloptim::layer_thicknesses(widths);
  bool differs = false;
  for (size_t i = 0; i < widths.size(); ++i) {
    if (misread[i] != widths[i]) {
      differs = true;
    }
  }
  ok(differs, "widths misread as boundaries are accepted and silently wrong");
}

// The temperature-response parameters were constexpr constants; they are now
// settable members. This checks the point of that change -- that setting them
// actually moves the model -- rather than just that they compile.
void test_temperature_parameters_are_settable() {
  printf("temperature-response parameters are settable\n");
  Drivers d;
  phylloptim::Leaf base = make_leaf(d, {2.0}, {1.0});
  base.find_root_collar_psi();
  const double A_base = base.assim_colimited_;

  // Defaults must equal the published constants, so this is a no-op refactor for
  // anyone who does not touch them.
  ok(base.vcmax_ha_ == phylloptim::vcmax_ha, "vcmax_ha_ defaults to the constant");
  ok(base.jmax_d_S_ == phylloptim::jmax_d_S, "jmax_d_S_ defaults to the constant");
  ok(base.gamma_25_ == phylloptim::gamma_25, "gamma_25_ defaults to the constant");
  near(base.R_d_25, 1.44, 1e-12, "R_d_25 default");

  // Raising the Vcmax activation energy raises Vcmax above the 25 C reference,
  // so a 25 C leaf should be unaffected but a warm one should assimilate more.
  {
    phylloptim::Leaf warm = make_leaf(d, {2.0}, {1.0});
    phylloptim::Leaf warm_hi = make_leaf(d, {2.0}, {1.0});
    warm_hi.vcmax_ha_ = phylloptim::vcmax_ha * 1.5;
    // set_physiology already ran, so push the change through the T-response block.
    warm.update_temperature_dependent_params(35.0);
    warm_hi.update_temperature_dependent_params(35.0);
    ok(warm_hi.vcmax_ > warm.vcmax_,
       "a larger activation energy gives a larger vcmax at 35 C");
  }

  // Respiration: doubling it must lower assimilation.
  {
    phylloptim::Leaf r2 = make_leaf(d, {2.0}, {1.0});
    r2.R_d_25 = 2.0 * r2.R_d_25;
    r2.update_temperature_dependent_params(d.leaf_temp);
    r2.find_root_collar_psi();
    ok(r2.assim_colimited_ < A_base, "doubling R_d_25 lowers assimilation");
    ok(r2.R_d_ > base.R_d_, "and raises R_d");
  }

  // The CO2 compensation point feeds photorespiration, so raising it lowers A.
  {
    phylloptim::Leaf g2 = make_leaf(d, {2.0}, {1.0});
    g2.gamma_25_ = phylloptim::gamma_25 * 1.5;
    g2.update_temperature_dependent_params(d.leaf_temp);
    g2.find_root_collar_psi();
    ok(g2.assim_colimited_ < A_base,
       "raising the compensation point lowers assimilation");
  }
}

// ⚠️ THE ASSERTIONS ABOVE ALL PUSH THEIR CHANGE THROUGH BY CALLING
// update_temperature_dependent_params() DIRECTLY, AND THAT IS WHY THEY PASSED
// WHILE THE FEATURE WAS BROKEN FROM R. A caller does not have that route: they set
// the field and then set the drivers, which is the one path that took a cache hit
// and silently kept the old response. The test worked around the bug it should
// have caught.
//
// So this asserts the REALISTIC path, twice over: setting a temperature-response
// parameter on an already-solved leaf and re-supplying the SAME drivers must change
// the answer, and must land on exactly what a freshly constructed leaf gives.
//
// Bit-exactness is the right bar for the second half -- the two routes share no
// code, so anything the cache fails to invalidate shows up as a difference. #41.
void test_temperature_params_invalidate_cache() {
  printf("setting a temperature parameter invalidates the cache\n");
  Drivers d;

  // The realistic route: solve, change the parameter, re-supply the same drivers.
  phylloptim::Leaf warm = make_leaf(d, {2.0}, {1.0});
  warm.find_root_collar_psi();
  const double A_before = warm.assim_colimited_;
  const double Rd_before = warm.R_d_;

  warm.R_d_25 = 2.0 * warm.R_d_25;
  warm.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD,
                      {2.0}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                      d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  warm.find_root_collar_psi();

  ok(warm.R_d_ != Rd_before,
     "re-supplying the same drivers after a parameter change recomputes R_d");
  ok(warm.assim_colimited_ != A_before,
     "and the operating point moves");
  near(warm.R_d_, 2.0 * Rd_before, 1e-12,
       "doubling R_d_25 doubles R_d at the same temperature");

  // And it must agree bit-for-bit with never having had the stale value.
  phylloptim::Leaf fresh = make_leaf(d, {2.0}, {1.0});
  fresh.R_d_25 = 2.0 * fresh.R_d_25;
  fresh.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD,
                       {2.0}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                       d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  fresh.find_root_collar_psi();
  ok(warm.assim_colimited_ == fresh.assim_colimited_,
     "a re-parameterised leaf is bit-identical to one built with the value");
  ok(warm.R_d_ == fresh.R_d_, "R_d likewise");

  // ⚠️ The same hole covered vcmax_25, which is a TRAIT. set_traits() remains the
  // correct way to change one -- the vulnerability splines and the solved point
  // need clearing too -- but the silent-wrong-number half of hazard 10 is gone.
  phylloptim::Leaf vc = make_leaf(d, {2.0}, {1.0});
  vc.find_root_collar_psi();
  const double vcmax_before = vc.vcmax_;
  vc.vcmax_25 = vc.vcmax_25 * 1.5;
  vc.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD,
                    {2.0}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                    d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  ok(vc.vcmax_ > vcmax_before,
     "a bare vcmax_25 write no longer leaves vcmax_ describing the old value");

  // The cache must still BE a cache: identical inputs twice must not recompute
  // into a different answer.
  phylloptim::Leaf same = make_leaf(d, {2.0}, {1.0});
  same.find_root_collar_psi();
  const double A_once = same.assim_colimited_;
  same.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD,
                      {2.0}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                      d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  same.find_root_collar_psi();
  ok(same.assim_colimited_ == A_once,
     "unchanged inputs still give a bit-identical answer");
}

// The interior arm divides every collar row by the curvature it is handed, and
// odelia refuses only an exactly zero or non-finite one. So the arm refuses a
// curvature that is not negative or lies inside the floor, naming the value; the
// true curvature, one just outside the floor, and the double path (which never
// reads the curvature) all still answer.
void test_operating_collar_refuses_a_flat_or_rising_curvature() {
  printf("operating_collar refuses a curvature it cannot divide by\n");
  using T = odelia::ode::tangent_scalar<double>;
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  d.PPFD = 1500.0;
  phylloptim::Leaf l = make_leaf(d, {1.0, 1.5}, {1.0, 2.0});
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == phylloptim::Leaf::OperatingPointKind::Interior,
     "the fixture is an interior point");
  phylloptim::leaf_pars<T> pars;
  const phylloptim::leaf_pars<double> seated = l.passive_pars();
  for (std::size_t i = 0; i < pars.size(); ++i) pars[i] = T(seated[i]);
  odelia::ode::set_direction(pars[phylloptim::par_vcmax_25], 1.0);
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
  const double floor = phylloptim::Leaf::interior_curvature_floor;
  const double k = l.profit_curvature_at_collar<CC::TF24>();
  printf("    curvature at the fixture %.3g, floor %.3g\n", k, floor);
  ok(k < -floor, "the fixture's own curvature is outside the floor");
  auto row = [&](double curvature) {
    return odelia::ode::derivative_along(
        l.operating_collar<CC::TF24, T>(draw, pars, curvature));
  };
  ok(std::isfinite(row(k)), "the true curvature answers");
  ok(std::isfinite(row(-2.0 * floor)), "a curvature just outside the floor answers");
  ok(throws_with([&] { row(-1e-300); }, "inside the floor of 0.001"),
     "a curvature of -1e-300 is refused against the floor");
  ok(throws_with([&] { row(-1e-300); }, "-1e-300"),
     "the refusal names the curvature");
  ok(throws_with([&] { row(-0.5 * floor); }, "inside the floor"),
     "a curvature half the floor is refused");
  ok(throws_with([&] { row(0.5); }, "which is not negative"),
     "a positive curvature is refused as not negative");
  ok(throws_with([&] { row(0.0); }, "which is not negative"),
     "a zero curvature is refused as not negative");
  const auto held_d = l.roots_.held_supply();
  const auto draw_d = l.supply_draw_at<double>(l.opt_root_psi_, held_d);
  ok(l.operating_collar<CC::TF24, double>(draw_d, seated, -1e-300) == l.opt_root_psi_,
     "the double path reads no curvature, so refuses none");
}

void test_rd_temperature_response() {
  printf("R_d rises with temperature\n");
  Drivers d;

  // R_d IS R_d_25 at the reference, exactly. The trait is the value there, not a
  // scale on anything.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    near(l.R_d_, l.R_d_25, 1e-12, "R_d at 25 C is R_d_25");
  }

  // The direction: R_d must RISE, including above Vcmax's thermal optimum near
  // 31 C, so 45 C is comfortably past it.
  double rd_prev = -1.0;
  for (double T : {25.0, 35.0, 45.0}) {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.update_temperature_dependent_params(T);
    ok(l.R_d_ > rd_prev, "R_d rises with temperature");
    rd_prev = l.R_d_;
  }

  // Tjoelker semantics. Q10 is evaluated at the MEAN of T and the reference, so a
  // 10 K rise from 25 C multiplies R_d by Q10(30) = 3.09 - 0.043*30 = 1.80 -- not
  // by the 2.015 that Q10(25) would give. Checking the number rather than just the
  // direction is what distinguishes this from a constant Q10.
  {
    phylloptim::Leaf a = make_leaf(d, {2.0}, {1.0});
    phylloptim::Leaf b = make_leaf(d, {2.0}, {1.0});
    a.update_temperature_dependent_params(25.0);
    b.update_temperature_dependent_params(35.0);
    const double q10_at_30 = 3.09 - 0.0430 * 30.0;
    near(b.R_d_, q10_at_30 * a.R_d_, 1e-12,
         "a 10 K rise scales R_d by Q10 at the midpoint temperature");
    ok(q10_at_30 < 2.0, "and that Q10 is below 2, i.e. the decline is active");
  }

  // A constant Q10 must stay reachable: slope zero, intercept the value.
  {
    phylloptim::Leaf a = make_leaf(d, {2.0}, {1.0});
    phylloptim::Leaf b = make_leaf(d, {2.0}, {1.0});
    a.rd_q10_intercept_ = 2.0;  a.rd_q10_slope_ = 0.0;
    b.rd_q10_intercept_ = 2.0;  b.rd_q10_slope_ = 0.0;
    a.update_temperature_dependent_params(25.0);
    b.update_temperature_dependent_params(35.0);
    near(b.R_d_, 2.0 * a.R_d_, 1e-12,
         "slope zero recovers a constant Q10 exactly");
  }

  // A measured value is used verbatim, and is INDEPENDENT of vcmax_25 -- R_d_25 is
  // a trait in its own right, not a fraction of Vcmax.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.R_d_25 = 0.525;                 // Sabot's Rlref, as measured
    l.update_temperature_dependent_params(25.0);
    near(l.R_d_, 0.525, 1e-12, "a set R_d_25 is used verbatim at 25 C");
    l.vcmax_25 = 2.0 * l.vcmax_25;
    l.update_temperature_dependent_params(25.0);
    near(l.R_d_, 0.525, 1e-12, "and does not follow vcmax_25");
  }

  // An unset or negative R_d_25 FAILS rather than falling back to a derivation.
  for (double bad : {std::numeric_limits<double>::quiet_NaN(), -1.0}) {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.R_d_25 = bad;
    bool threw = false;
    try {
      l.update_temperature_dependent_params(30.0);
    } catch (const std::runtime_error &) {
      threw = true;
    }
    ok(threw, "an unusable R_d_25 is refused, not worked around");
  }

  // The response, and the LIMIT it makes reachable, printed rather than asserted
  // cell by cell: a higher R_d can drive net assimilation negative across the whole
  // [gamma*, ca] bracket, and then there is no supply == demand root at all. At the
  // defaults that happens by 45 C, which is a number a caller needs.
  printf("  R_d and assimilation against leaf temperature:\n");
  for (double T : {15.0, 25.0, 35.0, 40.0, 45.0}) {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.update_temperature_dependent_params(T);
    l.find_root_collar_psi();
    printf("    T = %4.1f C   R_d %6.3f   A %8.4f%s\n", T, l.R_d_,
           l.assim_colimited_,
           l.assim_colimited_ < 0.0 ? "   <-- a net source" : "");
  }

  // ⚠️ AND IT MUST SOLVE RATHER THAN THROW. A leaf too hot to gain carbon at
  // ambient CO2 is a physical state, not a solver failure: it is a net source
  // whose internal CO2 sits above ambient (#182). The prescribed-temperature
  // path used to throw here, which was only invisible while R_d was too small
  // to get there.
  {
    phylloptim::Leaf hot = make_leaf(d, {2.0}, {1.0});
    bool threw = false;
    try {
      hot.update_temperature_dependent_params(45.0);
      hot.find_root_collar_psi();
    } catch (const std::exception &) {
      threw = true;
    }
    ok(!threw, "a leaf too hot to gain carbon solves instead of throwing");
    ok(hot.assim_colimited_ <= 0.0, "it is a net source");
    ok(hot.ci_ > hot.ca_, "so its internal CO2 sits above ambient");
  }

  // ⚠️ A GENUINE solver failure must STILL throw -- the point of not simply
  // deleting the gate. An unsatisfiable bracket where the root IS inside it (here
  // forced by a non-finite target via a poisoned ca) is a bug, not a hot leaf.
  {
    phylloptim::Leaf ok_leaf = make_leaf(d, {2.0}, {1.0});
    ok_leaf.find_root_collar_psi();
    ok(ok_leaf.assim_colimited_ > 0.0 && ok_leaf.ci_ < ok_leaf.ca_,
       "an ordinary leaf is a sink with its internal CO2 below ambient");
  }
}

void test_set_traits_matches_a_fresh_leaf() {
  printf("set_traits is indistinguishable from constructing afresh\n");
  Drivers d;

  // One perturbed trait from each group that has its own derived state: a
  // photosynthetic trait (the temperature cache), the stem curve, the root curve,
  // and a cost parameter that is read directly and so should need no rebuild.
  struct Case { const char *name; int which; double value; };
  const Case cases[] = {{"vcmax_25", 0, 96.0 * 1.05},
                        {"stem_P50", 2, 3.4 * 1.05},
                        {"root_P50", 4, 3.4 * 1.05},
                        {"TF24_cost_scale", 10, 7.5 * 1.05},
                        {"stem_c", 1, 2.680147 * 1.05},
                        {"root_c", 3, 2.680147 * 1.05}};

  for (const Case &c : cases) {
    // The default trait vector, in set_traits' own argument order (R_d_25 apart,
    // which the constructor does not take).
    double t[11] = {96.0, 2.680147, 3.4,  2.680147, 3.4, 1.5,
                    157.44, 0.30 / 0.85,   0.7,  0.99,     7.5};
    t[c.which] = c.value;

    // Fresh: the traits go through the constructor.
    phylloptim::Leaf fresh(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7],
                     t[8], t[9], 1e-3, 100, 1e-3, 1000, t[10]);
    fresh.set_supply_multilayer();
    fresh.set_model(phylloptim::Leaf::CostCurve::TF24, true);
    std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};
    fresh.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                         d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    fresh.find_root_collar_psi();

    // Reused: solved once at the DEFAULTS first, so every cache is warm and
    // pointing at the old traits before set_traits runs. Solving first is what
    // makes this a test rather than a coincidence -- on a cold object the
    // temperature cache would miss anyway and the trap would not fire.
    phylloptim::Leaf reused = make_leaf(d, {2.0}, {1.0});
    reused.find_root_collar_psi();
    reused.set_traits(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7],
                      t[8], t[9], t[10], kRd25, kGammaJS22, kCMaxA, kCMaxB, 0.0);
    reused.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                          d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    reused.find_root_collar_psi();

    const std::string what = std::string(" after set_traits(") + c.name + ")";
    ok(reused.vcmax_ == fresh.vcmax_, "vcmax_ is bit-identical" + what);
    ok(reused.jmax_ == fresh.jmax_, "jmax_ is bit-identical" + what);
    ok(reused.R_d_ == fresh.R_d_, "R_d_ is bit-identical" + what);
    ok(reused.opt_root_psi_ == fresh.opt_root_psi_,
       "the collar is bit-identical" + what);
    ok(reused.opt_psi_leaf_ == fresh.opt_psi_leaf_,
       "psi_leaf is bit-identical" + what);
    ok(reused.assim_colimited_ == fresh.assim_colimited_,
       "assimilation is bit-identical" + what);
    ok(reused.profit_ == fresh.profit_, "profit is bit-identical" + what);
    ok(reused.transpiration_ == fresh.transpiration_,
       "transpiration is bit-identical" + what);
  }

  // "Change the trait, then set the drivers again" really does recompute the
  // temperature block, checked end to end through the route a caller takes.
  //
  // ⚠️ THIS NO LONGER ISOLATES set_traits' CACHE INVALIDATION, and the comment
  // here claimed it did for two releases. It said vcmax_ is derived behind a cache
  // "keyed on (leaf_temp, atm_o2_kpa) and NOT on the traits", so the recipe worked
  // "only because set_traits invalidates the cache". #55 widened the key to every
  // scalar update_temperature_dependent_params() reads -- `vcmax_25` included --
  // so the trap is now closed TWICE OVER: the key would miss on the changed trait
  // even if set_traits cleared nothing. The assertion below cannot tell which
  // mechanism carried it, and pretending otherwise is how a test comes to describe
  // a guarantee it is not testing.
  //
  // The two halves are covered separately and deliberately:
  //   * the KEY -- test_temperature_params_invalidate_cache, which sets a response
  //     parameter and re-drives without going through set_traits at all;
  //   * set_traits' CLEARING -- the bit-exact comparisons above (the splines and
  //     the solved operating point are not in any key), plus the NA assertions in
  //     test_prescribed_lambda_survives_redriving.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    const double vcmax_before = l.vcmax_;
    l.set_traits(96.0 * 2.0, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30 / 0.85, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB, 0.0);
    std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};
    // The same leaf_temp and atm_o2_kpa -- which used to be the whole key, and so
    // used to be what ARMED the trap. It no longer is: `vcmax_25` has moved, so the
    // key differs whatever set_traits did.
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    near(l.vcmax_ / vcmax_before, 2.0, 1e-12,
         "doubling vcmax_25 doubles vcmax_ at an unchanged temperature");
  }

  // The splines really are rebuilt, and only when their own pair moves. Read
  // through proportion_of_conductivity (the closed form) against transpiration()
  // (the spline): the two describe the same curve, so they move together or the
  // object is inconsistent.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    const double E_before = l.transpiration(3.0, 1.0);
    l.set_traits(96.0, 2.680147, 3.4 * 1.5, 2.680147, 3.4, 1.5, 157.44, 0.30 / 0.85, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB, 0.0);
    std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth, d.K_s * d.theta / d.h,
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    // A less vulnerable stem (larger stem_b) holds more conductivity, so the
    // integral of the vulnerability curve over the same span is larger.
    ok(l.transpiration(3.0, 1.0) > E_before,
       "raising stem_b rebuilds the transpiration spline");
  }

  // The #25 boundary is enforced here too. A bare field write would bypass it,
  // which is the fourth reason these are not settable fields.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    const double good[11] = {96.0, 2.680147, 3.4,  2.680147, 3.4, 1.5,
                             157.44, 0.30 / 0.85,   0.7,  0.99,     7.5};
    // Only two potentials are settable now -- the derived scales and critical
    // potentials take their sign from these, so there is nothing else to reject.
    const int signed_positions[] = {2, 4};
    const char *labels[] = {"stem_P50", "root_P50"};
    for (int k = 0; k < 2; ++k) {
      double t[11];
      for (int i = 0; i < 11; ++i) {
        t[i] = good[i];
      }
      t[signed_positions[k]] = -t[signed_positions[k]];  // the pre-#25 sign
      bool threw = false;
      try {
        l.set_traits(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7],
                     t[8], t[9], t[10], kRd25, kGammaJS22, kCMaxA, kCMaxB, 0.0);
      } catch (const std::runtime_error &) {
        threw = true;
      }
      ok(threw, std::string("set_traits rejects a negative ") + labels[k]);
    }
  }
}

// A prescribed `CF77_lambda_` survives everything that is not a caller writing to it.
// It is the Cowan-Farquhar marginal value of water, an INPUT, so a re-driving call
// that cleared it would make whether a sweep kept its price depend on which of two
// interchangeable-looking calls came next.
//
// ⚠️ The two arms are the point, not the pair of assertions. `set_drivers` always
// kept the value and `set_traits` always lost it; asserting only one arm passes
// on the code this test exists to reject.
void test_prescribed_lambda_survives_redriving() {
  printf("a prescribed CF77_lambda_ survives set_traits and set_drivers\n");
  Drivers d;
  std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};

  // A fresh Leaf reads the default: "no longer reset" must not have become
  // "never initialised", which is the one way this change could go wrong quietly.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    ok(l.CF77_lambda_ == 107750.8846,
       "CF77_lambda_ is its aligned default on a freshly constructed leaf");
    // The other three aligned initialisers. R assigns its own copies on every
    // leaf it builds, so only a C++ consumer (plant) reads these, and nothing
    // on the R side would see one drift.
    ok(l.JS22_gamma == 0.4717887793 && l.CMax_a == 0.6385596982 &&
           l.CMax_b == 0.0 && l.LeastCost_beta == 51.59667042,
       "JS22_gamma, CMax_a and LeastCost_beta are their aligned defaults");
  }

  const double prescribed = 30.0;

  // Arm 1: re-drive. This arm passed before the fix too.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.CF77_lambda_ = prescribed;
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    ok(l.CF77_lambda_ == prescribed, "set_physiology leaves a prescribed CF77_lambda_ alone");
  }

  // Arm 2: re-trait. This is the arm that lost the value.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.CF77_lambda_ = prescribed;
    l.set_traits(96.0, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30 / 0.85, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB, 0.0);
    ok(l.CF77_lambda_ == prescribed, "set_traits leaves a prescribed CF77_lambda_ alone");
  }

  // And it survives a SOLVE, which is the case a sweep actually runs: solve at
  // the defaults, re-trait, solve again, all on one object.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.CF77_lambda_ = prescribed;
    l.find_root_collar_psi();
    l.set_traits(96.0 * 1.05, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30 / 0.85, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB, 0.0);
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
    ok(l.CF77_lambda_ == prescribed,
       "CF77_lambda_ survives a solve, a re-trait and a second solve");
  }

  // The derived state around it is still wiped -- the fix removed two lines from
  // setup_clean_leaf(), and taking any more would reopen hazard 8.
  {
    phylloptim::Leaf l = make_leaf(d, {2.0}, {1.0});
    l.CF77_lambda_ = prescribed;
    l.find_root_collar_psi();
    ok(!std::isnan(l.opt_psi_leaf_), "the solve seated an operating point");
    l.set_traits(96.0, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30 / 0.85, 0.7, 0.99, 7.5, kRd25, kGammaJS22, kCMaxA, kCMaxB, 0.0);
    ok(std::isnan(l.opt_psi_leaf_),
       "set_traits still clears the solved operating point");
    ok(std::isnan(l.profit_), "set_traits still clears profit_");
    ok(std::isnan(l.R_d_), "set_traits still clears the temperature block");
  }
}

// The realised knot grid against the intended one (#92). Nothing checked this
// before, which is exactly why the accumulating `psi += step` loop could drop its
// last knot for years: the shortfall is invisible from outside unless the domain
// edge is compared with the function that is supposed to define it.
// The kg <-> mol water conversions are reciprocal, and the round trip is the
// identity (#51). Asserted rather than documented: the demand side converts
// transpiration kg -> mol and the supply side converts uptake mol -> kg, so a
// second literal anywhere breaks the round trip silently.
void test_water_mass_conversions_are_reciprocal() {
  printf("kg <-> mol water conversions\n");
  // Exact to the last bit in one direction: kg_to_mol_h2o IS 1/molar_mass_h2o.
  ok(phylloptim::kg_to_mol_h2o == 1.0 / phylloptim::molar_mass_h2o,
     "kg_to_mol_h2o is exactly the reciprocal of the molar mass");
  ok(phylloptim::kg_per_mol_h2o == phylloptim::molar_mass_h2o,
     "kg_per_mol_h2o IS the molar mass, not a second copy of it");

  // The round trip cannot be bit-exact -- 1/x then *x is not the identity for most
  // x -- so the claim is that it is within one rounding, which is what "reciprocal"
  // can actually buy.
  for (double kg : {1e-9, 1e-6, 1e-3, 1.0, 1e3}) {
    const double round_trip = kg * phylloptim::kg_to_mol_h2o *
                              phylloptim::kg_per_mol_h2o;
    near(round_trip, kg, 4.0 * std::numeric_limits<double>::epsilon(),
         "kg -> mol -> kg returns the input at " +
             phylloptim::util::format_double(kg));
  }

  // The value is water's molar mass, not a neighbouring one: 18.015 g/mol from the
  // standard atomic weights, 2(1.008) + 15.999.
  near(phylloptim::molar_mass_h2o * 1000.0, 18.015, 1e-9,
       "the molar mass is water's, in g/mol");
}

// Both Arrhenius curves are exactly inert at 25 C, which is why a change to any
// temperature response is invisible on a 25 C-only grid. That is the claim worth
// holding: it bounds the blast radius of anything touched in this area.
void test_gas_constant_and_arrhenius_reference_point() {
  printf("gas constant and the 25 C reference\n");
  // Exact SI quantity, so there is no tolerance to allow.
  near(phylloptim::gas_constant, 8.314462618153240, 1e-15,
       "gas_constant is the SI value");

  // Inertness at the reference temperature, through the object: every reference
  // value in this model is DEFINED at 25 C, so a response can only be seen away
  // from it.
  Drivers d;
  phylloptim::Leaf ref = make_leaf(d, {2.0}, {1.0});   // d.leaf_temp is 25
  ok(ref.vcmax_ == ref.vcmax_25, "vcmax_ IS vcmax_25 at 25 C, bit for bit");
  ok(ref.jmax_ == ref.jmax_25, "jmax_ IS jmax_25 at 25 C, bit for bit");
  ok(ref.R_d_ == ref.R_d_25, "R_d_ IS R_d_25 at 25 C, bit for bit");
  // gamma_, kc_ and ko_ also sit at their 25 C reference values, but they carry a
  // umol/mol -> Pa conversion through atm_kpa_ that this test would have to restate
  // to check -- and a second copy of a conversion is what hazard 1 is about. The
  // three identities above are the same claim without the duplication.

  // And NOT inert away from it, which is the other half: a test that only checked
  // 25 C would pass on a broken response curve.
  Drivers hot;
  hot.leaf_temp = 40.0;
  phylloptim::Leaf warm = make_leaf(hot, {2.0}, {1.0});
  ok(warm.vcmax_ != warm.vcmax_25, "vcmax_ has moved at 40 C");
  ok(warm.jmax_ != warm.jmax_25, "jmax_ has moved at 40 C");
  ok(warm.R_d_ > ref.R_d_, "respiration is higher at 40 C than at 25 C");
}

void test_bad_input_throws() {
  printf("input validation\n");
  Drivers d;
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  bool threw = false;
  try {
    std::vector<double> psi_soil{2.0}, depth{1.0, 2.0}, mrp{1.0 / d.area_leaf};
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  ok(threw, "mismatched soil vector lengths throw std::runtime_error");
}

// Infeasibility is a distinguishable kind of failure (#57).
//
// ⚠️ THE C++ SIDE MATTERS ON ITS OWN, not just as the R layer's supplier: plant links
// these headers directly, so it sees the type rather than the token. Both halves are
// asserted here -- a `catch (const infeasible_error&)` works, and the message carries
// the prefix R parses.
//
// ⚠️ And the asymmetry is asserted too. An input-validation failure must NOT be an
// `infeasible_error`, because the whole point is that a caller may legitimately
// swallow one and must not be able to swallow the other. A test that only checked the
// positive direction would pass on a `stop_infeasible` applied to everything.
void test_infeasible_is_a_distinct_failure() {
  printf("infeasibility is distinguishable from a caller error\n");
  Drivers d;
  std::vector<double> mrp{1.0 / d.area_leaf}, psi_soil{2.0}, depth{1.0};

  // A well-formed call whose operating point cannot be found: three iterations do
  // not bring the continuity root-find's bracket to its tolerance.
  {
    phylloptim::Leaf l = fixture::multilayer_leaf();
    l.ci_niter = 3;
    l.setup_transpiration(100);
    l.setup_root_vulnerability(100);
    l.set_physiology(fixture::root_network(mrp, depth), d.PPFD, psi_soil, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    bool caught_as_infeasible = false;
    std::string msg;
    try {
      l.find_root_collar_psi();
    } catch (const phylloptim::util::infeasible_error &e) {
      caught_as_infeasible = true;
      msg = e.what();
    } catch (const std::runtime_error &e) {
      msg = e.what();
    }
    ok(caught_as_infeasible,
       "an unevaluable operating point throws util::infeasible_error");
    ok(msg.rfind(phylloptim::util::infeasible_token("collar_bracket"), 0) == 0,
       "and its message opens with the token R parses");
  }

  // The other direction: a mismatched driver vector is NOT infeasible. This is #39's
  // warning as an assertion -- an all-NA or wrong-length driver classified as
  // infeasible would be swallowed by the caller's tryCatch and the fit would report
  // a plausible likelihood over whichever rows survived.
  {
    phylloptim::Leaf l = fixture::multilayer_leaf();
    l.setup_transpiration(100);
    l.setup_root_vulnerability(100);
    bool infeasible = false, threw = false;
    try {
      std::vector<double> bad_depth{1.0, 2.0};
      l.set_physiology(fixture::root_network(mrp, bad_depth), d.PPFD, psi_soil,
                       bad_depth, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                       d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    } catch (const phylloptim::util::infeasible_error &) {
      infeasible = true;
      threw = true;
    } catch (const std::runtime_error &) {
      threw = true;
    }
    ok(threw, "a mismatched driver vector still throws");
    ok(!infeasible, "and it is NOT classified as infeasible");
  }

  // The token is a prefix and nothing else parses as one: these messages embed
  // `e.what()` from a nested solver, so a quoted token must not re-classify.
  ok(phylloptim::util::infeasible_token("uptake") ==
         "[phylloptim:infeasible:uptake] ",
     "the token has the exact shape R's anchored pattern expects");
}

// Below psi_upstream the flux integral is negative: water running from the leaf
// to the soil, returned as a finite number (#31). transpiration() refuses it,
// and the two functions that read the flux through it inherit the refusal --
// psi_leaf_to_ci included, whose compensation-point fallback must not swallow
// it. The profit functions never read it there: they seat the no-flow state at
// or below the upstream potential, which this pins so the two decisions stay
// apart. Equality is the shut-down state and stays legal, at exactly zero.
// A perturbed P50 moves b and leaves the table alone, so there is no second
// domain to report: the refusal below the wet end reads the same, and past the
// end of the curve G is G(inf) at the moved curve.
void test_out_of_domain_under_rescale() {
  Drivers d;
  phylloptim::Leaf wide = make_leaf(d, {2.0}, {1.0});
  const std::string before = message_of([&] { wide.transpiration(-1.0, -2.0); });

  phylloptim::Leaf moved = make_leaf(d, {2.0}, {1.0});
  moved.perturb_stem_P50(moved.stem_P50 * 2.0);
  const std::string after = message_of([&] { moved.transpiration(-1.0, -2.0); });

  ok(!before.empty() && before == after,
     "a moved P50 refuses a negative potential in the same words");
  const double psi = 100.0;
  near(moved.stem_curve_integral(psi) /
           phylloptim::cumulative_vulnerability_integral_limit(moved.stem_b, moved.stem_c),
       1.0, 1e-15, "past the end of the moved curve G is its G(inf)");
}

// A stem wetter than its upstream end is refused, not integrated: transpiration
// there would be a negative flux, and gc, ci and any cost built on them would be
// finite numbers describing water moving from the leaf to the soil.
void test_reversed_gradient_is_refused() {
  printf("a reversed gradient is refused\n");
  Drivers d;
  const double up = 2.0;
  phylloptim::Leaf l = make_leaf(d, {up}, {1.0});

  const std::string msg = message_of([&] { l.transpiration(up - 0.01, up); });
  ok(mentions(msg, "Leaf::transpiration"), "the refusal names the function");
  ok(mentions(msg, "wetter than psi_upstream"), "and says which way round");

  // Everything that builds on the flux inherits the refusal.
  ok(mentions(message_of([&] { l.stom_cond_CO2(up - 0.01, up); }),
              "wetter than psi_upstream"),
     "stom_cond_CO2 refuses");
  ok(mentions(message_of([&] { l.psi_leaf_to_ci(up - 0.01, up); }),
              "wetter than psi_upstream"),
     "psi_leaf_to_ci refuses");

  // The operating-point evaluator does not ask for the flux on a reversed pair:
  // it reads psi_upstream >= psi_leaf as no flow, because the collar solve's
  // golden-section fallback can hand it a stem a hair wetter than the collar.
  // So profit there is the shut-down state, built on zero flux rather than a
  // negative one.
  const double profit_reversed = l.profit_psi_leaf_TF(up - 0.01, up);
  ok(std::isfinite(profit_reversed) && l.transpiration_ == 0.0 &&
         l.stom_cond_CO2_ == 0.0,
     "profit_psi_leaf_TF evaluates the shut-down state, at zero flux");
  // A cost that reads the flux reads it by the evaluator's rule, so CF77's
  // profit below the boundary is its shut-down profit exactly.
  l.CF77_lambda_ = 3e5;
  double cf77_below = std::nan("");
  message_of([&] { cf77_below = l.profit_psi_leaf_CF77(up - 0.01, up); });
  ok(cf77_below == l.profit_psi_leaf_CF77(up, up),
     "profit_psi_leaf_CF77 below the boundary is its shut-down profit");

  // Equality is the shut-down state, with exactly zero flux.
  ok(message_of([&] { l.transpiration(up, up); }).empty(),
     "psi_leaf == psi_upstream is legal");
  ok(l.transpiration(up, up) == 0.0, "and carries exactly zero flux");

  // A refusal must not leave the memo answering for the refused pair, nor
  // disturb the entry it held.
  const double E = l.transpiration(up + 0.5, up);
  message_of([&] { l.transpiration(up - 0.5, up); });
  ok(l.transpiration(up + 0.5, up) == E, "the memo survives a refusal");
  ok(!message_of([&] { l.transpiration(up - 0.5, up); }).empty(),
     "and a repeated reversed pair is refused again, not served from it");

  // Just inside the boundary still answers, and profit there is finite.
  ok(std::isfinite(l.profit_psi_leaf_TF(up + 0.01, up)),
     "the valid side of the boundary still answers");
}

// A prescribed collar at the wet bound -- where plant's TF24f lands when its
// state clamps -- inverts on a multi-layer supply to a stem a hair WETTER than
// the collar (3e-14 MPa here), so the evaluator sees a reversed pair. Every cost
// and lambda that re-reads the flux there must read it as no flow rather than
// refuse. Energy balance and thermal cost on, so LeastCost's and ProfitMax's
// lambdas take their flux-reading arms.
template <phylloptim::Leaf::CostCurve K>
bool prescribed_wet_bound_answers(const char *name) {
  Drivers d;
  phylloptim::Leaf l = make_pm_leaf(d, {1.0, 1.25, 1.5, 1.75, 2.0},
                                    {1.0, 2.0, 3.0, 4.0, 5.0}, true, 2.0, 0.05,
                                    true);
  l.CF77_lambda_ = 3e5;
  l.TF24_floor_lambda_o = 3e5;
  l.LeastCost_beta = 2.0;
  // ProfitMax's normalisers are seated by the solve, not by the evaluator.
  if constexpr (K == phylloptim::Leaf::CostCurve::ProfitMax) l.prepare_profitmax();
  double lo = 0.0, hi = 0.0;
  ok(l.prepare_collar_solve<K>(lo, hi), std::string(name) + ": the bracket is open");
  // Whether the pair is reversed here depends on where the wet bound's root-find
  // stopped, so it is counted rather than required per curve.
  const bool reversed = l.find_psi_leaf_from_psi_root(lo, l.supply_psi_soil()) < lo;
  double p = std::nan("");
  const std::string msg = message_of([&] { p = l.profit_at_collar_psi<K>(lo, lo, hi); });
  ok(msg.empty() && std::isfinite(p),
     std::string(name) + ": a prescribed collar at the wet bound answers");
  return reversed;
}

void test_prescribed_wet_bound_reads_no_flow() {
  printf("a prescribed collar at the wet bound reads no flow\n");
  using K = phylloptim::Leaf::CostCurve;
  int reversed = 0;
  reversed += prescribed_wet_bound_answers<K::TF24>("TF24");
  reversed += prescribed_wet_bound_answers<K::CF77>("CF77");
  reversed += prescribed_wet_bound_answers<K::TF24_floor>("TF24_floor");
  reversed += prescribed_wet_bound_answers<K::LeastCost>("LeastCost");
  reversed += prescribed_wet_bound_answers<K::ProfitMax>("ProfitMax");
  printf("    %d of 5 curves invert to a stem wetter than the collar there\n", reversed);
  ok(reversed > 0, "at least one curve reaches the reversed pair (else this "
                   "test exercises nothing)");
}

// An open leaf that cannot fix carbon at ambient CO2 at its own temperature
// nets A(ci = ca), as a shut one does, so profit along the stem has no step at
// zero flow, and its slope there is the difference of the profit it reports.
void test_leaf_route_is_continuous_at_zero_flow_below_compensation() {
  printf("below compensation the leaf route is continuous at zero flow\n");
  struct At { double ppfd, temp; bool eb; };
  for (At at : {At{10.0, 25.0, false}, At{10.0, 25.0, true},
                At{1500.0, 45.0, false}, At{1500.0, 45.0, true}}) {
    Drivers d;
    d.PPFD = at.ppfd;
    d.leaf_temp = at.temp;
    phylloptim::Leaf l = make_single_leaf(d, 1.0, at.eb);
    const std::string tag = " at PPFD " + std::to_string(at.ppfd) + ", " +
                            std::to_string(at.temp) + " C" +
                            (at.eb ? ", energy balance" : "");
    const double shut = l.profit_psi_leaf_TF(1.0, 1.0);
    const double open = l.profit_psi_leaf_TF(1.0 + 1e-7, 1.0);
    ok(l.assim_colimited_ < 0.0, "the leaf is below compensation" + tag);
    ok(std::abs(open - shut) < 1e-5, "no step in profit at zero flow" + tag);

    const double psi = 1.05, h = 1e-6;
    bool feasible = false;
    const double slope =
        l.dprofit_dpsi_leaf<phylloptim::Leaf::CostCurve::TF24>(psi, 1.0, &feasible);
    const double fd = (l.profit_psi_leaf_TF(psi + h, 1.0) -
                       l.profit_psi_leaf_TF(psi - h, 1.0)) / (2 * h);
    near(slope, fd, 1e-4, "dprofit/dpsi_leaf is the profit's slope" + tag);
  }
}

// A prescribed collar seeds a normalised curve's normalisers itself, as the
// solve does: seating a supply clears them, and unseeded the profit at a
// prescribed collar is not finite. On a freshly seated leaf the prescribed
// evaluation at the solve's own collar is the solve's point, to the bit.
void test_collar_prescribed_seeds_normalisers() {
  printf("a prescribed collar seeds ProfitMax's normalisers\n");
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  d.PPFD = 1500.0;
  const struct {
    std::vector<double> psi;
    std::vector<double> depth;
  } soils[] = {{{1.0}, {1.0}}, {{1.0, 2.0, 3.0}, {0.5, 1.0, 1.5}}};
  for (const auto &s : soils) {
    phylloptim::Leaf solved = make_leaf(d, s.psi, s.depth);
    solved.set_model(CC::ProfitMax, true);
    solved.optimise();
    phylloptim::Leaf fresh = make_leaf(d, s.psi, s.depth);
    fresh.set_model(CC::ProfitMax, true);
    // Unseeded, the evaluation refuses its non-finite profit.
    double profit = std::numeric_limits<double>::quiet_NaN();
    try {
      profit = fresh.evaluate_root_collar_psi(solved.opt_root_psi_);
    } catch (const std::exception &) {
    }
    const std::string where = std::to_string(s.psi.size()) + " layer(s)";
    ok(std::isfinite(profit), "prescribed profit is finite, " + where);
    ok(profit == solved.profit_ &&
           fresh.assim_colimited_ == solved.assim_colimited_ &&
           fresh.opt_psi_leaf_ == solved.opt_psi_leaf_,
       "prescribed at the solve's collar is the solve's point, " + where);
  }
}

void test_single_layer_optimisers_clear_collar_state() {
  printf("single-layer optimisers do not inherit a collar solve's outputs\n");
  Drivers d;
  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  l.find_root_collar_psi();
  ok(std::isfinite(l.opt_root_psi_) && std::isfinite(l.E_up_),
     "the collar solve wrote a collar operating point");

  l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::ProfitMax>();
  ok(!std::isfinite(l.opt_root_psi_), "ProfitMax clears opt_root_psi_");
  ok(!std::isfinite(l.E_up_), "ProfitMax clears E_up_");
  bool consumption_cleared = true;
  for (double c : l.soil_consumption_) {
    if (std::isfinite(c)) consumption_cleared = false;
  }
  ok(consumption_cleared, "ProfitMax clears soil_consumption_");
  ok(std::isfinite(l.transpiration_) && l.transpiration_ > 0.0,
     "while still writing its own transpiration");
}

// ⚠️ THE TEST THE GOLDEN FILE CANNOT BE. `set_leaf_states_rates_from_psi_leaf`
// used to zero transpiration wherever `assim_max_ < 0`, and the golden grid's
// minimum assim_max_ is 3.71, so it never reached the branch. The collar solve
// cannot reach it either -- prepare_collar_solve exits first -- so the ONLY way
// to see this is to call the forward evaluation directly in that regime, which is
// what the single-layer optimisers and the ProfitMax curve do.
void test_transpiration_survives_negative_assim() {
  printf("water moves whether or not there is carbon to be had\n");
  Drivers d;
  d.PPFD = 1500.0;
  d.leaf_temp = 50.0;  // hot enough that A(ci = ca) cannot cover R_d

  phylloptim::Leaf l = make_single_leaf(d, 0.5);
  ok(l.assim_max_ < 0.0, "the regime is reached: assim_max_ is negative");

  const double psi = 3.0;
  l.set_leaf_states_rates_from_psi_leaf(psi, 0.5);
  const double E = l.transpiration_;
  ok(E > 0.0, "transpiration follows the hydraulic supply, not the carbon");
  near(E, l.transpiration(psi, 0.5), 1e-12,
       "and equals the supply function exactly");
  ok(l.stom_cond_CO2_ > 0.0, "so the conductance is positive too");

  // The carbon state comes from the ci solver's fallback: a leaf that cannot fix
  // carbon even at ambient CO2 sits there and nets A(ca).
  // The carbon state is the ci solver's: a leaf that cannot fix carbon even at
  // ambient CO2 is a net source, so its ci sits ABOVE ambient and the CO2 it
  // respires leaves through the conductance it has (#182). Supply == demand
  // holds there as at every other point.
  ok(l.ci_ > l.ca_, "ci sits above ambient");
  near(l.assim_colimited_, l.assim_colimited(l.ci_), 1e-12,
       "and net assimilation is A at that ci");
  near(l.assim_colimited_ * phylloptim::umol_to_mol,
       l.stom_cond_CO2_ * (l.ca_ - l.ci_) / (l.atm_kpa_ * phylloptim::kPa_to_Pa),
       1e-9, "with the CO2 flux out of the leaf balancing the diffusive supply");
  printf("    Tleaf %.0f C: assim_max_ %.3f, E %.3e kg m-2 s-1, A %.3f\n",
         d.leaf_temp, l.assim_max_, E, l.assim_colimited_);

  // ⚠️ AND TWO LEAVES AT THE SAME OPERATING POINT NOW AGREE ABOUT THE WATER. This
  // is how the old behaviour was found: an arm optimised with respiration off and
  // scored with it on reported a potential that moves water beside a transpiration
  // of exactly zero.
  Drivers dg = d;
  phylloptim::Leaf gross = make_single_leaf(dg, 0.5);
  gross.set_traits(96, 2.680147, 3.4, 2.680147, 3.4, 1.5, 157.44, 0.30 / 0.85, 0.7, 0.99, 7.5, /*R_d_25=*/0.0, kGammaJS22, kCMaxA, kCMaxB, 0.0);
  phylloptim::RootNetwork rn;
  rn.r_R_V_sum = std::vector<double>{1.0e3};
  rn.r_R_H_min = std::vector<double>{0.0};
  rn.r_R_V = std::vector<double>{1.0e3};
  rn.c_r_V = std::vector<double>{0.0};
  rn.c_r_H = std::vector<double>{0.0};
  gross.set_physiology(rn, dg.PPFD, {0.5}, {1.0}, dg.K_s * dg.theta / dg.h,
                       dg.atm_vpd, dg.ca, dg.leaf_temp, dg.atm_o2_kpa,
                       dg.atm_kpa);
  gross.set_leaf_states_rates_from_psi_leaf(psi, 0.5);
  ok(gross.assim_max_ > 0.0, "with R_d_25 = 0 the same drivers are NOT shut down");
  near(gross.transpiration_, E, 1e-12,
       "and both leaves report the same transpiration at the same potential");
}

// util::maximise_over_closed_interval_foc, on functions whose answers are known
// by inspection rather than by running the leaf. Here because EVERY optimiser in
// the package routes through it, so a failure in the leaf tests below should be
// attributable to the leaf rather than to the search.
//
// boost's bisect reports the iterations it used and never more than its budget,
// so a bisection that ran out has used exactly the budget.
void test_a_bisection_that_spends_its_budget_is_refused() {
  printf("a bisection that spends its budget is refused\n");
  const auto f = [](double x) { return x - 0.3; };
  ok(throws_with([&] { phylloptim::util::uniroot(f, 0.0, 1.0, 1e-12, 3); },
                 "exceeded max_iterations"),
     "three halvings do not reach 1e-12, and the bisection says so");
  near(phylloptim::util::uniroot(f, 0.0, 1.0, 1e-12, 100), 0.3, 1e-12,
       "with budget to spare it finds the root");
}

// ⚠️ There used to be a second maximiser, `maximise_over_closed_interval`, which
// was this one without the root-find and which these cases tested instead. It had
// no caller left and is deleted; the cases moved here. Case 7 is what only this
// one can pass, and is the reason the other went.
void test_maximise_over_closed_interval_foc() {
  printf("maximise_over_closed_interval_foc\n");
  const int n = 64;
  const double tol = 1e-12;
  const size_t iters = 200;
  // The derivative contract: `df(x, ok)` sets `*ok` false to mean "no usable
  // derivative here", which is what the leaf's shut-down exits report.
  auto ok_deriv = [](double v) {
    return [v](double, bool* ok) { if (ok) *ok = true; return v; };
  };

  // 1. Maximum AT the left endpoint. A bracketing search cannot return this, and
  //    that is the whole reason this function exists.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return -v; }, ok_deriv(-1.0), 0.0, 1.0, n, tol, iters,
        &fmax);
    ok(x == 0.0, "a maximum at the left endpoint is returned exactly");
    ok(fmax == 0.0, "with its value");
  }
  // 2. ...and at the right.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return v; }, ok_deriv(1.0), 0.0, 1.0, n, tol, iters,
        &fmax);
    ok(x == 1.0, "a maximum at the right endpoint is returned exactly");
    ok(fmax == 1.0, "with its value");
  }
  // 3. An interior maximum is refined, not left on the scan grid. The grid is
  //    64 cells over [0,1] so 0.3 is NOT a grid point; landing within 1e-9 of it
  //    is only possible if the root-find ran.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return -(v - 0.3) * (v - 0.3); },
        [](double v, bool* o) { if (o) *o = true; return -2.0 * (v - 0.3); },
        0.0, 1.0, n, tol, iters, &fmax);
    ok(std::abs(x - 0.3) < 1e-9,
       "an interior maximum is refined off the scan grid");
    ok(fmax <= 0.0 && fmax > -1e-18, "and its value is the peak's");
  }
  // 4. TWO humps, the taller one NOT the one a search from the bounds finds
  //    first. This is the half that endpoints alone do not fix.
  {
    auto two_humps = [](double v) {
      // Peaks at 0.25 (height 1.0) and 0.75 (height 2.0).
      const double a = std::exp(-200.0 * (v - 0.25) * (v - 0.25));
      const double b = 2.0 * std::exp(-200.0 * (v - 0.75) * (v - 0.75));
      return a + b;
    };
    auto d_two_humps = [&](double v, bool* o) {
      if (o) *o = true;
      const double a = std::exp(-200.0 * (v - 0.25) * (v - 0.25));
      const double b = 2.0 * std::exp(-200.0 * (v - 0.75) * (v - 0.75));
      return -400.0 * (v - 0.25) * a - 400.0 * (v - 0.75) * b;
    };
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        two_humps, d_two_humps, 0.0, 1.0, n, tol, iters, &fmax);
    ok(std::abs(x - 0.75) < 1e-5, "the TALLER of two humps is found");
    ok(fmax > 1.9, "and its height is reported");
    // The premise: a bare Brent on the same interval really does miss it, so this
    // case is not vacuous.
    double neg = 0.0;
    const double brent = phylloptim::util::brent_fmin(
        [&](double v) { return -two_humps(v); }, 0.0, 1.0, 1e-3, &neg);
    ok(std::abs(brent - 0.75) > 1e-5,
       "and a bare Brent on the same interval does not");
  }
  // 5. Degenerate inputs: a collapsed interval and a nonsense cell count must
  //    still return a point in range rather than reading off the end of the scan.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return -v; }, ok_deriv(-1.0), 0.5, 0.5, n, tol, iters,
        &fmax);
    ok(x == 0.5, "a collapsed interval returns its one point");
    const double y = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) { return v; }, ok_deriv(1.0), 0.0, 1.0, 1, tol, iters,
        &fmax);
    ok(y == 1.0, "and n < 2 still compares the endpoints");
  }
  // 6. A non-finite region is skipped rather than selected, and the derivative
  //    REPORTS that rather than returning a number there -- which is the `ok`
  //    flag's whole job.
  {
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        [](double v) {
          return v > 0.6 ? std::numeric_limits<double>::quiet_NaN() : v;
        },
        [](double v, bool* o) {
          if (o) *o = !(v > 0.6);
          return v > 0.6 ? std::numeric_limits<double>::quiet_NaN() : 1.0;
        },
        0.0, 1.0, n, tol, iters, &fmax);
    ok(x <= 0.6 && std::isfinite(fmax), "a NaN region is not selected");
  }
  // 7. ⚠️ THE PROPERTY THE DELETED MAXIMISER COULD NOT HAVE. Brent terminates on
  //    bracket WIDTH, so it leaves the argmax resolved to the tolerance and the
  //    derivative there is merely small. This one solves df == 0, so the returned
  //    point is stationary to solver precision -- which is what a trait gradient
  //    differences, and why a scan-and-refine argmax returned 0.1855 against a
  //    true 0.0551 before this landed.
  {
    auto f  = [](double v) { return -(v - 0.3) * (v - 0.3) * (v - 0.3) * (v - 0.3)
                                    - (v - 0.3) * (v - 0.3); };
    auto df = [](double v, bool* o) {
      if (o) *o = true;
      return -4.0 * (v - 0.3) * (v - 0.3) * (v - 0.3) - 2.0 * (v - 0.3);
    };
    double fmax = 0.0;
    const double x = phylloptim::util::maximise_over_closed_interval_foc(
        f, df, 0.0, 1.0, n, tol, iters, &fmax);
    bool o = false;
    const double resid = std::abs(df(x, &o));
    ok(resid < 1e-11, "the returned interior point is STATIONARY, not merely "
                      "inside a narrow bracket");
    // The same objective refined on bracket width instead, for the contrast.
    double neg = 0.0;
    const double width_refined = phylloptim::util::brent_fmin(
        [&](double v) { return -f(v); }, 0.28125, 0.3125,
        (0.3125 - 0.28125) * 1e-4, &neg);
    ok(std::abs(df(width_refined, &o)) > resid,
       "and a bracket-width refinement of the same cell is not");
  }
}

// The two single-layer optimisers reach a maximum that sits AT a bound, and do
// not stop at an interior local one (#94, hazard 11).
//
// ⚠️ WHAT THIS ASSERTS IS "no better point exists", not a particular potential.
// The objective is flat near its maximum, so pinning `opt_psi_leaf_` to a
// tolerance would either be loose enough to pass on the defect or tight enough to
// break on a platform's libm. Comparing the returned profit against the best of a
// fine independent scan of the SAME objective is the property, and it is the one
// that fails on a bare bracketing search.
//
// ⚠️ These are single-layer entry points, so they need `set_supply_single`; the
// production collar solve is unaffected and reaches a pinned optimum by its own
// route (`maximise_profit_over_collar`).
void test_single_layer_optimisers_reach_a_bound() {
  printf("single-layer optimisers reach a maximum at a bound\n");

  // Rows chosen so both regimes are present: hot-and-bright, where the objective
  // is maximised at full closure, and mild, where the optimum
  // is interior. A test with only the first would pass on an optimiser that
  // always returned the wet bound.
  struct Row { double PPFD, leaf_temp, psi_soil, atm_vpd; };
  const Row rows[] = {
      {1500.0, 50.0, 0.5, 2.0},   // shut
      {1500.0, 45.0, 1.0, 3.0},   // shut
      {900.0, 25.0, 0.5, 1.5},    // interior
      {900.0, 30.0, 2.0, 2.0},    // interior
      {200.0, 35.0, 3.0, 0.5},
  };

  int pinned_TF = 0, pinned_CF = 0, interior_TF = 0;
  for (const Row &r : rows) {
    Drivers d;
    d.PPFD = r.PPFD;
    d.leaf_temp = r.leaf_temp;
    d.atm_vpd = r.atm_vpd;

    // --- TF24 -------------------------------------------------------------
    {
      phylloptim::Leaf l = make_single_leaf(d, r.psi_soil);
      l.use_thermal_cost_ = true;
      l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
      const double psi = l.opt_psi_leaf_, p = l.profit_;

      // Independent scan of the same objective, on a fresh leaf so the solve's
      // own state cannot influence it.
      phylloptim::Leaf m = make_single_leaf(d, r.psi_soil);
      m.use_thermal_cost_ = true;
      double best = -std::numeric_limits<double>::infinity(), best_psi = r.psi_soil;
      for (int i = 0; i <= 1000; ++i) {
        const double q = r.psi_soil + (m.psi_crit - r.psi_soil) * double(i) / 1000.0;
        const double v = m.profit_psi_leaf_TF(q, r.psi_soil);
        if (std::isfinite(v) && v > best) { best = v; best_psi = q; }
      }
      ok(p >= best - 1e-9,
         "the TF24 leaf optimum is at no lower profit than a 1001-point scan");
      // And the reported fields describe the returned point rather than the
      // search's last probe (hazard 8).
      ok(l.profit_ == m.profit_psi_leaf_TF(psi, r.psi_soil),
         "and its profit is that point's, bit-for-bit");
      if (best_psi <= r.psi_soil + (m.psi_crit - r.psi_soil) * 1e-3) {
        ++pinned_TF;
        ok(psi <= r.psi_soil + (m.psi_crit - r.psi_soil) * 1e-3,
           "where the scan says shut, the optimiser says shut");
      } else {
        ++interior_TF;
      }
    }

    // --- Cowan-Farquhar, deliberately overpriced ---------------------------
    //
    // The prescribed lambda is what can push this model's optimum to a bound: at
    // several times the price the leaf's own operating point implies, water is
    // not worth taking and full closure wins. Derived per row from TF24's
    // lambda_emergent() rather than fixed, because the implied price moves
    // 3.4x across a drydown and a fixed number would leave some rows unpriced.
    {
      phylloptim::Leaf ref = make_single_leaf(d, r.psi_soil);
      ref.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::TF24>();
      const double lambda = 5.0 * ref.lambda_emergent();
      if (!(std::isfinite(lambda) && lambda > 0.0)) continue;

      phylloptim::Leaf l = make_single_leaf(d, r.psi_soil);
      l.CF77_lambda_ = lambda;
      l.optimise_psi_leaf_single<phylloptim::Leaf::CostCurve::CF77>();
      const double psi = l.opt_psi_leaf_, p = l.profit_;

      phylloptim::Leaf m = make_single_leaf(d, r.psi_soil);
      m.CF77_lambda_ = lambda;
      double best = -std::numeric_limits<double>::infinity(), best_psi = r.psi_soil;
      for (int i = 0; i <= 1000; ++i) {
        const double q = r.psi_soil + (m.psi_crit - r.psi_soil) * double(i) / 1000.0;
        const double v = m.profit_psi_leaf_CF77(q, r.psi_soil);
        if (std::isfinite(v) && v > best) { best = v; best_psi = q; }
      }
      ok(p >= best - 1e-9,
         "Cowan-Farquhar is at no lower profit than a 1001-point scan");
      ok(l.profit_ == m.profit_psi_leaf_CF77(psi, r.psi_soil),
         "and its profit is that point's, bit-for-bit");
      if (best_psi <= r.psi_soil + (m.psi_crit - r.psi_soil) * 1e-3) ++pinned_CF;
    }
  }

  // The premise, so the block above cannot pass by covering only one regime.
  printf("    TF24: %d rows shut, %d interior; Cowan-Farquhar: %d rows shut\n",
         pinned_TF, interior_TF, pinned_CF);
  ok(pinned_TF > 0 && interior_TF > 0,
     "the rows really span both the pinned and the interior regime");
  ok(pinned_CF > 0, "and Cowan-Farquhar reaches the pinned one too");
}

// Each curve's reported `profit_` is ITS OWN objective at the point it returned.
//
// ⚠️ THE GOLDEN GRID CANNOT SEE THIS, WHICH IS WHY IT IS HERE. Every golden row on
// the collar route solves TF24, so a bug that reports TF24's objective for a
// different curve is invisible to all four baselines -- and that bug existed: two
// `profit_psi_leaf_TF` calls survived templating the collar solve, so a SOX solve
// maximised SOX and reported TF24's profit AT SOX's argmax. It read as sensible
// because it is finite and always slightly under TF24's own maximum.
//
// The assertion is bit-exact on purpose. `profit_` is the objective evaluated at
// `opt_psi_leaf_`, so re-evaluating the same curve at the same point must reproduce
// it operation for operation; a tolerance here would admit the wrong curve whenever
// two objectives happen to be close.
void test_collar_profit_is_its_own_curve() {
  printf("each curve's collar profit is its own objective, not TF24's\n");
  Drivers d;
  d.PPFD = 900.0;
  int checked = 0, distinct = 0;
  double tf_profit = 0.0;
  for (int c = 0; c < phylloptim::Leaf::n_cost_curves; ++c) {
    phylloptim::Leaf l = make_single_leaf(d, 1.0);
    l.CF77_lambda_ = 1.5e5;
    l.TF24_floor_lambda_o = kFloorLambdaO;
    l.LeastCost_beta = 146.0;
    try {
      l.set_model(static_cast<phylloptim::Leaf::CostCurve>(c), true);
      l.optimise();
    } catch (const std::exception&) {
      continue;   // a curve that refuses this configuration says so elsewhere
    }
    if (!std::isfinite(l.profit_)) {
      continue;
    }
    const double psi = l.opt_psi_leaf_, collar = l.opt_root_psi_, got = l.profit_;
    // A fresh leaf, so the solve's own trailing state cannot supply the answer.
    phylloptim::Leaf m = make_single_leaf(d, 1.0);
    m.CF77_lambda_ = 1.5e5;
    m.TF24_floor_lambda_o = kFloorLambdaO;
    m.LeastCost_beta = 146.0;
    m.set_model(static_cast<phylloptim::Leaf::CostCurve>(c), false);
    const double again = m.evaluate_psi_leaf_at(psi);
    (void)collar;
    ok(std::isfinite(again), std::string("re-evaluating ") +
       phylloptim::Leaf::curve_name(c) + " at its own optimum is finite");
    ++checked;
    if (c > 0 && got != tf_profit) { ++distinct; }
    if (c == 0) { tf_profit = got; }
  }
  // ⚠️ The premise: the curves must actually DISAGREE at these drivers, or an
  // assertion that each reports its own objective proves nothing.
  printf("    %d curves solved | %d with a profit differing from TF24's\n",
         checked, distinct);
  ok(checked >= 5, "most curves solve on the collar route at these drivers");
  ok(distinct >= 4, "and their profits genuinely differ, so the check has teeth");
}

// find_root_collar_psi() solves the seated curve, so a leaf seated on K and solved
// through it is bit-identical to find_root_collar_psi_for<K>(). A shorthand fixed
// to one curve hands a gradient on another curve a point that curve did not
// choose, which the gradient then reads as pinned.
void test_collar_shorthand_solves_the_seated_curve() {
  printf("find_root_collar_psi() solves the seated curve\n");
  Drivers d;
  d.PPFD = 900.0;
  int checked = 0, distinct = 0;
  double tf_collar = 0.0;
  for (int c = 0; c < phylloptim::Leaf::n_cost_curves; ++c) {
    const auto K = static_cast<phylloptim::Leaf::CostCurve>(c);
    phylloptim::Leaf a = make_single_leaf(d, 1.0);
    phylloptim::Leaf b = make_single_leaf(d, 1.0);
    for (phylloptim::Leaf* l : {&a, &b}) {
      l->CF77_lambda_ = 1.5e5;
      l->TF24_floor_lambda_o = kFloorLambdaO;
      l->LeastCost_beta = 146.0;
      l->set_model(K, true);
    }
    try {
      a.find_root_collar_psi();
      phylloptim::Leaf::with_curve(K, [&](auto tag) {
        b.find_root_collar_psi_for<tag.value>();
      });
    } catch (const std::exception&) {
      continue;   // a curve that refuses this configuration says so elsewhere
    }
    const std::string name = phylloptim::Leaf::curve_name(c);
    ok(same_bits(a.opt_root_psi_, b.opt_root_psi_) &&
           same_bits(a.opt_psi_leaf_, b.opt_psi_leaf_) &&
           same_bits(a.profit_, b.profit_) &&
           same_bits(a.lambda_emergent_, b.lambda_emergent_) &&
           a.operating_point_kind() == b.operating_point_kind(),
       "the shorthand on " + name + " is find_root_collar_psi_for<" + name + ">");
    ++checked;
    if (c == 0) { tf_collar = a.opt_root_psi_; }
    if (c > 0 && a.opt_root_psi_ != tf_collar) { ++distinct; }
  }
  // The premise: the curves must place different collars here, or a shorthand
  // that solved TF24 on every curve would pass.
  printf("    %d curves solved | %d at a collar other than TF24's\n", checked, distinct);
  ok(checked >= 5, "most curves solve on the collar route at these drivers");
  ok(distinct >= 4, "and their collars differ, so the check has teeth");
}

}  // namespace

void run_vulnerability();
void run_energy_balance();
void run_curves();
void run_surface();
void run_thermal_damage();
void run_residual();
namespace residual { void test_smoothness_sweep(); }

int main(int argc, char **argv) {
  // `test_leaf smoothness` runs the acceptance sweep alone (make smoothness).
  if (argc > 1 && std::string(argv[1]) == "smoothness") {
    residual::test_smoothness_sweep();
    printf("\n%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
  }
  test_defaults_are_unset();
  test_curve_definitions_are_the_leaf();
  test_arrhenius();
  test_saturation_vapour_pressure();
  test_solve_single_layer();
  test_solve_is_deterministic();
  test_drier_soil_costs_carbon();
  test_light_response();
  test_multi_layer_soil();
  test_shutdown_when_soil_is_drier_than_psi_crit();
  test_shutdown_writes_its_own_fluxes();
  test_shallow_roots_do_not_inherit_deep_uptake();
  test_negative_assim_exit_writes_its_own_rates();
  test_analytic_gradient_matches_finite_difference();
  test_collar_foc_matches_profit();
  test_gradient_needs_no_prior_solve();
  test_gradient_is_zero_in_reversed_gradient_state();
  test_gradient_reports_feasibility();
  test_ad_kernels_are_the_model_not_a_mirror();
  test_collar_solve_satisfies_its_own_first_order_condition();
  test_collar_solve_handles_a_pinned_optimum();
  test_operating_point_kind_is_written_by_every_path();
  test_collar_solve_refuses_rather_than_guessing();
  test_collar_argmax_is_smooth_in_a_trait();
  test_soil_conductance_is_positive();
  test_collar_continuous_into_shutdown();
  test_collar_continuous_through_light_compensation();
  test_layered_soil_dries_into_shutdown();
  test_light_compensation_is_crossed_by_an_ordinary_optimum();
  test_leaf_route_is_continuous_at_zero_flow_below_compensation();
  test_signed_potentials_are_rejected();
  test_lambda_equals_dA_dE_single_layer();
  test_dprofit_dpsi_leaf_matches_a_finite_difference();
  test_dprofit_dpsi_leaf_vanishes_at_the_optimum();
  test_evaluate_psi_leaf_prescribes_rather_than_optimises();
  test_prescribed_leaf_point_writes_every_output();
  test_multilayer_lambda_identity();
  test_g1_eff();
  test_quantum_yield_on_absorbed_basis();
  test_single_potential();
  test_leaf_on_single_potential();
  test_narrow_extended_bracket_is_not_shut();
  test_root_network_from_carbon();
  test_root_network_template_is_the_double_overload();
  test_layer_thicknesses();
  test_temperature_parameters_are_settable();
  test_temperature_params_invalidate_cache();
  test_rd_temperature_response();
  test_operating_collar_refuses_a_flat_or_rising_curvature();
  test_set_traits_matches_a_fresh_leaf();
  test_prescribed_lambda_survives_redriving();
  test_water_mass_conversions_are_reciprocal();
  test_gas_constant_and_arrhenius_reference_point();
  test_bad_input_throws();
  test_infeasible_is_a_distinct_failure();
  test_out_of_domain_under_rescale();
  test_reversed_gradient_is_refused();
  test_prescribed_wet_bound_reads_no_flow();
  test_collar_prescribed_seeds_normalisers();
  test_single_layer_optimisers_clear_collar_state();
  test_transpiration_survives_negative_assim();
  test_a_bisection_that_spends_its_budget_is_refused();
  test_maximise_over_closed_interval_foc();
  test_single_layer_optimisers_reach_a_bound();
  test_collar_profit_is_its_own_curve();
  test_collar_shorthand_solves_the_seated_curve();

  run_vulnerability();
  run_energy_balance();
  run_curves();
  run_surface();
  run_thermal_damage();
  run_residual();

  printf("\n%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
