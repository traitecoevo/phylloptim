// Leaf temperature, the deficit at the leaf, the boundary layer and the longwave balance.

#include "leaf_harness.hpp"

namespace {

// The same check with the energy balance ON, which the one above cannot see: its
// leaf is at a prescribed temperature, so the energy-balance term of the
// first-order condition, dA/dTleaf and dgc/dT = -gc * Delta(Tleaf) / D_leaf never
// enter. With the gate on they all do.
//
// ⚠️ THE ASSERTION IS ON THE MEDIAN, AND THAT IS MEASURED. At h = 3e-4 a wrong
// slope of the saturation curve (FAO-56's 4098, 4.2e-5 below 17.27 * 237.3) put
// the median relative error at 4.0e-06 over these 34 points, never below 2.2e-07;
// the exact slope puts it at 1.5e-08. But single points of the exact build reach
// 1.3e-06, where the ci root-find's ~1e-10 noise in profit is divided by 2h on a
// hot leaf with small A, so no per-point bound separates the two. A median does,
// by 270x. The per-point bound only catches a missing or mis-signed term, which
// is an O(1) error.
void test_dprofit_dpsi_leaf_matches_a_finite_difference_energy_balance() {
  printf("dprofit/dpsi_leaf with the energy balance on: analytic vs a central difference\n");
  Drivers d;
  d.PPFD = 900.0;
  const double h = 3e-4;
  std::vector<double> errs;
  for (double leaf_temp : {20.0, 30.0, 40.0}) {
    for (double psi_soil : {0.5, 1.0, 2.0}) {
      phylloptim::Leaf l = fixture::multilayer_leaf();
      l.setup_transpiration(400);
      l.setup_root_vulnerability(400);
      l.use_energy_balance_ = true;
      l.NIR_int_ = d.PPFD / phylloptim::umol_par_per_joule;  // daylight
      l.set_supply_single(0.0);
      l.set_physiology(fixture::series_resistance(1e3), d.PPFD, {psi_soil}, {1.0},
                       d.K_s * d.theta / d.h, d.atm_vpd, d.ca, leaf_temp,
                       d.atm_o2_kpa, d.atm_kpa);
      const double lo = psi_soil, hi = l.psi_crit;
      for (double frac : {0.15, 0.35, 0.55, 0.75}) {
        const double psi = lo + frac * (hi - lo);
        bool feasible = false;
        const double an =
            l.dprofit_dpsi_leaf<phylloptim::Leaf::CostCurve::TF24>(psi, psi_soil,
                                                                   &feasible);
        if (!feasible) continue;
        const double up = l.profit_psi_leaf_TF24(psi + h, psi_soil);
        const double dn = l.profit_psi_leaf_TF24(psi - h, psi_soil);
        const double fd = (up - dn) / (2.0 * h);
        const double err = std::abs(an - fd) / std::max(std::abs(an), 1.0);
        errs.push_back(err);
        near(err, 0.0, 1e-5,
             "TF24 dprofit, energy balance on, T=" + std::to_string(leaf_temp) +
                 " psi_soil=" + std::to_string(psi_soil) +
                 " frac=" + std::to_string(frac));
      }
    }
  }
  ok(errs.size() >= 30, "enough feasible rows reached the comparison");
  std::vector<double> sorted = errs;
  std::sort(sorted.begin(), sorted.end());
  const double median = sorted[sorted.size() / 2];
  ok(median < 1e-7, "median relative error " + std::to_string(median) +
                        " is at the difference's floor, not a wrong slope's");
  printf("    %zu points | median relative error %.3e, worst %.3e\n", errs.size(),
         median, sorted.back());
}

// Tleaf as a reported OUTPUT, on both paths and out of every exit.
//
// The gap this closes: on the PM path the leaf's temperature is solved per
// operating point, was used to re-derive the whole Farquhar block, and was then
// discarded -- so the one quantity that path exists to produce was the one thing a
// caller could not read. `leaf_temp_` is not it there, because set_physiology has
// reinterpreted that driver as AIR temperature.
//
// ⚠️ The shut-down cases are the point of this test, not an afterthought. They do
// not go through set_leaf_states_rates_from_psi_leaf, so they are where hazard 8
// bites: an output a branch declines to write becomes the PREVIOUS solve's value,
// and plant drives every individual in a patch through one persistent Leaf.
void test_leaf_temperature_is_reported() {
  printf("Tleaf is reported, on both paths and from every exit\n");
  Drivers d;

  // Unsolved: the NA sentinel, like every other output.
  {
    phylloptim::Leaf bare = fixture::multilayer_leaf();
    ok(!std::isfinite(bare.Tleaf_), "Tleaf is NA before a solve");
  }

  // Off the PM path Tleaf IS the driver, and exactly so -- a tolerance here would
  // pass on a value that had been round-tripped through the energy balance.
  for (double T : {15.0, 25.0, 40.0}) {
    Drivers dt = d;
    dt.leaf_temp = T;
    dt.atm_vpd = 1.0;  // below saturation at 15 C (1.71 kPa)
    phylloptim::Leaf l = make_leaf(dt, {2.0}, {1.0});
    l.find_root_collar_psi();
    ok(l.Tleaf_ == T,
       "off the energy-balance path Tleaf is the driver, at T = " +
           std::to_string(T));
    // ⚠️ NOT `.back()`: outputs follow `Tleaf`. Addressed by its position in
    // R's `.operating_point_names`, which test-surface.R pins.
    const std::vector<double> ov = l.operating_point_values();
    ok(ov[12] == l.Tleaf_,
       "operating_point_values() reports it, at T = " + std::to_string(T));
  }

  // On the PM path it is not the driver, and it is hotter: the leaf absorbs
  // radiation and sheds only part of it as latent heat.
  {
    phylloptim::Leaf eb = make_pm_leaf(d, {2.0}, {1.0}, true);
    eb.find_root_collar_psi();
    ok(eb.Tleaf_ != eb.leaf_temp_,
       "on the energy-balance path Tleaf is not the leaf_temp driver");
    ok(eb.Tleaf_ == eb.leaf_temp_from_E(eb.transpiration_),
       "and it is the temperature the solve's own transpiration implies");
    ok(eb.Tleaf_ > eb.Tair_, "the leaf runs hotter than the air here");
  }

  // ONE leaf, driven twice, second time into shut-down: the case a fresh-object
  // test cannot see. Both paths, because the two exits differ.
  for (bool gate : {false, true}) {
    const std::string what = gate ? " (energy balance on)" : " (prescribed)";
    phylloptim::Leaf l = make_pm_leaf(d, {2.0}, {1.0}, gate);
    l.find_root_collar_psi();
    const double wet = l.Tleaf_;
    ok(std::isfinite(wet), "a wet solve reports a finite Tleaf" + what);

    // Drier than the end of the stem curve, so the collar solve takes a
    // shut-down exit.
    std::vector<double> root{1.0 / d.area_leaf}, dry{17.0}, depth{1.0};
    l.set_physiology(fixture::root_network(root, depth), d.PPFD, dry, depth,
                     d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
    l.find_root_collar_psi();
    ok(l.transpiration_ == 0.0, "the second solve really did shut down" + what);
    // The claim is not "it changed" -- it is that the value describes THIS
    // operating point. At zero transpiration there is no latent cooling, so the
    // PM answer is leaf_temp_from_E(0); off that path it is still the driver.
    const double want = gate ? l.leaf_temp_from_E(0.0) : l.leaf_temp_;
    ok(l.Tleaf_ == want,
       "a shut-down leaf reports its own temperature, not the last solve's" +
           what);
    if (gate) {
      ok(l.Tleaf_ > wet,
         "and a leaf that has stopped transpiring is hotter than one that had"
         " not" + what);
    }
  }
}

// One temperature per reported operating point, at the exits that transpire
// nothing (#105).
//
// The gap: on the PM path both shut-down exits reported a leaf temperature and a
// respiration rate belonging to DIFFERENT temperatures. They do not go through
// set_leaf_states_rates_from_psi_leaf, so they inherited the Tair baseline
// set_physiology derives and then described a state whose transpiration is zero.
//
// ⚠️ `operating_points.tsv` cannot see any of this -- the golden grid runs with
// the gate off, where Tleaf IS leaf_temp_ and there is nothing to disagree with.
// A bit-identical golden run is not evidence here, which is why this test exists
// rather than a regenerated baseline.
//
// ⚠️ The two cases move in OPPOSITE directions, so a test that only checked "R_d
// went up" would pass on half the fix. Rn is absorbed shortwave (proportional to
// PPFD) less the isothermal longwave loss, so a hydraulically shut leaf in full
// sun is hotter than the air and a dark one is cooler. The dark leaf is no
// longer an exit of its own (#182): it is a wet-pinned optimum of the ordinary
// solve, and the identities below hold there for the same reason they hold at
// every evaluated point.
void test_shutdown_reports_one_temperature() {
  printf("a shut-down leaf reports one temperature, not two\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  Drivers d;
  d.leaf_temp = 30.0;  // AIR temperature on this path

  // (psi_soil, PPFD, expected kind) -- one row per case. Reached through the
  // drivers rather than by calling the exits directly, so the test breaks if a
  // future bracket change stops routing here.
  struct Case {
    const char *name;
    double psi_soil;
    double PPFD;
    Kind kind;
    bool hotter_than_air;
  };
  const Case cases[] = {
      // Soil drier than psi_crit: no collar potential both moves water and stays
      // inside the critical potentials. Full sun, so the leaf runs hot.
      {"hydraulic shutdown", 17.0, 900.0,
       Kind::HydraulicShutdown, true},
      // Wet soil, no light: gross assimilation at ci = ca cannot cover R_d, so
      // the leaf declines to open. Rn is NEGATIVE here, so this leaf is cooler
      // than the air.
      {"dark leaf", 1.0, 1.0, Kind::BoundarySoil, false},
  };

  for (const Case &c : cases) {
    Drivers dc = d;
    dc.PPFD = c.PPFD;
    phylloptim::Leaf l = make_pm_leaf(dc, {c.psi_soil}, {1.0}, true);
    l.find_root_collar_psi();
    const std::string what = std::string(" (") + c.name + ")";

    ok(l.operating_point_kind() == c.kind,
       "the drivers really reach this case" + what);
    // The shut leaf transpires nothing; the pinned dark leaf sits a step
    // inside the wet bound at g_min = 0 and moves a rounding-level flux.
    ok(l.transpiration_ < 1e-9, "and it transpires (almost) nothing" + what);
    ok(l.Tleaf_ == l.leaf_temp_from_E(l.transpiration_),
       "Tleaf is the temperature its own flux implies" + what);

    // The premise: the two temperatures differ, so there is something to get
    // wrong. Without this the assertions below would pass on a leaf at Tair.
    ok(l.Tleaf_ != l.Tair_, "Tleaf is not air temperature" + what);
    ok(c.hotter_than_air ? (l.Tleaf_ > l.Tair_) : (l.Tleaf_ < l.Tair_),
       "and it sits on the expected side of it" + what);

    // The fix, stated as the identity it restores: the respiration in force is
    // the one the model's OWN curve gives at the temperature being reported.
    // Bit-exact -- same function, same argument, so a tolerance here would pass
    // on a value derived at some third temperature.
    const double rd_reported = l.R_d_;
    const double gamma_reported = l.gamma_;
    const double Tleaf = l.Tleaf_, Tair = l.Tair_;

    l.update_temperature_dependent_params(Tleaf);
    ok(l.R_d_ == rd_reported, "R_d is the curve's value at the reported Tleaf" + what);
    ok(l.gamma_ == gamma_reported,
       "and so is the compensation point" + what);

    l.update_temperature_dependent_params(Tair);
    ok(l.R_d_ != rd_reported,
       "and it is NOT the Tair value -- the two really differ here" + what);
  }

  // The accounting identity every branch is supposed to satisfy, re-checked at
  // these points because the shut one used to form profit by hand rather than
  // through profit_psi_leaf_TF24. `ci_` is in here because an exit once did not
  // write it at all: it reported the previous solve's internal CO2, or the NA
  // sentinel on a cold object.
  for (const Case &c : cases) {
    Drivers dc = d;
    dc.PPFD = c.PPFD;
    phylloptim::Leaf l = make_pm_leaf(dc, {c.psi_soil}, {1.0}, true);
    l.find_root_collar_psi();
    const std::string what = std::string(" (") + c.name + ")";

    ok(std::isfinite(l.ci_), "ci is written, and finite" + what);
    // The block's value at the reported ci -- or exactly zero where nothing
    // conducts, since the ci solve's residual there IS A (#182).
    ok(l.assim_colimited_ ==
           (l.stom_cond_CO2_ == 0.0 ? 0.0 : l.assim_colimited(l.ci_)),
       "net assimilation is the block's value at the reported ci" + what);
    near(l.profit_, l.assim_colimited_ - l.hydraulic_cost_TF24(l.opt_psi_leaf_),
         1e-12, "profit == assimilation - hydraulic cost" + what);
    // ci is the root of supply == demand at the reported conductance, which
    // with no conductance at all is the compensation point of this block.
    near(l.assim_colimited(l.ci_) * phylloptim::umol_to_mol -
             l.stom_cond_CO2_ * (l.ca_ - l.ci_) / (l.atm_kpa_ * phylloptim::kPa_to_Pa),
         0.0, 1e-18, "ci balances supply and demand at the reported temperature" + what);
    // #93's leaf-to-air deficit is on the same footing: reported, read by
    // g1_eff(), and derived from a temperature -- so it has to be THIS
    // temperature and not set_physiology's Tair baseline.
    ok(l.vpd_leaf_ == l.atm_vpd_ + (l.saturation_vapour_pressure(l.Tleaf_) -
                                    l.saturation_vapour_pressure(l.Tair_)),
       "the leaf-to-air deficit is at the reported Tleaf" + what);
    ok(l.vpd_leaf_ != l.atm_vpd_,
       "and it is not the air deficit -- there is something to get wrong" + what);
  }

  // Off the PM path nothing moves: leaf_temp_ IS the temperature the block was
  // derived at, so there was never a disagreement to fix. Asserted because the
  // change adds a conditional recompute, and the gate-off arm has to stay exactly
  // as it was -- this is the arm the golden file covers.
  {
    Drivers dc = d;
    dc.PPFD = 900.0;
    phylloptim::Leaf l = make_pm_leaf(dc, {17.0}, {1.0}, false);
    l.find_root_collar_psi();
    ok(l.operating_point_kind() == Kind::HydraulicShutdown,
       "the prescribed path shuts down here too");
    ok(l.Tleaf_ == l.leaf_temp_, "and Tleaf is still the driver");
    const double rd = l.R_d_;
    l.update_temperature_dependent_params(l.leaf_temp_);
    ok(l.R_d_ == rd, "R_d is the driver temperature's value, bit-for-bit");
  }
}

// `set_leaf_states_rates_from_psi_leaf`'s zero-transpiration branch is
// SELF-CONTAINED: every quantity it reports is derived from its own `Tleaf_` and
// from nothing the caller left behind.
//
// ⚠️ WHY THE ASSERTIONS ARE ORDER-INDEPENDENCE RATHER THAN VALUES. The transpiring
// branch re-derives the temperature block per candidate by design, so this branch
// runs with the PREVIOUS candidate's block in place. A wrong value here is
// therefore plausible at any single history and only shows up as disagreement
// BETWEEN histories -- and because this function is the objective the collar solve
// maximises, that disagreement is an objective whose value depends on evaluation
// order, which is what hazard 3 forbids. Add a value assertion here if it helps,
// but the one that has teeth is "the same call twice, from two histories, agrees".
//
// ⚠️ `operating_points.tsv` cannot see any of this: the grid runs gate-off, where
// `set_leaf_vpd` returns `atm_vpd_` exactly and the block is never re-derived. A
// bit-identical golden run is not evidence about this test's subject.
void test_zero_E_branch_derives_its_own_block() {
  printf("the zero-E branch inside the objective reads its own temperature\n");
  Drivers d;
  d.leaf_temp = 30.0;  // AIR temperature on this path
  d.PPFD = 900.0;      // full sun, so a leaf that stops transpiring runs hot

  // The same (psi_leaf, psi_upstream) from two histories: one object whose block
  // still holds set_physiology's baseline, one where a transpiring candidate has
  // moved it. Bit-exact, since the branch derives all of these from its own
  // `Tleaf_`.
  {
    phylloptim::Leaf cold = make_pm_leaf(d, {1.0}, {1.0}, true);
    cold.set_leaf_states_rates_from_psi_leaf(1.0, 1.0);

    phylloptim::Leaf warm = make_pm_leaf(d, {1.0}, {1.0}, true);
    warm.set_leaf_states_rates_from_psi_leaf(5.0, 1.0);  // transpiring
    const double T_transpiring = warm.Tleaf_;
    warm.set_leaf_states_rates_from_psi_leaf(1.0, 1.0);  // the SAME zero-E point

    // The premise, asserted rather than assumed: the two histories must really
    // differ, or the comparisons below pass on a block that never moved. Re-check
    // these two if the drivers here are ever retuned.
    ok(T_transpiring != cold.Tleaf_,
       "the transpiring candidate sits at a different temperature");
    ok(cold.Tleaf_ > cold.Tair_,
       "and the zero-E leaf is hotter than the air -- no latent cooling");

    ok(warm.Tleaf_ == cold.Tleaf_, "Tleaf does not depend on the history");
    ok(warm.assim_colimited_ == cold.assim_colimited_,
       "and neither does assimilation");
    ok(warm.R_d_ == cold.R_d_, "nor respiration");
    ok(warm.vpd_leaf_ == cold.vpd_leaf_, "nor the leaf-to-air deficit");
    ok(warm.ci_ == cold.ci_, "nor the internal CO2");
  }

  // Order-independence alone would also hold if every history were wrong the same
  // way, so pin the value too: it is this branch's own temperature's. Same
  // identities the exits outside the solve carry, checked again here because this
  // branch reaches them by a different route.
  {
    phylloptim::Leaf l = make_pm_leaf(d, {1.0}, {1.0}, true);
    l.set_leaf_states_rates_from_psi_leaf(1.0, 1.0);
    const double Tleaf = l.Tleaf_, Tair = l.Tair_;

    ok(l.transpiration_ == 0.0, "the branch really transpires nothing");
    ok(l.Tleaf_ == l.leaf_temp_from_E(0.0),
       "and its temperature is the E = 0 one");
    ok(l.ci_ == l.solve_ci(0.0), "ci sits at the compensation point of this block");
    ok(l.stom_cond_CO2_ == 0.0 && l.assim_colimited_ == 0.0,
       "nothing conducts, so net assimilation is the root's zero exactly");
    // ci is the root of net assimilation, to the root-find's 1e-13 in ci. The
    // bit-exact statement that pins the TEMPERATURE is the next block.
    near(l.assim_colimited(l.ci_), 0.0, 1e-8,
         "and this block's assimilation at that ci is zero");
    ok(l.vpd_leaf_ == l.atm_vpd_ + (l.saturation_vapour_pressure(Tleaf) -
                                    l.saturation_vapour_pressure(Tair)),
       "and the deficit is the leaf-to-air one at that temperature");
    ok(l.vpd_leaf_ != l.atm_vpd_,
       "which is not the air deficit -- there is something to get wrong");

    // Bit-exact against the model's own curves at the reported temperature. A
    // tolerance here would pass on a block derived at some third temperature.
    const double rd = l.R_d_, gamma = l.gamma_, ci = l.ci_;
    l.update_temperature_dependent_params(Tleaf);
    ok(l.R_d_ == rd, "R_d is the curve's value at the reported Tleaf");
    ok(l.gamma_ == gamma, "and so is the compensation point");
    ok(ci == l.solve_ci(0.0), "so ci follows it");
    (void)gamma;
    l.update_temperature_dependent_params(Tair);
    ok(l.R_d_ != rd, "and it is NOT the Tair value -- the two really differ");
  }

  // Closure is still a profit CANDIDATE the optimiser can select (hazard 11) even
  // though it is excluded from |A|max, so its carbon and its temperature have to
  // stay mutually consistent. This used to be asserted through
  // `profitmax_scan_A_[0]`; |A|max is found by optimisation now and those arrays
  // are gone, so the property is asserted where it actually lives.
  {
    phylloptim::Leaf l = make_pm_leaf(d, {1.0}, {1.0}, true);
    l.use_thermal_cost_ = true;
    l.set_leaf_states_rates_from_psi_leaf(1.0, 1.0);      // psi_leaf == psi_soil
    const double a0 = l.assim_colimited_;
    ok(l.transpiration_ == 0.0, "at psi_leaf == psi_soil no water moves");
    l.set_leaf_states_rates_from_psi_leaf(1.0, 1.0);
    ok(l.assim_colimited_ == a0, "and the closure candidate is reproducible");
    l.update_temperature_dependent_params(l.leaf_temp_from_E(0.0));
    // Zero exactly: nothing conducts, so ci is where A is zero, at the
    // temperature E = 0 implies -- which the root itself confirms to its tolerance.
    ok(a0 == 0.0, "its carbon is the root's zero exactly");
    near(l.assim_colimited(l.solve_ci(0.0)), 0.0, 1e-8,
         "and A at the compensation point of the E = 0 temperature is zero");

    // And |A|max EXCLUDES it: a positive normaliser cannot be a respiration rate.
    l.prepare_profitmax();
    ok(l.profitmax_A_max_ > 0.0, "|A|max is a positive carbon gain");
    ok(l.profitmax_A_max_ != std::abs(a0),
       "and it is not the closure point's magnitude");
  }

  // Gate off, which is the arm the golden file covers and which must stay exactly
  // as it is: the recompute is inside `if (use_energy_balance_)`, and
  // `set_leaf_vpd` returns `atm_vpd_` there. This is the guard on that.
  {
    phylloptim::Leaf a = make_pm_leaf(d, {1.0}, {1.0}, false);
    a.set_leaf_states_rates_from_psi_leaf(1.0, 1.0);
    phylloptim::Leaf b = make_pm_leaf(d, {1.0}, {1.0}, false);
    b.set_leaf_states_rates_from_psi_leaf(5.0, 1.0);
    b.set_leaf_states_rates_from_psi_leaf(1.0, 1.0);
    ok(a.Tleaf_ == a.leaf_temp_, "gate off: Tleaf is still the driver");
    ok(a.vpd_leaf_ == a.atm_vpd_, "gate off: the deficit is still the air one");
    ok(b.assim_colimited_ == a.assim_colimited_ && b.ci_ == a.ci_ &&
           b.vpd_leaf_ == a.vpd_leaf_,
       "gate off: the branch was already history-independent");
  }
}

void test_energy_balance_path_runs() {
  printf("Penman-Monteith energy-balance path\n");
  Drivers d;
  phylloptim::Leaf l = make_pm_leaf(d, {2.0}, {1.0}, false);
  l.find_root_collar_psi();
  const double A_prescribed = l.assim_colimited_;

  phylloptim::Leaf eb = make_pm_leaf(d, {2.0}, {1.0}, true);
  eb.find_root_collar_psi();
  ok(std::isfinite(eb.profit_), "energy-balance profit is finite");
  ok(std::isfinite(eb.assim_colimited_), "energy-balance assimilation is finite");
  ok(eb.assim_colimited_ != A_prescribed,
     "energy balance changes the operating point");
  const double Tleaf = eb.leaf_temp_from_E(eb.transpiration_);
  ok(Tleaf > eb.Tair_ - 20.0 && Tleaf < eb.Tair_ + 20.0,
     "leaf temperature is within 20 K of the air");
  // This line used to be the ONLY way to get at the leaf's temperature: derive it
  // again, outside the model, from an output. It is now reported, and the two
  // agreeing bit-for-bit is what says the reported value is the one the solve
  // actually ran at rather than a plausible recomputation.
  ok(eb.Tleaf_ == Tleaf, "the reported Tleaf is the solve's own");

  // The wind model really is what set ra_, rather than the fixed fallback.
  // ⚠️ Assert that rather than assuming it: setting wind_speed_/d_ and then
  // reassigning the leaf discards them, and the assertion still passes because
  // 2.0 / 0.05 are also the defaults.
  ok(std::isfinite(eb.ra_) && eb.ra_ > 0.0, "ra is finite and positive");
  near(eb.ra_, phylloptim::aerodynamic_resistance_coef * std::sqrt(0.05 / 2.0), 1e-12,
       "ra comes from the wind model, not the fixed fallback");
  ok(std::isfinite(eb.Rn_), "net radiation is finite");
}

// Isaac Towers' review of plant #567: the ra fallback accepted a non-finite wind
// speed as well as a zero one. Zero wind is physically ra -> infinity and a
// legitimate fallback; NA is a broken driver or an unset trait and should fail
// rather than silently produce a plausible number. Fixed in plant `76df7169` and
// carried into this package -- but the test lived only in plant, in a demo smoke
// test, while the code now lives here. Ported so the package that owns the
// contract also guards it.
void test_pm_wind_speed_validation() {
  printf("PM path fails fast on a non-finite wind speed (review: itowers1)\n");
  Drivers d;
  const double nan_v = std::numeric_limits<double>::quiet_NaN();

  // 1. gate ON + non-finite wind: fail fast.
  bool threw = false;
  try {
    make_pm_leaf(d, {2.0}, {1.0}, true, nan_v);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "a non-finite wind speed throws on the energy-balance path");

  // ... and the same for the leaf dimension, which enters the same formula.
  threw = false;
  try {
    make_pm_leaf(d, {2.0}, {1.0}, true, 2.0, nan_v);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "a non-finite leaf dimension throws on the energy-balance path");

  // 2. gate OFF + non-finite wind: fine, the wind model is never read.
  threw = false;
  try {
    phylloptim::Leaf off = make_pm_leaf(d, {2.0}, {1.0}, false, nan_v);
    off.find_root_collar_psi();
    ok(std::isfinite(off.profit_), "and still solves");
  } catch (const std::exception &) {
    threw = true;
  }
  ok(!threw, "a non-finite wind speed is ignored with the gate off");

  // 3. gate ON + ZERO wind: legitimate (ra -> infinity), falls back to the fixed
  // ra rather than erroring or producing an infinity.
  threw = false;
  try {
    phylloptim::Leaf zero = make_pm_leaf(d, {2.0}, {1.0}, true, 0.0);
    near(zero.ra_, phylloptim::aerodynamic_resistance_fixed, 1e-12,
         "zero wind falls back to the fixed ra");
    zero.find_root_collar_psi();
    ok(std::isfinite(zero.profit_), "and still solves");
  } catch (const std::exception &) {
    threw = true;
  }
  ok(!threw, "zero wind is a legitimate case, not an error");
}

// Light arrives INTERCEPTED and the leaf absorbs it: net radiation is each band
// times its own absorptance, and the NIR band has no default because the right
// one depends on the light source.
void test_pm_absorbed_shortwave() {
  printf("PM net radiation is the absorbed shortwave of two intercepted bands\n");
  Drivers d;

  phylloptim::Leaf l = make_pm_leaf(d, {2.0}, {1.0}, true);
  // The longwave half is the isothermal sky term of #97, written out from the
  // drivers here so the shortwave half is asserted against something independent
  // of Rn_ itself; test_longwave_and_radiative_conductance owns the longwave.
  const double Tair_K = l.Tair_ + phylloptim::C_to_K;
  const double ea_Pa = 1000.0 * (l.saturation_vapour_pressure(l.Tair_) - l.atm_vpd_);
  // Well below 1 at these drivers, where the soft minimum is exactly the fit.
  const double ema =
      phylloptim::atmos_emissivity_coef *
      std::pow(ea_Pa / Tair_K, phylloptim::atmos_emissivity_exponent);
  const double expected =
      l.alpha_PAR_ * d.PPFD / phylloptim::umol_par_per_joule +
      l.alpha_NIR_ * l.NIR_int_ -
      (1.0 - ema) * phylloptim::stefan_boltzmann * Tair_K * Tair_K * Tair_K *
          Tair_K;
  near(l.Rn_, expected, 1e-12,
       "Rn_ = alpha_PAR*PAR + alpha_NIR*NIR - isothermal longwave");

  // The NIR band heats the leaf and drives no photosynthesis.
  phylloptim::Leaf hot = make_pm_leaf(d, {2.0}, {1.0}, true);
  hot.NIR_int_ = 2.0 * l.NIR_int_;
  hot.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}), d.PPFD,
                     {2.0}, {1.0}, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                     d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  near(hot.Rn_ - l.Rn_, l.alpha_NIR_ * l.NIR_int_, 1e-12,
       "doubling NIR adds alpha_NIR * NIR to Rn_");
  near(hot.electron_transport_, l.electron_transport_, 0.0,
       "and leaves electron transport untouched");

  // No default: a missing NIR band throws on the energy-balance path...
  phylloptim::Leaf unset = fixture::multilayer_leaf();
  unset.setup_transpiration(100);
  unset.setup_root_vulnerability(100);
  unset.use_energy_balance_ = true;
  bool threw = false;
  try {
    unset.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}),
                         d.PPFD, {2.0}, {1.0}, d.K_s * d.theta / d.h,
                         d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  } catch (const std::exception &) {
    threw = true;
  }
  ok(threw, "an unset NIR_int_ throws on the energy-balance path");

  // ...and is never read off it.
  unset.use_energy_balance_ = false;
  threw = false;
  try {
    unset.set_physiology(fixture::root_network({1.0 / d.area_leaf}, {1.0}),
                         d.PPFD, {2.0}, {1.0}, d.K_s * d.theta / d.h,
                         d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
    unset.find_root_collar_psi();
    ok(std::isfinite(unset.profit_), "and solves without it off the path");
  } catch (const std::exception &) {
    threw = true;
  }
  ok(!threw, "an unset NIR_int_ is ignored with the gate off");
}

// The behavioural content of plant's PM demo smoke test, which asserted these
// through the R shim over a Fick-vs-PM grid. The implementation is here now, so
// the sign of the effect is asserted here too.
void test_pm_leaf_temperature_response() {
  printf("PM leaf temperature: Tleaf == Tair off, departs from it on\n");
  for (double ppfd : {400.0, 2000.0}) {
    for (double tair : {20.0, 40.0}) {
      // Both deficits below saturation at 20 C (2.34 kPa).
      for (double vpd : {1.0, 2.0}) {
        Drivers d;
        d.PPFD = ppfd; d.leaf_temp = tair; d.atm_vpd = vpd;
        const std::string at = " at PPFD=" + std::to_string(int(ppfd)) +
                               " Tair=" + std::to_string(int(tair)) +
                               " VPD=" + std::to_string(int(vpd));

        phylloptim::Leaf fick = make_pm_leaf(d, {2.0}, {1.0}, false);
        fick.find_root_collar_psi();
        ok(std::isfinite(fick.profit_) && std::isfinite(fick.assim_colimited_),
           "Fick outputs are finite" + at);

        phylloptim::Leaf pm = make_pm_leaf(d, {2.0}, {1.0}, true);
        pm.find_root_collar_psi();
        ok(std::isfinite(pm.profit_) && std::isfinite(pm.assim_colimited_),
           "PM outputs are finite" + at);

        // Hot and bright is where PM matters: the leaf runs warmer than the air.
        if (ppfd == 2000.0 && tair == 40.0) {
          const double Tleaf = pm.leaf_temp_from_E(pm.transpiration_);
          ok(Tleaf > pm.Tair_, "the leaf is warmer than the air when hot and bright");
        }
      }
    }
  }
}

// With the energy balance on, the collar solve finds the maximum of the objective
// it reports: interior optima satisfy dprofit == 0 to solver precision, and a
// derivative-free scan of profit, which cannot share a defect of the derivative,
// finds the same collar to its own resolution and no more profit.
void test_energy_balance_collar_solve_finds_the_maximum() {
  printf("energy balance: the collar solve against a derivative-free scan\n");

  double worst_dpsi = 0.0, worst_resid = 0.0, worst_dprofit = 0.0,
         worst_pinned_shortfall = 0.0;
  int rows = 0, pinned_rows = 0;

  for (double ppfd : {500.0, 1500.0}) {
    for (double tair : {20.0, 30.0, 40.0}) {
      for (double vpd : {1.0, 2.0}) {
        Drivers d;
        d.PPFD = ppfd; d.leaf_temp = tair; d.atm_vpd = vpd;

        phylloptim::Leaf eb = make_pm_leaf(d, {2.0}, {1.0}, true);
        eb.find_root_collar_psi();
        const double psi_solver = eb.opt_root_psi_;
        const double profit_solver = eb.profit_;

        // The residual the solver believes it drove to zero. On the non-EB path
        // this is ~5.6e-15; here it is whatever staleness leaves.
        bool feasible = true;
        const double resid = eb.dprofit_droot_collar_psi(psi_solver, &feasible);
        if (!feasible) continue;

        // ORACLE: scan profit over the same feasible interval. A fresh leaf,
        // because dprofit_droot_collar_psi above has just re-seated the
        // temperature parameters at yet another point.
        phylloptim::Leaf scan = make_pm_leaf(d, {2.0}, {1.0}, true);
        double lo = 0.0, hi = 0.0;
        if (!scan.prepare_collar_solve<phylloptim::Leaf::CostCurve::TF24>(lo, hi)) continue;
        const int N = 20001;
        double best_psi = lo, best_profit = -std::numeric_limits<double>::max();
        for (int i = 0; i < N; ++i) {
          const double psi = lo + (hi - lo) * double(i) / double(N - 1);
          const double p = scan.profit_at_collar_psi<phylloptim::Leaf::CostCurve::TF24>(psi, lo, hi);
          if (p > best_profit) { best_profit = p; best_psi = psi; }
        }

        const double dpsi = std::abs(best_psi - psi_solver);
        const double dprofit = best_profit - profit_solver;
        worst_dpsi = std::max(worst_dpsi, dpsi);

        // ⚠️ INTERIOR ROWS ONLY for the residual. A CONSTRAINED optimum sits on
        // a bracket bound, where dprofit is genuinely non-zero and "the gradient
        // should vanish" is simply the wrong statement -- the same distinction
        // test_collar_optimum_pinned_to_its_constraint makes on the non-EB path.
        // Lumping the two together reports a pinned row's honest gradient as if
        // it were a solver failure.
        //
        // ⚠️ AND INTERIOR ROWS ONLY for the shortfall, at g_min = 0. A wet-pinned
        // leaf is reported a 1e-6 step inside the zero-uptake collar, because
        // the bound itself has no conductance; the scan evaluates the bound,
        // where a net source nets its saturated assimilation. The gap between
        // them is the boundary-soil artefact #183 records (measured 9.2e-05
        // here), which a residual conductance removes and this grid, at
        // g_min = 0, keeps. It is printed, not asserted away.
        const double span = std::max(hi - lo, 1e-12);
        const bool pinned = (psi_solver - lo) < 1e-3 * span ||
                            (hi - psi_solver) < 1e-3 * span;
        if (pinned) {
          ++pinned_rows;
          worst_pinned_shortfall = std::max(worst_pinned_shortfall, dprofit);
        } else {
          worst_resid = std::max(worst_resid, std::abs(resid));
          worst_dprofit = std::max(worst_dprofit, dprofit);
        }
        ++rows;
      }
    }
  }

  printf("    %d feasible rows (%d pinned) | worst |dprofit|, interior only: %.3e\n",
         rows, pinned_rows, worst_resid);
  printf("    worst |psi_scan - psi_solver|: %.3e MPa | worst profit shortfall: %.3e "
         "(pinned rows, the step-in artefact at g_min = 0: %.3e)\n",
         worst_dpsi, worst_dprofit, worst_pinned_shortfall);

  ok(rows > 0, "the EB grid has feasible rows to measure");

  ok(worst_resid < 1e-9,
     "interior EB optima satisfy dprofit == 0 to solver precision");
  // The scan resolves the argmax only to (hi-lo)/20000, ~5e-5 MPa on this grid,
  // so agreement AT that scale is agreement -- asking for more would be asking
  // the oracle for precision it does not have.
  ok(worst_dpsi < 1e-3,
     "the collar agrees with a derivative-free scan to the scan's resolution");
  // The scan can only ever find a profit >= the solver's, up to its own grid
  // resolution. A NEGATIVE shortfall beyond that would mean the scan is broken,
  // which is the one way this measurement could mislead.
  ok(worst_dprofit > -1e-6, "the scan never finds LESS profit than the solver");
  ok(worst_dprofit < 1e-6, "and never finds MORE: the solver is at the maximum");
}

// Gate-off inertness, proved by POISONING rather than by a recorded value. The
// golden file already pins values; what it cannot say is that the energy-balance
// inputs are never *read* when the gate is off. Setting them to NaN and demanding
// bit-identical results says exactly that: if any of the new code path touched
// Rn_, g_rad_, ra_ or Tair_ off-gate, a NaN would propagate and nothing would
// compare equal. Same shape as test_pm_wind_speed_validation uses for wind_speed_.
void test_energy_balance_gate_off_is_inert() {
  printf("gate off: the EB inputs are never read (NaN poisoning)\n");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (double tair : {20.0, 40.0}) {
    Drivers d;
    d.leaf_temp = tair;
    phylloptim::Leaf clean = make_pm_leaf(d, {2.0}, {1.0}, false);
    phylloptim::Leaf poisoned = make_pm_leaf(d, {2.0}, {1.0}, false);
    poisoned.Rn_ = nan; poisoned.g_rad_ = nan; poisoned.ra_ = nan;
    poisoned.Tair_ = nan;

    clean.find_root_collar_psi();
    poisoned.find_root_collar_psi();
    const std::string at = " at Tair=" + std::to_string(int(tair));
    ok(clean.profit_ == poisoned.profit_, "profit is bit-identical" + at);
    ok(clean.assim_colimited_ == poisoned.assim_colimited_,
       "assimilation is bit-identical" + at);
    ok(clean.opt_root_psi_ == poisoned.opt_root_psi_,
       "the collar is bit-identical" + at);
    // And the derivative, which is where the new code lives.
    const double psi = clean.opt_root_psi_;
    ok(clean.dprofit_droot_collar_psi(psi) ==
           poisoned.dprofit_droot_collar_psi(psi),
       "dprofit is bit-identical" + at);
  }
}

// THE SCIENTIFIC ACCEPTANCE TEST. Jones et al. (2026, Global Change Biology
// 32:e70972) report that stomatal conductance keeps rising above the
// photosynthetic thermal optimum while assimilation falls, and that an
// optimality model reproduces it only when leaf temperature is inside the
// objective. This asserts that the corrected first-order condition produces that
// behaviour, and -- the half that matters -- that the gate-off arm does NOT, at
// the same drivers. Without the contrast the test could pass on the Arrhenius
// optimum alone, which has nothing to do with energy balance.
//
// ⚠️ THE SIGNATURE IS IN TRANSPIRATION, NOT IN CONDUCTANCE. Conductance rising
// where assimilation falls is what you see when stom_cond_CO2 divides by the
// prescribed AIR deficit however hot the leaf gets: gs is then E rescaled by a
// constant and inherits E's shape exactly. With the leaf-to-air deficit wired in
// the two come apart, because
//
//     gs = P*E / (1.6 * D_leaf),      D_leaf = atm_vpd + esat(Tleaf) - esat(Tair)
//
// and over the decoupled window D_leaf grows FASTER than E does. Measured on this
// sweep, between the thermal optimum and E's own peak: E x1.034, D_leaf x1.089,
// so gs x0.950 -- falling throughout. The gs counts went from a positive number
// to gate-on 0, gate-off 1.
//
// That is not a loss of the mechanism, and the assertions below are rewritten to
// say which is which:
//
//   * the ENERGY-BALANCE mechanism is intact and is what the gate buys. E keeps
//     rising above the thermal optimum, and it does so over a wider window with
//     the gate on (4 points) than off (1).
//   * the CONDUCTANCE signature, which is what a gas-exchange dataset reports,
//     does not survive the correction at these drivers. Anyone comparing this
//     model to a measured dgs/dT needs to know that, and the third assertion pins
//     the ratio that causes it rather than the count it produces.
//
// The sweep holds AIR VPD fixed, which is the protocol of the paper's Figure S2
// rather than its Figure 2 -- and note that fixing the air deficit no longer
// fixes the deficit the leaf actually sees.
void test_energy_balance_stomatal_decoupling() {
  printf("decoupling: E rises where A falls; gs does not, once D moves with Tleaf\n");
  // ⚠️ THE FIXTURE DEFAULTS CANNOT SHOW THIS, and the reason is physical rather
  // than numerical. kmax = K_s*theta/h is 3.14e-5 at the defaults, and a leaf
  // that conductive simply cannot move enough water to cool itself: measured in
  // R against Jones's own drivers, Tleaf - Tair floors at +1.9 K where their
  // leaf reaches -2.0 K, so the evaporative-cooling benefit never becomes large
  // enough to reverse the conductance response. Shortening the path length to
  // 1.5 m (kmax ~ 1.05e-4) puts the leaf in the range their measured species
  // occupy. This is choosing a regime where the mechanism is active, not tuning
  // until a test passes -- at the defaults the honest answer is "no decoupling",
  // and that is what the probe found.
  auto sweep = [&](bool gate) {
    std::vector<double> Ta, A, gs, E, D;
    for (double t = 15.0; t <= 45.0; t += 1.0) {
      Drivers d;
      d.PPFD = 1200.0; d.leaf_temp = t; d.atm_vpd = 1.0; d.h = 1.5;
      phylloptim::Leaf l = make_pm_leaf(d, {0.5}, {1.0}, gate);
      // ⚠️ OVERRIDE THE DERIVED RADIATION BUDGET, and this is the second thing
      // the defaults cannot express. At PPFD 1200 the derived isothermal Rn is
      // 420-525 W m^-2 over this sweep, and with ra from the wind model (31.6 s
      // m^-1) the leaf runs 10-12 K above air at zero transpiration, so its
      // thermal optimum falls off the bottom of any sensible sweep. These are
      // the values Jones et al. drive with (their eq. 25 and 29 at Iabs = 500,
      // ra = 10), and the fields are settable precisely so a caller can supply
      // its own budget. ⚠️ ra_ = 12 is their r_HR, which ALREADY folds in the
      // radiative conductance, so g_rad_ is zeroed or that path counts twice.
      l.Rn_ = 400.0;
      l.ra_ = 12.0;
      l.g_rad_ = 0.0;
      l.find_root_collar_psi();
      if (!std::isfinite(l.assim_colimited_) || !std::isfinite(l.stom_cond_CO2_)) continue;
      Ta.push_back(t); A.push_back(l.assim_colimited_); gs.push_back(l.stom_cond_CO2_);
      E.push_back(l.transpiration_); D.push_back(l.vpd_leaf_);
    }
    return std::make_tuple(Ta, A, gs, E, D);
  };

  // Points above the thermal optimum where the response variable RISES while
  // assimilation falls. `y` is E for the mechanism and gs for the signature.
  auto count_decoupled = [](const std::vector<double>& A,
                            const std::vector<double>& y) {
    if (A.size() < 5) return 0;
    const size_t peak = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
    int n = 0;
    for (size_t i = peak + 2; i + 1 < A.size(); ++i) {
      if ((A[i + 1] - A[i - 1]) < 0.0 && (y[i + 1] - y[i - 1]) > 0.0) ++n;
    }
    return n;
  };

  const auto on = sweep(true);
  const auto off = sweep(false);
  const auto& A_on = std::get<1>(on);
  const auto& gs_on = std::get<2>(on);
  const auto& E_on = std::get<3>(on);
  const auto& D_on = std::get<4>(on);

  const int nE_on = count_decoupled(A_on, E_on);
  const int nE_off = count_decoupled(std::get<1>(off), std::get<3>(off));
  const int ng_on = count_decoupled(A_on, gs_on);
  const int ng_off = count_decoupled(std::get<1>(off), std::get<2>(off));
  printf("    decoupled points -- in E: gate on %d, gate off %d;"
         " in gs: gate on %d, gate off %d\n", nE_on, nE_off, ng_on, ng_off);

  ok(A_on.size() > 5, "the decoupling sweep produced a curve");
  const size_t peak = std::distance(A_on.begin(),
                                    std::max_element(A_on.begin(), A_on.end()));
  ok(peak > 0 && peak + 1 < A_on.size(),
     "assimilation has an interior thermal optimum");
  ok(nE_on > 0, "with the gate ON, transpiration rises where assimilation falls");
  ok(nE_on > nE_off,
     "and it does so over a wider window than with the gate off -- "
     "the energy balance is why");

  // WHY the conductance signature does not follow, stated as the ratio rather
  // than as the count it happens to produce here. Over the window from the
  // thermal optimum to E's own peak, the deficit outgrows the flux, and
  // gs = P*E/(1.6*D) therefore falls even though E is rising.
  const size_t Epeak = std::distance(E_on.begin(),
                                     std::max_element(E_on.begin(), E_on.end()));
  ok(Epeak > peak, "transpiration peaks ABOVE the assimilation optimum");
  const double E_growth = E_on[Epeak] / E_on[peak];
  const double D_growth = D_on[Epeak] / D_on[peak];
  printf("    optimum -> E peak: E x%.4f, D_leaf x%.4f, gs x%.4f\n",
         E_growth, D_growth, gs_on[Epeak] / gs_on[peak]);
  ok(D_growth > E_growth,
     "the leaf-to-air deficit outgrows the flux, so gs falls where E rises");
}

// The boundary layer's split of the CO2 path (#122): the plain series wherever
// the total conductance is below the boundary layer's own, which is every
// physical point; finite and rising past it, where only a candidate on
// vpd_leaf_min goes; and with the slope the hand-written partials are scaled by.
void test_boundary_layer_split_of_the_co2_path() {
  printf("the boundary layer's split of the CO2 path (#122)\n");
  Drivers d;
  phylloptim::Leaf l = make_pm_leaf(d, {2.0}, {1.0}, true);
  const double k = l.boundary_layer_co2_excess();
  const double g_bw = l.boundary_layer_cond_H2O();
  const double r_s = l.H2O_CO2_stom_diff_ratio_;
  near(k, (r_s - std::pow(r_s, 2.0 / 3.0)) / g_bw, 1e-14,
       "k is (r_s - r_s^(2/3)) / g_bw");

  // Every physical point: total water conductance from zero to g_bw.
  int exact = 0, rows = 0;
  double worst_slope = 0.0;
  for (double frac : {0.0, 1e-3, 0.05, 0.25, 0.5, 0.9, 0.999}) {
    const double gc_hat = frac * g_bw / r_s;
    const auto split = l.boundary_layer_split_CO2(gc_hat);
    exact += split.first == gc_hat / (1.0 - k * gc_hat);
    const double h = 1e-6 * std::max(gc_hat, 1e-3);
    const double fd = (l.boundary_layer_split_CO2(gc_hat + h).first -
                       l.boundary_layer_split_CO2(gc_hat - h).first) /
                      (2.0 * h);
    worst_slope = std::max(worst_slope, std::abs(split.second / fd - 1.0));
    if (gc_hat > 0.0) {
      near(l.boundary_layer_split_slope_at(split.first), split.second, 1e-12,
           "the slope at a solved conductance is the split's own");
    }
    ++rows;
  }
  ok(exact == rows, "at every physical point the split is the plain series, "
                    "bit for bit (" + std::to_string(exact) + " of " +
                    std::to_string(rows) + ")");
  ok(worst_slope < 1e-7, "and its slope is the derivative of its value");

  // Past what the boundary layer can carry the plain series has a pole at
  // gc_hat = 1/k; the floor keeps the split finite, positive and rising there.
  double prev = 0.0;
  bool rising = true, finite = true;
  for (double mult : {0.5, 0.9, 0.99, 1.0, 1.01, 2.0, 10.0, 1e3}) {
    const double g = l.boundary_layer_split_CO2(mult / k).first;
    finite = finite && std::isfinite(g) && g > 0.0;
    rising = rising && g > prev;
    prev = g;
  }
  ok(finite, "past the boundary layer's capacity the split stays finite");
  ok(rising, "and rises with the total conductance, through the old pole");
  const double fd_floor =
      (l.boundary_layer_split_CO2(1.0 / k * (1.0 + 1e-6)).first -
       l.boundary_layer_split_CO2(1.0 / k * (1.0 - 1e-6)).first) /
      (2e-6 / k);
  near(l.boundary_layer_split_CO2(1.0 / k).second / fd_floor, 1.0, 1e-6,
       "and its slope is the derivative of its value on the floor too");
  printf("    %d physical rows, worst slope error %.2e\n", rows, worst_slope);
}

// The reported conductance is a TOTAL leaf-to-air one on the energy-balance path,
// and #122's two accessors separate the stoma from the boundary layer in series.
// Off the path there is no boundary layer, so the stomatal conductance must be the
// reported one bit for bit; on it, the two pieces must recombine into the total.
void test_stomatal_conductance_excludes_the_boundary_layer() {
  printf("stomatal conductance with the boundary layer taken out in series (#122)\n");
  Drivers d;
  d.PPFD = 1500.0;
  d.leaf_temp = 30.0;

  phylloptim::Leaf off = make_pm_leaf(d, {2.0}, {1.0}, false);
  off.find_root_collar_psi();
  ok(off.stom_cond_stomatal_H2O() == off.stom_cond_H2O(),
     "gate off: the stomatal conductance IS the reported one, bit for bit");
  ok(std::isinf(off.boundary_layer_cond_H2O()),
     "gate off: the boundary layer is not modelled, so its conductance is Inf");

  // Low wind, low deficit, low light and wet soil, so the boundary layer is a
  // quarter of the resistance and the series correction is far from negligible.
  // The worst row of the grid below.
  Drivers wet = d;
  wet.atm_vpd = 0.5;
  wet.PPFD = 500.0;
  wet.leaf_temp = 20.0;
  phylloptim::Leaf on = make_pm_leaf(wet, {0.5}, {1.0}, true, 0.5);
  on.find_root_collar_psi();
  const double g = on.stom_cond_H2O();
  const double gb = on.boundary_layer_cond_H2O();
  const double gs = on.stom_cond_stomatal_H2O();
  ok(g / gb > 0.1, "gate on: the boundary layer is a material share (hazard 16)");
  // The boundary layer from its definition: the heat resistance ra_ in molar
  // form at air temperature, scaled to vapour.
  const double rho_mol = on.atm_kpa_ * phylloptim::kPa_to_Pa /
                         (phylloptim::gas_constant *
                          (on.Tair_ + phylloptim::C_to_K));
  near(gb, phylloptim::bl_vapour_heat_ratio * rho_mol / on.ra_, 1e-15,
       "g_b is 1.075 * P/(R*Tair) / ra_");
  near(1.0 / gs + 1.0 / gb, 1.0 / g, 1e-14,
       "the stoma and the boundary layer recombine in series to the total");
  ok(gs > g, "the stomatal conductance exceeds the total, as a series part must");
  printf("    g_total %.4f, g_b %.4f, g_s %.4f mol m^-2 s^-1 (g_total/g_b %.3f)\n",
         g, gb, gs, g / gb);

  // A boundary layer too thin to carry the solved E: no finite g_s exists, and
  // the accessor says why rather than returning a negative or infinite number.
  phylloptim::Leaf thin = on;
  thin.ra_ = 1e3 * on.ra_;
  ok(throws_with([&] { thin.stom_cond_stomatal_H2O(); },
                 "no stomatal conductance exists"),
     "a total conductance at or above g_b is refused, with the reason");

  // Shut: no conductance on either side of the series.
  phylloptim::Leaf shut = make_pm_leaf(d, {17.0}, {1.0}, true);
  shut.find_root_collar_psi();
  ok(shut.stom_cond_H2O() == 0.0 && shut.stom_cond_stomatal_H2O() == 0.0,
     "a shut leaf reports zero stomatal conductance, not a refusal");
}

// How far a realistic driver grid stays from g_total = g_b, which is where the
// accessor above refuses. Hydraulics set E, so nothing in the solve keeps E below
// what the boundary layer can carry; this measures the margin instead.
void test_boundary_layer_bounds_the_total_conductance() {
  printf("MEASUREMENT: g_total / g_b over a driver grid, energy balance on\n");
  double worst = 0.0, worst_at_default_wind = 0.0;
  int rows = 0;
  for (double wind : {0.5, 2.0, 5.0})
    for (double tair : {20.0, 30.0, 40.0})
      for (double vpd : {0.5, 1.0, 2.0, 4.0})
        for (double ppfd : {500.0, 1500.0})
          for (double psi_soil : {0.5, 1.0, 3.0}) {
            // A deficit past saturation (4 kPa at 20 C) is refused, not a row.
            if (!(vpd < phylloptim::Leaf().saturation_vapour_pressure(tair))) {
              continue;
            }
            Drivers d;
            d.PPFD = ppfd;
            d.leaf_temp = tair;
            d.atm_vpd = vpd;
            phylloptim::Leaf l = make_pm_leaf(d, {psi_soil}, {1.0}, true, wind);
            l.find_root_collar_psi();
            const double g = l.stom_cond_H2O();
            if (!(g > 0.0)) {
              continue;
            }
            const double share = g / l.boundary_layer_cond_H2O();
            worst = std::max(worst, share);
            if (wind == 2.0) {
              worst_at_default_wind = std::max(worst_at_default_wind, share);
            }
            ok(std::isfinite(l.stom_cond_stomatal_H2O()),
               "a finite stomatal conductance at wind " + std::to_string(wind) +
                   " Tair " + std::to_string(tair) + " vpd " +
                   std::to_string(vpd));
            ++rows;
          }
  printf("    %d transpiring rows | worst g_total/g_b %.3f (at wind 2 m/s: %.3f)\n",
         rows, worst, worst_at_default_wind);
  ok(rows >= 150, "enough transpiring rows were measured (hazard 16)");
  ok(worst < 1.0, "the total conductance stays below the boundary layer's");
}

// ===========================================================================
// Leaf-to-air VPD (#7)
// ---------------------------------------------------------------------------
void test_leaf_to_air_vpd() {
  printf("leaf-to-air VPD: the deficit Fick's law divides by\n");
  Drivers d;
  d.atm_vpd = 1.5;

  // Gate off, the leaf IS at air temperature, so the two deficits must agree --
  // and EXACTLY, not nearly. That is what keeps the prescribed-temperature path,
  // and therefore the golden file, untouched by this change.
  phylloptim::Leaf off = make_pm_leaf(d, {2.0}, {1.0}, false);
  ok(off.vpd_leaf_ == off.atm_vpd_,
     "gate off: vpd_leaf_ is bit-identical to the driver");
  off.find_root_collar_psi();
  ok(off.vpd_leaf_ == off.atm_vpd_,
     "gate off: and still is after a solve");

  // Gate on, with the leaf running hot: the deficit it sees is larger than the
  // air's, so the same water flux implies a SMALLER conductance.
  phylloptim::Leaf on = make_pm_leaf(d, {2.0}, {1.0}, true);
  on.Rn_ = 400.0;
  on.ra_ = 12.0;
  on.find_root_collar_psi();
  const double Tleaf = on.leaf_temp_from_E(on.transpiration_);
  ok(Tleaf > on.Tair_, "the gate-on leaf is hotter than the air here");
  ok(on.vpd_leaf_ > on.atm_vpd_,
     "gate on: a hotter leaf sees a larger deficit than the air's");
  // And the value is the definition, not an approximation of it.
  near(on.vpd_leaf_,
       on.atm_vpd_ + on.saturation_vapour_pressure(Tleaf) -
           on.saturation_vapour_pressure(on.Tair_),
       1e-12, "vpd_leaf_ = atm_vpd + esat(Tleaf) - esat(Tair)");
  printf("    Tair %.1f C, Tleaf %.2f C: D_air %.3f kPa -> D_leaf %.3f kPa (x%.2f)\n",
         on.Tair_, Tleaf, on.atm_vpd_, on.vpd_leaf_, on.vpd_leaf_ / on.atm_vpd_);

  // The consequence, stated as an identity rather than a direction: the total
  // conductance to water is E rescaled by the deficit, so getting the deficit
  // wrong scales it by the ratio.
  const double g_w = on.atm_kpa_ * on.transpiration_ *
                     phylloptim::kg_to_mol_h2o / on.vpd_leaf_;
  near(on.stom_cond_H2O(), g_w, 1e-12,
       "the total conductance is the transpiration divided by the deficit it sees");
  // And to CO2 the boundary layer's share of that series diffuses at
  // r_s^(2/3), not r_s (#122): 1/g_c = r_s/g_w - (r_s - r_s^(2/3))/g_bw.
  const double r_s = phylloptim::H2O_CO2_stom_diff_ratio;
  const double g_bw = on.boundary_layer_cond_H2O();
  near(1.0 / on.stom_cond_CO2_,
       r_s / g_w - (r_s - std::pow(r_s, 2.0 / 3.0)) / g_bw, 1e-12,
       "1/g_c is the series of the stomatal and boundary-layer resistances");
  ok(on.stom_cond_CO2_ > g_w / r_s,
     "which conducts more CO2 than the whole path at the stomatal ratio");
}

// The objective whose derivative dprofit_dpsi_leaf<K> IS: the log of the product
// for the two log links, the reported profit for everything else. NaN where the
// log form does not exist, which the callers skip.
template <phylloptim::Leaf::CostCurve K>
double linked_objective(phylloptim::Leaf &l, double psi_leaf, double psi_soil) {
  const double v = l.template profit_psi_leaf_for<K>(psi_leaf, psi_soil);
  if constexpr (phylloptim::Leaf::benefit_link<K>() ==
                phylloptim::Leaf::BenefitLink::Log) {
    return v > 0.0 ? std::log(v) : std::numeric_limits<double>::quiet_NaN();
  } else {
    return v;
  }
}

// Seeds ProfitMax's normalisers, which its derivative reads through h'; a no-op
// on every other curve. False where the curve has no objective to differentiate.
template <phylloptim::Leaf::CostCurve K>
bool seed_linked_curve(phylloptim::Leaf &l) {
  if constexpr (K == phylloptim::Leaf::CostCurve::ProfitMax) {
    l.template optimise_psi_leaf_single<K>();
    return l.profitmax_A_max() > 0.0 && l.profitmax_k_span_ > 0.0;
  } else {
    return true;
  }
}

template <phylloptim::Leaf::CostCurve K>
int check_linked_eb_derivative_fd(const Drivers &d, bool tc, double psi_soil,
                                  double &worst, int &thermal_rows,
                                  int &assigned_rows, bool print_v) {
  const std::string name =
      phylloptim::Leaf::curve_name(static_cast<int>(K)) + (tc ? "+TC" : "") +
      " at Tair " + std::to_string(static_cast<int>(d.leaf_temp));
  phylloptim::Leaf l = make_linked_eb_leaf(d, psi_soil, tc);
  if (!seed_linked_curve<K>(l)) {
    return 0;
  }
  int rows = 0;
  // The error is read against the derivative's own scale across the bracket,
  // not against its value at the point: near the optimum dJ/dpsi passes through
  // zero, and a relative error there grows without the absolute one moving.
  const double fracs[] = {0.15, 0.35, 0.55, 0.75};
  double scale = 0.0;
  for (double frac : fracs) {
    bool f = false;
    const double v = l.template dprofit_dpsi_leaf<K>(
        psi_soil + frac * (l.psi_crit - psi_soil), psi_soil, &f);
    if (std::isfinite(v)) {
      scale = std::max(scale, std::abs(v));
    }
  }
  for (double frac : fracs) {
    const double psi = psi_soil + frac * (l.psi_crit - psi_soil);
    bool feasible = false;
    const double an = l.template dprofit_dpsi_leaf<K>(psi, psi_soil, &feasible);
    // Read before the objective below re-solves ci: a net source (A < 0, ci
    // above ambient) is a row the old assigned-ci branch used to take, and is
    // counted so the sweep is known to reach it. A log link has no derivative
    // there: that is the 0.0 sentinel, and it is not comparable.
    const bool source = l.ci_ > l.ca_;
    if (!std::isfinite(an) || !feasible) {
      continue;
    }
    if (source) {
      ++assigned_rows;
    }
    // ⚠️ THE STEP IS SWEPT, NOT PICKED (hazard 15). The differenced objective
    // runs through the ci root-find AND the leaf-temperature balance, so its
    // noise floor is not the gate-off test's, and only the floor of the V is a
    // derivative. The V is printed for one configuration so the floor is read
    // rather than assumed.
    double best = std::numeric_limits<double>::infinity();
    std::string v;
    for (double h : {1e-5, 1e-4, 1e-3, 1e-2}) {
      const double fd = (linked_objective<K>(l, psi + h, psi_soil) -
                         linked_objective<K>(l, psi - h, psi_soil)) /
                        (2.0 * h);
      if (!std::isfinite(fd)) {
        continue;
      }
      const double err = std::abs(an - fd) / scale;
      best = std::min(best, err);
      char buf[48];
      snprintf(buf, sizeof buf, "  %.0e: %.1e", h, err);
      v += buf;
    }
    if (!std::isfinite(best)) {
      continue;
    }
    if (print_v) {
      printf("    %s, frac %.2f, error/scale by step:%s\n", name.c_str(), frac,
             v.c_str());
    }
    ok(best < 2e-4, name + " dprofit at frac=" + std::to_string(frac) +
                        " matches a difference at the floor of its V");
    worst = std::max(worst, best);
    // Tleaf_ is written by every objective evaluation, so this reads the last
    // point differenced; a live thermal cost there means the term was exercised.
    if (tc && l.thermal_cost_at(l.Tleaf_) > 1e-2) {
      ++thermal_rows;
    }
    ++rows;
  }
  return rows;
}

void test_linked_dprofit_with_energy_balance_matches_a_finite_difference() {
  printf("linked curves, energy balance on: dprofit/dpsi_leaf against a "
         "swept-step difference of the objective\n");
  using CC = phylloptim::Leaf::CostCurve;
  double worst = 0.0;
  int rows = 0, thermal_rows = 0, assigned_rows = 0;
  // 25 and 30 C air give interior rows; at 45 C every row under these drivers
  // has ci assigned at ca, which is the other branch the link reaches.
  for (double tair : {25.0, 30.0, 45.0}) {
    Drivers d;
    d.PPFD = 1500.0;
    d.leaf_temp = tair;
    for (bool tc : {false, true}) {
      const bool print_v = (tair == 25.0 && tc);
      rows += check_linked_eb_derivative_fd<CC::SOX>(
          d, tc, 1.0, worst, thermal_rows, assigned_rows, print_v);
      rows += check_linked_eb_derivative_fd<CC::JW26>(
          d, tc, 1.0, worst, thermal_rows, assigned_rows, false);
      rows += check_linked_eb_derivative_fd<CC::ProfitMax>(
          d, tc, 1.0, worst, thermal_rows, assigned_rows, print_v);
      rows += check_linked_eb_derivative_fd<CC::LeastCost>(
          d, tc, 1.0, worst, thermal_rows, assigned_rows, print_v);
    }
  }
  printf("    %d rows (%d on the compensation branch) | worst error/scale at "
         "the floor %.3e | %d rows with a live thermal cost\n",
         rows, assigned_rows, worst, thermal_rows);
  ok(rows - assigned_rows >= 36, "enough interior rows were compared");
  ok(assigned_rows >= 4,
     "the compensation-point branch was compared too (hazard 16)");
  ok(thermal_rows >= 12,
     "the thermal-cost term was exercised, not only compiled (hazard 16)");
}

void test_linked_dprofit_with_energy_balance_is_the_identity_link_scaled() {
  printf("linked curves, energy balance on: dJ/dpsi is h'(A) dA/dpsi - dC/dpsi, "
         "dA/dpsi from a zero-cost curve\n");
  using CC = phylloptim::Leaf::CostCurve;
  const double psi_soil = 1.0;
  int rows = 0;
  double worst = 0.0;
  for (double tair : {25.0, 30.0}) {
    Drivers d;
    d.PPFD = 1500.0;
    d.leaf_temp = tair;
    for (bool tc : {false, true}) {
      phylloptim::Leaf l = make_linked_eb_leaf(d, psi_soil, tc);
      if (!seed_linked_curve<CC::ProfitMax>(l)) {
        continue;
      }
      // A zero unit cost: dC/dpsi is exactly 0, so JS22's derivative IS dA/dpsi,
      // the temperature path included -- the vignette's link-identity device.
      l.JS22_gamma = 0.0;
      for (double frac : {0.15, 0.45, 0.75}) {
        const double psi = psi_soil + frac * (l.psi_crit - psi_soil);
        bool fa = false;
        const double dA = l.dprofit_dpsi_leaf<CC::JS22>(psi, psi_soil, &fa);
        if (!fa) {
          continue;
        }
        // Seat the operating point the derivative was taken at: A, E and Tleaf.
        l.profit_psi_leaf_for<CC::JS22>(psi, psi_soil);
        const double A = l.assim_colimited_;
        if (!(A > 0.0)) {
          continue;
        }
        // The thermal cost's own slope along psi, for ProfitMax's cost side.
        double dT_dE = 0.0;
        l.leaf_temp_from_E(l.transpiration_, &dT_dE);
        const double dTC_dpsi = l.thermal_cost_temp_deriv(l.Tleaf_) * dT_dE *
                                l.leaf_specific_conductance_max_ *
                                l.stem_curve_integral_deriv(psi);
        auto check = [&](auto tag, double h_prime, double dC, const char *nm) {
          constexpr CC K = decltype(tag)::value;
          bool f = false;
          const double lhs = l.template dprofit_dpsi_leaf<K>(psi, psi_soil, &f);
          if (!f) {
            return;
          }
          const double rhs = h_prime * dA - dC;
          // Against the scale of the two terms rather than their sum, which
          // cancels at the optimum; the points are away from it regardless.
          const double scale = std::max(std::abs(h_prime * dA), std::abs(dC));
          const double err = std::abs(lhs - rhs) / scale;
          worst = std::max(worst, err);
          ok(err <= 1e-10, std::string(nm) + (tc ? "+TC" : "") + " at Tair " +
                               std::to_string(static_cast<int>(tair)) +
                               " frac " + std::to_string(frac) +
                               ": dJ/dpsi is h'(A) dA/dpsi - dC/dpsi");
          ++rows;
        };
        check(std::integral_constant<CC, CC::SOX>{}, 1.0 / A,
              l.cost_deriv<CC::SOX>(psi, psi_soil), "SOX");
        check(std::integral_constant<CC, CC::JW26>{}, 1.0 / A,
              l.cost_deriv<CC::JW26>(psi, psi_soil), "JW26");
        check(std::integral_constant<CC, CC::ProfitMax>{},
              1.0 / l.profitmax_A_max(),
              l.cost_deriv<CC::ProfitMax>(psi, psi_soil) + dTC_dpsi,
              "ProfitMax");
        // LeastCost's capacity term moves with Tleaf through Vcmax:
        // d(log W)/dpsi gains beta*Vcmax_mol'(T)*dT/dpsi / W, its slope taken
        // by vcmax_temp_deriv rather than by the term under test.
        const double dVcmax_dpsi = l.vcmax_temp_deriv(l.Tleaf_) * dT_dE *
                                   l.leaf_specific_conductance_max_ *
                                   l.stem_curve_integral_deriv(psi);
        check(std::integral_constant<CC, CC::LeastCost>{}, 1.0 / A,
              l.cost_deriv<CC::LeastCost>(psi, psi_soil) +
                  l.LeastCost_beta * phylloptim::umol_to_mol * dVcmax_dpsi /
                      l.leastcost_flux(psi, psi_soil),
              "LeastCost");
      }
    }
  }
  printf("    %d rows | worst relative error %.3e\n", rows, worst);
  ok(rows >= 30, "enough rows reached the identity (hazard 16)");
}

void test_linked_optimum_with_energy_balance_is_stationary() {
  printf("linked curves, energy balance on: the solved leaf optimum is a root of "
         "dJ/dpsi\n");
  using CC = phylloptim::Leaf::CostCurve;
  int rows = 0, pinned = 0;
  double worst = 0.0;
  auto one = [&](auto tag, const Drivers &d, double psi_soil, bool tc) {
    constexpr CC K = decltype(tag)::value;
    phylloptim::Leaf l = make_linked_eb_leaf(d, psi_soil, tc);
    if constexpr (K == CC::LeastCost) {
      // At 146 every row here pins at psi_crit, so the root is never asked
      // for; at 2 every row is interior.
      l.LeastCost_beta = 2.0;
    }
    l.template optimise_psi_leaf_single<K>();
    if constexpr (K == CC::ProfitMax) {
      if (!(l.profitmax_A_max() > 0.0)) {
        return;
      }
    }
    const double psi = l.opt_psi_leaf_;
    if (psi <= psi_soil + 1e-6 || psi >= l.psi_crit - 1e-6) {
      ++pinned;  // a constrained optimum: the gradient is honestly non-zero
      return;
    }
    bool f = false;
    const double at = l.template dprofit_dpsi_leaf<K>(psi, psi_soil, &f);
    if (!f) {
      return;
    }
    // Scaled by the derivative's own magnitude across the bracket, so the bound
    // means "small compared with what this function returns here".
    const double width = l.psi_crit - psi_soil;
    const double scale = std::max(
        std::abs(l.template dprofit_dpsi_leaf<K>(psi_soil + 0.1 * width,
                                                 psi_soil, &f)),
        std::abs(l.template dprofit_dpsi_leaf<K>(psi_soil + 0.9 * width,
                                                 psi_soil, &f)));
    const double resid = std::abs(at) / scale;
    worst = std::max(worst, resid);
    ok(resid < 1e-9, phylloptim::Leaf::curve_name(static_cast<int>(K)) +
                         (tc ? "+TC" : "") + " at Tair " +
                         std::to_string(static_cast<int>(d.leaf_temp)) +
                         " psi_soil " + std::to_string(psi_soil) +
                         ": |dJ/dpsi| at the optimum is at the solver floor");
    ++rows;
  };
  for (double tair : {25.0, 30.0}) {
    Drivers d;
    d.PPFD = 1500.0;
    d.leaf_temp = tair;
    for (double psi_soil : {0.5, 1.0, 2.0}) {
      for (bool tc : {false, true}) {
        one(std::integral_constant<CC, CC::SOX>{}, d, psi_soil, tc);
        one(std::integral_constant<CC, CC::JW26>{}, d, psi_soil, tc);
        one(std::integral_constant<CC, CC::ProfitMax>{}, d, psi_soil, tc);
        one(std::integral_constant<CC, CC::LeastCost>{}, d, psi_soil, tc);
      }
    }
  }
  printf("    %d interior rows (%d pinned) | worst |dJ/dpsi| relative to its "
         "scale %.3e\n",
         rows, pinned, worst);
  ok(rows >= 12, "enough interior optima were checked (hazard 16)");
}

// The longwave half of the energy balance (#97 / #28): an isothermal net radiation
// from a temperature- and humidity-dependent atmospheric emissivity, plus a
// radiative conductance for the leaf's own departure from air temperature.
//
// ⚠️ THE CONTRACT THAT MATTERS IS THAT THE BALANCE STAYS EXPLICIT. A leaf's own
// emission goes as Tleaf^4, so a faithful longwave would make Rn a function of the
// temperature the balance solves for -- an inner iteration, and no analytic
// dTleaf/dE for the collar solve's implicit-function theorem to stand on.
// Linearising the leaf's emission about Tair moves that dependence into the
// DENOMINATOR, so Tleaf stays affine in E. The first block below is what asserts
// that, and it is the assertion to keep if any of the rest is ever rewritten.
void test_longwave_and_radiative_conductance() {
  printf("longwave: isothermal Rn and the radiative conductance\n");
  Drivers d;
  d.PPFD = 900.0;
  d.atm_vpd = 2.0;

  // 1. STILL AFFINE IN E, with the slope the model reports. Three transpirations:
  //    the second difference of an affine function is zero.
  {
    d.leaf_temp = 30.0;
    phylloptim::Leaf l = make_pm_leaf(d, {1.0}, {1.0}, true);
    const double E0 = 0.0, E1 = 1e-5, E2 = 2e-5;
    const double t0 = l.leaf_temp_from_E(E0), t1 = l.leaf_temp_from_E(E1),
                 t2 = l.leaf_temp_from_E(E2);
    near((t2 - t1) - (t1 - t0), 0.0, 1e-9,
         "leaf temperature is still AFFINE in transpiration");
    double dT_dE = 0.0;
    l.leaf_temp_from_E(E1, &dT_dE);
    near(dT_dE, (t2 - t0) / (E2 - E0), 1e-6,
         "and the reported slope is that line's");

    // The slope's closed form, which is the whole point: one conductance, made of
    // the sensible and the radiative parts.
    const double g_total = phylloptim::vol_heat_cap_air / l.ra_ + l.g_rad_;
    near(dT_dE, -phylloptim::latent_heat_vap / g_total, 1e-9,
         "dTleaf/dE is -lambda / (g_sensible + g_radiative)");
  }

  // 2. The radiative conductance is the linearisation it claims to be, and it is
  //    NOT negligible -- an assertion, because a term worth 13% is exactly the
  //    size that gets dropped in a rewrite and not noticed.
  {
    d.leaf_temp = 25.0;
    phylloptim::Leaf l = make_pm_leaf(d, {1.0}, {1.0}, true);
    const double Tair_K = l.Tair_ + phylloptim::C_to_K;
    ok(l.g_rad_ == 4.0 * phylloptim::leaf_emissivity *
                       phylloptim::stefan_boltzmann * Tair_K * Tair_K * Tair_K,
       "g_rad is 4*eps*sigma*Tair^3, bit-for-bit");
    const double g_sens = phylloptim::vol_heat_cap_air / l.ra_;
    ok(l.g_rad_ / (l.g_rad_ + g_sens) > 0.05,
       "and it is more than 5% of the total conductance -- not negligible");

    // Radiative coupling pulls the leaf TOWARD air temperature: the same Rn over a
    // larger conductance is a smaller departure.
    const double with_rad = l.leaf_temp_from_E(0.0);
    const double without_rad = l.Tair_ + l.Rn_ / g_sens;
    ok(std::abs(with_rad - l.Tair_) < std::abs(without_rad - l.Tair_),
       "a radiatively coupled leaf sits closer to air temperature");
  }

  // 3. The atmospheric emissivity, against the form it is usually quoted in.
  //    ⚠️ THIS IS THE UNIT TRAP: `0.642*(ea/T)^(1/7)` needs ea in PASCALS, and the
  //    same relation is more often written `1.24*(ea/T)^(1/7)` with ea in hPa.
  //    Using 0.642 with kPa gives ~0.31 where the answer is ~0.82 -- finite,
  //    plausible, and it triples the longwave loss.
  {
    d.leaf_temp = 25.0;
    d.atm_vpd = 1.0;
    phylloptim::Leaf l = make_pm_leaf(d, {1.0}, {1.0}, true);
    const double Tair_K = l.Tair_ + phylloptim::C_to_K;
    const double ea_Pa =
        1000.0 * (l.saturation_vapour_pressure(l.Tair_) - l.atm_vpd_);
    const double ema_ours =
        phylloptim::atmos_emissivity_coef *
        std::pow(ea_Pa / Tair_K, phylloptim::atmos_emissivity_exponent);
    // The hPa form, independently written out.
    const double ema_hPa = 1.2395 * std::pow((ea_Pa / 100.0) / Tair_K,
                                             1.0 / 7.0);
    near(ema_ours, ema_hPa, 1e-4, "the Pa and hPa forms of Brutsaert agree");
    ok(ema_ours > 0.7 && ema_ours < 0.95,
       "and a clear-sky emissivity lands in 0.7-0.95, not 0.3");

    // Rn_ is the ISOTHERMAL net radiation: absorbed shortwave minus that term.
    const double sw = l.alpha_PAR_ * l.PPFD_ / phylloptim::umol_par_per_joule +
                      l.alpha_NIR_ * l.NIR_int_;
    near(l.Rn_,
         sw - (1.0 - ema_ours) * phylloptim::stefan_boltzmann * Tair_K * Tair_K *
                  Tair_K * Tair_K,
         1e-9, "Rn_ is the isothermal net radiation");
    ok(l.Rn_ - sw < -40.0,
       "and at 25 C the net longwave is stronger than a fixed -40 offset");
  }

  // 4. THE BOUND, which is physical and not padding. Brutsaert's fit runs past
  //    emissivity 1 at high air temperature, which would make the sky heat the
  //    leaf. A soft minimum holds it below 1, so the isothermal term goes to zero
  //    from below -- smoothly, where a min() at 1 put a corner in Rn's slope.
  {
    d.leaf_temp = 60.0;
    d.atm_vpd = 2.0;
    phylloptim::Leaf hot = make_pm_leaf(d, {1.0}, {1.0}, true);
    const double Tair_K = hot.Tair_ + phylloptim::C_to_K;
    const double ea_Pa =
        1000.0 * (hot.saturation_vapour_pressure(hot.Tair_) - hot.atm_vpd_);
    const double unclamped =
        phylloptim::atmos_emissivity_coef *
        std::pow(ea_Pa / Tair_K, phylloptim::atmos_emissivity_exponent);
    ok(unclamped > 1.0,
       "the premise: the fit really does exceed 1 here, so the bound is live");
    const double sw = hot.alpha_PAR_ * hot.PPFD_ / phylloptim::umol_par_per_joule +
                      hot.alpha_NIR_ * hot.NIR_int_;
    ok(hot.Rn_ - sw <= 0.0 && hot.Rn_ - sw > -1e-6,
       "past the fit's unity the isothermal longwave is zero to 1e-6, never positive");

    // Across the fit's unity (43.9 C at 2 kPa) the longwave's slope in Tair is
    // continuous: refined to 0.01 K, the largest jump between neighbouring
    // slopes is a small fraction of the largest slope. The min() gave 1.0.
    std::vector<double> lw;
    const double h = 0.01;
    for (double t = 40.0; t <= 48.0 + 1e-9; t += h) {
      Drivers dt = d;
      dt.leaf_temp = t;
      phylloptim::Leaf x = make_pm_leaf(dt, {1.0}, {1.0}, true);
      lw.push_back(x.Rn_ - (x.alpha_PAR_ * x.PPFD_ / phylloptim::umol_par_per_joule +
                            x.alpha_NIR_ * x.NIR_int_));
    }
    double mx = 0.0, jump = 0.0;
    for (size_t i = 1; i < lw.size(); ++i) mx = std::max(mx, std::abs(lw[i] - lw[i - 1]) / h);
    for (size_t i = 2; i < lw.size(); ++i)
      jump = std::max(jump, std::abs((lw[i] - lw[i - 1]) - (lw[i - 1] - lw[i - 2])) / h / mx);
    printf("    longwave across the emissivity bound: largest slope jump %.3e of scale at 0.01 K\n", jump);
    ok(jump < 0.05, "the longwave's slope in Tair has no corner at the emissivity bound");
    for (double v : lw) {
      if (v > 0.0) { ok(false, "the sky never heats the leaf across the bound"); break; }
    }
  }

  // 5. A deficit at or past saturation is not a physical atmosphere, and it is
  //    refused on both paths rather than floored.
  for (bool eb : {true, false}) {
    d.leaf_temp = 5.0;      // esat ~ 0.87 kPa
    d.atm_vpd = 5.0;        // deficit far beyond it
    const std::string msg = message_of([&] { make_pm_leaf(d, {1.0}, {1.0}, eb); });
    ok(msg.find("atm_vpd") != std::string::npos && msg.find("e_s") != std::string::npos,
       std::string("a deficit past saturation is refused, naming atm_vpd and e_s, energy balance ") +
           (eb ? "on: " : "off: ") + msg);
  }

  // 6. Gate off, nothing reads any of it: Tleaf is the driver. This is the arm the
  //    golden file covers.
  {
    d.leaf_temp = 25.0;
    d.atm_vpd = 2.0;
    phylloptim::Leaf off = make_pm_leaf(d, {1.0}, {1.0}, false);
    off.find_root_collar_psi();
    ok(off.Tleaf_ == off.leaf_temp_, "gate off: Tleaf is still the driver");
  }

  // 7. THE ESCAPE HATCH FOR A COMBINED RESISTANCE IS A CONTRACT. A caller that
  //    injects a convective-plus-radiative resistance into ra_ (Jones et al.'s
  //    r_HR) sets g_rad_ = 0, and must then get the single-resistance balance
  //    back, value AND slope. To rounding, not bit-for-bit: dividing by
  //    vol_heat_cap_air/ra_ and multiplying by ra_/vol_heat_cap_air round apart.
  {
    d.leaf_temp = 30.0;
    d.atm_vpd = 2.0;
    phylloptim::Leaf l = make_pm_leaf(d, {1.0}, {1.0}, true);
    l.Rn_ = 400.0;
    l.ra_ = 12.0;
    l.g_rad_ = 0.0;
    for (double E : {0.0, 5e-5, 1e-4}) {
      double dT_dE = 0.0;
      const double got = l.leaf_temp_from_E(E, &dT_dE);
      const double want = l.Tair_ + (l.Rn_ - phylloptim::latent_heat_vap * E) *
                                        l.ra_ / phylloptim::vol_heat_cap_air;
      near(got, want, 1e-14, "g_rad_ = 0: Tleaf is Tair + (Rn - lambda*E)*ra/rho_cp");
      near(dT_dE, -phylloptim::latent_heat_vap * l.ra_ / phylloptim::vol_heat_cap_air,
           1e-14, "g_rad_ = 0: and dTleaf/dE is -lambda*ra/rho_cp");
    }
  }
}

}  // namespace

// A short, well-plumbed leaf whose stem could carry enough water to cool it past
// the Tetens pole: no operating point gets there, since the leaf condenses below
// the dew point first. The solve's flux ends at the dew-point flux on both
// routes, so set_physiology accepts the leaf and every solve stays at or above
// the dew point.
void test_energy_balance_flux_ends_at_the_dew_point() {
  printf("energy balance: the solve's flux ends at the dew point\n");
  Drivers d;
  d.h = 0.05;
  d.theta = 3.0 * d.theta;
  d.PPFD = 50.0;
  using K = phylloptim::Leaf::CostCurve;
  for (bool collar : {true, false}) {
    phylloptim::Leaf l = collar ? make_pm_leaf(d, {0.5}, {1.0}, true)
                                : make_linked_eb_leaf(d, 0.5, false);
    const char* route = collar ? "collar" : "leaf";
    const double E_most =
        l.leaf_specific_conductance_max_ * l.stem_curve_limit();
    ok(l.leaf_temp_from_E(E_most) < -237.3,
       "the stem's limit alone would cool the leaf past the Tetens pole");
    const double E_dew = l.dew_point_flux();
    ok(E_dew > 0.0 && E_dew < E_most, "the dew-point flux is the nearer bound");
    near(l.leaf_vpd_at(l.leaf_temp_from_E(E_dew)), phylloptim::vpd_leaf_min,
         1e-9, "at the dew-point flux the deficit is at its floor");
    if (collar) {
      l.find_root_collar_psi();
    } else {
      l.template optimise_psi_leaf_single<K::TF24_floor>();
    }
    ok(std::isfinite(l.profit_), route);
    ok(l.transpiration_ <= E_dew * (1.0 + 1e-9),
       "the operating point's flux is at most the dew-point flux");
  }
  phylloptim::Leaf off = make_pm_leaf(d, {0.5}, {1.0}, false);
  ok(std::isinf(off.flux_ceiling()),
     "with the energy balance off there is no flux ceiling");
}

void run_energy_balance() {
  test_energy_balance_flux_ends_at_the_dew_point();
  test_dprofit_dpsi_leaf_matches_a_finite_difference_energy_balance();
  test_leaf_temperature_is_reported();
  test_shutdown_reports_one_temperature();
  test_zero_E_branch_derives_its_own_block();
  test_energy_balance_path_runs();
  test_longwave_and_radiative_conductance();
  test_pm_wind_speed_validation();
  test_pm_absorbed_shortwave();
  test_pm_leaf_temperature_response();
  test_energy_balance_collar_solve_finds_the_maximum();
  test_energy_balance_gate_off_is_inert();
  test_energy_balance_stomatal_decoupling();
  test_boundary_layer_split_of_the_co2_path();
  test_stomatal_conductance_excludes_the_boundary_layer();
  test_boundary_layer_bounds_the_total_conductance();
  test_leaf_to_air_vpd();
  test_linked_dprofit_with_energy_balance_matches_a_finite_difference();
  test_linked_dprofit_with_energy_balance_is_the_identity_link_scaled();
  test_linked_optimum_with_energy_balance_is_stationary();
}
