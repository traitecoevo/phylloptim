// The thermal-damage layer, on the leaf and in its dynamics.

#include "leaf_harness.hpp"

namespace {

// Every input NaN: read anywhere, it reaches the outputs.
void poison_atls_inputs(phylloptim::Leaf &l) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  l.topt_offset_ = nan;
  l.k_i_ = nan;
  l.k_i_Ea_ = nan;
  l.k_i_Tiso_ = nan;
  l.k_i_baseline_ = nan;
  l.k_rec_ = nan;
  l.k_rec_Ea_ = nan;
  l.dTopt_max_ = nan;
  l.dTimpair_max_ = nan;
  l.K_A_ = nan;
  l.A_ = nan;
  l.I_ = nan;
  l.D_ = nan;
  l.k_mat_ = nan;
  l.k_mat_Ea_ = nan;
  l.acclim_alpha_ = nan;
  l.acclim_beta_ = nan;
  l.acclim_T_ = nan;
  l.acclim_s_ = nan;
  l.c_acclim_maint_ = nan;
  l.c_repair_maint_ = nan;
  l.c_protect_maint_ = nan;
  l.c_repair_flux_ = nan;
}

// make_pm_leaf with the layer's gate and inputs set BEFORE the drivers, which is
// the order a caller has to use: set_physiology is what seats the block.
phylloptim::Leaf make_thermal_leaf(const Drivers &d, std::vector<double> psi_soil,
                                   std::vector<double> soil_depth, bool eb,
                                   bool gate,
                                   void (*inputs)(phylloptim::Leaf &) = nullptr) {
  phylloptim::Leaf l = fixture::multilayer_leaf();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  l.use_energy_balance_ = eb;
  l.use_thermal_damage_ = gate;
  l.NIR_int_ = d.PPFD / phylloptim::umol_par_per_joule;  // daylight
  if (inputs != nullptr) {
    inputs(l);
  }
  std::vector<double> root(psi_soil.size(),
                           1.0 / double(psi_soil.size()) / d.area_leaf);
  l.set_physiology(fixture::root_network(root, soil_depth), d.PPFD, psi_soil,
                   soil_depth, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                   d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  return l;
}

// The drivers make_thermal_leaf seated, again: what a caller does each step after
// moving the layer's state.
void reseat_thermal_leaf(phylloptim::Leaf &l, const Drivers &d,
                         const std::vector<double> &psi_soil,
                         const std::vector<double> &soil_depth) {
  std::vector<double> root(psi_soil.size(),
                           1.0 / double(psi_soil.size()) / d.area_leaf);
  l.set_physiology(fixture::root_network(root, soil_depth), d.PPFD, psi_soil,
                   soil_depth, d.K_s * d.theta / d.h, d.atm_vpd, d.ca,
                   d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
}

// The block update_temperature_dependent_params writes, as one comparable value.
std::array<double, 12> thermal_block(const phylloptim::Leaf &l) {
  return {l.vcmax_, l.jmax_, l.gamma_, l.ko_, l.kc_, l.R_d_, l.km_,
          l.electron_transport_, l.k_imp_eff_, l.k_rec_eff_, l.functional_,
          l.photo_temp_};
}

// NaN in every input, gate off, against a leaf that never had them: bit-identical
// on both temperature paths and both soil profiles. Paired with the same poison
// reaching vcmax_ once the gate is on, so the identity is evidence about the
// layer rather than about a gate that does nothing.
void test_thermal_damage_gate_off_is_inert() {
  printf("thermal damage, gate off: no input is read (NaN poisoning)\n");
  int compared = 0;
  for (bool eb : {false, true}) {
    for (double temp : {25.0, 40.0}) {
      for (int layers : {1, 3}) {
        Drivers d;
        d.leaf_temp = temp;
        std::vector<double> psi, depth;
        for (int i = 0; i < layers; ++i) {
          psi.push_back(1.5 + 0.4 * i);
          depth.push_back(0.5 * (i + 1));
        }
        phylloptim::Leaf clean = make_thermal_leaf(d, psi, depth, eb, false);
        phylloptim::Leaf poisoned =
            make_thermal_leaf(d, psi, depth, eb, false, poison_atls_inputs);
        clean.find_root_collar_psi();
        poisoned.find_root_collar_psi();
        const std::string at = " at Tair=" + std::to_string(int(temp)) +
                               " eb=" + std::to_string(eb) +
                               " layers=" + std::to_string(layers);
        ok(clean.profit_ == poisoned.profit_, "profit is bit-identical" + at);
        ok(clean.assim_colimited_ == poisoned.assim_colimited_,
           "assimilation is bit-identical" + at);
        ok(clean.opt_root_psi_ == poisoned.opt_root_psi_,
           "the collar is bit-identical" + at);
        ok(clean.Tleaf_ == poisoned.Tleaf_, "Tleaf is bit-identical" + at);
        ok(thermal_block(clean) == thermal_block(poisoned),
           "the temperature block is bit-identical" + at);
        ok(poisoned.k_imp_eff_ == 0.0 && poisoned.functional_ == 1.0 &&
               poisoned.k_rec_eff_ == 0.0,
           "the outputs are zero" + at);
        const double p0 = clean.opt_root_psi_;
        ok(clean.dprofit_droot_collar_psi(p0) ==
               poisoned.dprofit_droot_collar_psi(p0),
           "dprofit is bit-identical" + at);
        if (layers == 1) {
          clean.set_model(phylloptim::Leaf::CostCurve::TF24, false);
          poisoned.set_model(phylloptim::Leaf::CostCurve::TF24, false);
          clean.optimise();
          poisoned.optimise();
          ok(clean.opt_psi_leaf_ == poisoned.opt_psi_leaf_ &&
                 clean.profit_ == poisoned.profit_,
             "the leaf route is bit-identical" + at);
        }
        ++compared;

        phylloptim::Leaf on =
            make_thermal_leaf(d, psi, depth, eb, true, poison_atls_inputs);
        ok(std::isnan(on.vcmax_) && std::isnan(on.k_imp_eff_) &&
               std::isnan(on.functional_),
           "the same poison reaches the block with the gate on" + at);
      }
    }
  }
  ok(compared == 8, "every configuration was compared");
}

// ⚠️ THE CACHE INVARIANT. The temperature cache's key does not hold A_, I_ or
// D_, so set_physiology must not serve it while the gate is on: re-seating at
// unchanged drivers after the caller moves the state has to recompute. And the
// gate itself is in the key, so switching it off at unchanged drivers must not be
// served the gate-on block.
void test_thermal_damage_bypasses_the_cache() {
  printf("thermal damage: with the gate on, a re-seat reads the current state\n");
  Drivers d;
  // Not 25 C: the peaked Arrhenius is normalised to its reference value there,
  // so no entropy shift moves vcmax at 25 C and the acclimation check would be
  // blind.
  d.leaf_temp = 35.0;
  const std::vector<double> psi{2.0}, depth{1.0};
  phylloptim::Leaf l = make_thermal_leaf(d, psi, depth, false, true);
  const double vcmax_undamaged = l.vcmax_;
  l.I_ = 0.25;
  reseat_thermal_leaf(l, d, psi, depth);
  near(l.vcmax_ / vcmax_undamaged, 0.75, 1e-14,
       "a new damage state at unchanged drivers reaches vcmax");
  const double vcmax_damaged = l.vcmax_;
  l.A_ = 2.0;
  reseat_thermal_leaf(l, d, psi, depth);
  ok(l.vcmax_ != vcmax_damaged,
     "a new acclimation state at unchanged drivers reaches vcmax");

  l.use_thermal_damage_ = false;
  reseat_thermal_leaf(l, d, psi, depth);
  l.find_root_collar_psi();
  phylloptim::Leaf never = make_thermal_leaf(d, psi, depth, false, false);
  never.find_root_collar_psi();
  ok(thermal_block(l) == thermal_block(never) && l.profit_ == never.profit_,
     "switched off at unchanged drivers, the leaf is one never switched on");
}

// The layer has no corners left: the exponentials are the functions themselves,
// A, I and protection are read as themselves, and F's floor is soft. Each is
// asserted where the old clamp would have bound, and inside, where nothing may
// have moved.
void test_thermal_damage_has_no_corners() {
  printf("thermal damage: no clamp, no max, a soft floor on F\n");
  using T = odelia::ode::tangent_scalar<double>;
  // 1. The softplus is its plain form inside [-30, 30] and the function past it,
  //    with a non-zero slope on both sides of each end.
  {
    bool inside_same = true;
    for (double z = -30.0; z <= 30.0; z += 0.5) {
      inside_same &= phylloptim::Leaf::softplus_(z / 2.0, 2.0) == std::log1p(std::exp(z)) / 2.0;
    }
    ok(inside_same, "softplus is bit-identical to log1p(exp(z))/s inside [-30, 30]");
    near(phylloptim::Leaf::softplus_(-20.0, 2.0), std::log1p(std::exp(-40.0)) / 2.0, 1e-30,
         "below -30 it is the function, not its value at -30");
    near(phylloptim::Leaf::softplus_(20.0, 2.0), 20.0, 1e-15, "above 30 it is x to rounding");
    T x_lo(-20.0), x_hi(20.0);
    odelia::ode::set_direction(x_lo, 1.0);
    odelia::ode::set_direction(x_hi, 1.0);
    const double lo = odelia::ode::derivative_along(phylloptim::damage::softplus(x_lo, T(2.0)));
    const double hi = odelia::ode::derivative_along(phylloptim::damage::softplus(x_hi, T(2.0)));
    ok(lo > 0.0 && std::abs(hi - 1.0) < 1e-15,
       "and its slope is the logistic's on both sides, not the clamp's zero");
  }
  // 2. A high-Ea species below 25 C: the impairment rate keeps falling, and its
  //    sensitivity to Ea is the exponential's, where the clamp held both.
  {
    phylloptim::Leaf l = fixture::multilayer_leaf();
    l.k_i_Ea_ = 619e3;
    const double r10 = l.damage_rates(10.0, 0.0, 0.0, 0.0)[1];
    const double r5 = l.damage_rates(5.0, 0.0, 0.0, 0.0)[1];
    ok(r10 > 0.0 && r5 > 0.0 && r5 < r10,
       "a 619 kJ species' impairment keeps falling below 10 C (argument < -30)");
    const std::vector<double> g =
        l.damage_trajectory_gradient({0.0, 1.0}, {5.0, 5.0}, 0.0, 0.0, 0.0, {1});
    ok(g[1 + 2] != 0.0, "and I's sensitivity to k_i_Ea_ there is not the clamp's zero");
    // A requested subset, reordered and repeated, is the full gradient's
    // columns: one pass carries every direction and the request picks them.
    std::vector<int> all(17);
    for (int k = 0; k < 17; ++k) all[std::size_t(k)] = k;
    const std::vector<double> full = l.damage_trajectory_gradient(
        {0.0, 0.5, 1.0}, {5.0, 30.0, 45.0}, 0.2, 0.1, 0.05, all);
    const std::vector<int> some{16, 5, 14, 1, 16};
    const std::vector<double> part = l.damage_trajectory_gradient(
        {0.0, 0.5, 1.0}, {5.0, 30.0, 45.0}, 0.2, 0.1, 0.05, some);
    bool same = part.size() == 3 * 3 * some.size();
    for (std::size_t k = 0; same && k < some.size(); ++k) {
      for (std::size_t j = 0; j < 9; ++j) {
        same = same && part[j + 9 * k] == full[j + 9 * std::size_t(some[k])];
      }
    }
    ok(same, "a reordered, repeated subset of inputs is those columns, bit for bit");
  }
  // 3. F's soft floor: the plain fraction wherever the correction is below an
  //    ulp, positive and falling through I + D = 1 (it underflows to zero only
  //    some 700 widths past it).
  {
    Drivers d;
    double prev = 0.0;
    bool rising = true, positive = true, exact = true;
    for (double sum : {0.0, 0.5, 0.99, 0.9999, 1.0 - 1e-5, 1.0, 1.0 + 1e-6, 1.0 + 1e-5}) {
      phylloptim::Leaf l = make_thermal_leaf(d, {1.0}, {1.0}, false, true,
                                             +[](phylloptim::Leaf &) {});
      l.I_ = sum;
      l.D_ = 0.0;
      reseat_thermal_leaf(l, d, {1.0}, {1.0});
      positive &= l.functional_ > 0.0;
      if (1.0 - sum > 1e-4) exact &= l.functional_ == 1.0 - sum;
      if (sum > 0.0) rising &= l.functional_ < prev;
      prev = l.functional_;
    }
    ok(exact, "F is exactly 1 - I - D away from the floor");
    ok(positive && rising, "and positive and strictly falling in I + D through 1");
  }
  // 4. Protection is signed: a leaf less protected than the baseline pays less
  //    respiration than one at it, which pays less than a more protected one.
  {
    Drivers d;
    auto rd = [&](double k_i) {
      phylloptim::Leaf l = make_thermal_leaf(d, {1.0}, {1.0}, false, true,
                                             +[](phylloptim::Leaf &x) { x.c_protect_maint_ = 0.01; });
      l.k_i_ = k_i;
      reseat_thermal_leaf(l, d, {1.0}, {1.0});
      return l.R_d_;
    };
    const double base = phylloptim::Leaf().k_i_baseline_;
    const double a = rd(base - 10.0), b = rd(base), c = rd(base + 10.0);
    near((a - b) - (b - c), 0.0, 1e-12, "protection's cost is linear through the baseline");
    ok(a > b && b > c, "and signed: less protection than the baseline is a saving");
  }
  // 5. Refusals: a negative fraction or the shift's pole, as infeasible where a
  //    state is seated; a shifted isokinetic temperature out of range, as an error.
  {
    Drivers d;
    auto seat = [&](double A, double I, double D) -> std::string {
      phylloptim::Leaf l = make_thermal_leaf(d, {1.0}, {1.0}, false, true,
                                             +[](phylloptim::Leaf &) {});
      l.A_ = A; l.I_ = I; l.D_ = D;
      try {
        reseat_thermal_leaf(l, d, {1.0}, {1.0});
      } catch (const phylloptim::util::infeasible_error &e) {
        return e.what();
      }
      return "";
    };
    ok(seat(0.0, -1e-9, 0.0).find("damage_state") != std::string::npos,
       "a negative I_ is refused as infeasible");
    ok(seat(0.0, 0.0, -1e-9).find("damage_state") != std::string::npos,
       "a negative D_ likewise");
    ok(seat(-phylloptim::Leaf().K_A_, 0.0, 0.0).find("damage_state") != std::string::npos,
       "and A_ at the shift's pole");
    ok(seat(-1e-9, 0.0, 0.0).empty(), "a hair below zero in A_ is read as itself");
    phylloptim::Leaf l = fixture::multilayer_leaf();
    l.k_i_Tiso_ = 90.0;
    ok(message_of([&] { l.damage_rates(25.0, 0.0, 0.0, 0.0); }).find("[250, 350]") !=
           std::string::npos,
       "an isokinetic temperature outside [250, 350] K is refused, not held");
  }
}

// The solve root-finds dprofit == 0, so with the layer on dprofit has to be the
// derivative of the profit the solve REPORTS with the layer on. It is, because
// both read the block update_temperature_dependent_params writes -- and on the
// energy-balance path both re-seat it at each candidate's own leaf temperature.
// Collar route against a difference of evaluate_root_collar_psi, leaf route
// against a difference of profit_psi_leaf_TF24, at every input non-default.
void test_thermal_damage_dprofit_is_the_solves_own() {
  printf("thermal damage, gate on: dprofit is the derivative of the reported "
         "profit\n");
  using CC = phylloptim::Leaf::CostCurve;
  int collar_rows = 0, leaf_rows = 0, stationary = 0;
  double worst_collar = 0.0, worst_leaf = 0.0, worst_resid = 0.0;
  for (bool eb : {false, true}) {
    for (double temp : {30.0, 40.0}) {
      Drivers d;
      d.PPFD = 1500.0;
      d.leaf_temp = temp;
      const struct {
        std::vector<double> psi;
        std::vector<double> depth;
      } soils[] = {{{1.0}, {1.0}}, {{1.0, 2.0, 3.0}, {0.5, 1.0, 1.5}}};
      for (const auto &s : soils) {
        phylloptim::Leaf l =
            make_thermal_leaf(d, s.psi, s.depth, eb, true, set_atls_inputs);
        l.find_root_collar_psi();
        const double p0 = l.opt_root_psi_;
        const std::string at = " at Tair=" + std::to_string(int(temp)) +
                               " eb=" + std::to_string(eb) + " layers=" +
                               std::to_string(s.psi.size());
        if (l.operating_point_kind() ==
            phylloptim::Leaf::OperatingPointKind::Interior) {
          const double r = std::abs(l.dprofit_droot_collar_psi(p0));
          worst_resid = std::max(worst_resid, r);
          ok(r < 1e-8, "the solved collar is stationary" + at);
          ++stationary;
        }
        for (double target : {p0 - 0.05, p0 - 0.1}) {
          const double eps = 1e-5;
          const double analytic = l.dprofit_droot_collar_psi_for<CC::TF24>(target);
          // A target wetter than the wet bound is refused, not a row.
          double up = NAN, dn = NAN;
          try {
            up = l.evaluate_root_collar_psi_for<CC::TF24>(target + eps);
            dn = l.evaluate_root_collar_psi_for<CC::TF24>(target - eps);
          } catch (const phylloptim::util::infeasible_error &) {
            continue;
          }
          if (!std::isfinite(up) || !std::isfinite(dn)) {
            continue;
          }
          const double fd = (up - dn) / (2 * eps);
          worst_collar = std::max(worst_collar, std::abs(analytic - fd));
          near(analytic, fd, 1e-4, "collar dprofit at " +
                                       std::to_string(target) + at);
          ++collar_rows;
        }
        if (s.psi.size() == 1) {
          // As test_dprofit_dpsi_leaf_matches_a_finite_difference: h = 1e-3 is
          // the floor of that difference's V.
          const double h = 1e-3, psi_soil = s.psi[0];
          for (double frac : {0.25, 0.5}) {
            const double psi = psi_soil + frac * (l.psi_crit - psi_soil);
            bool feasible = false;
            const double an =
                l.dprofit_dpsi_leaf<CC::TF24>(psi, psi_soil, &feasible);
            const double fd = (l.profit_psi_leaf_TF24(psi + h, psi_soil) -
                               l.profit_psi_leaf_TF24(psi - h, psi_soil)) /
                              (2.0 * h);
            if (!feasible) {
              continue;
            }
            worst_leaf = std::max(worst_leaf, std::abs(an / fd - 1.0));
            near(an / fd, 1.0, 5e-5, "leaf dprofit at frac=" +
                                         std::to_string(frac) + at);
            ++leaf_rows;
          }
        }
      }
    }
  }
  printf("    %d collar rows (worst |error| %.3e), %d leaf rows (worst relative "
         "%.3e), %d stationary (worst |dprofit| %.3e)\n",
         collar_rows, worst_collar, leaf_rows, worst_leaf, stationary,
         worst_resid);
  ok(collar_rows >= 12 && leaf_rows >= 6 && stationary >= 4,
     "the sweep compared enough rows of each kind");
}

// temp_deriv and dprofit_energy_balance_term perturb the block by
// 1e-3 degrees and put it back by assignment. The layer adds three outputs to
// the block, so each must put back eleven numbers, not eight.
void test_thermal_damage_fd_helpers_restore_the_block() {
  printf("thermal damage: the finite-difference helpers restore the block\n");
  Drivers d;
  d.leaf_temp = 30.0;
  phylloptim::Leaf l =
      make_thermal_leaf(d, {1.0, 2.0}, {0.5, 1.0}, true, true, set_atls_inputs);
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == phylloptim::Leaf::OperatingPointKind::Interior,
     "the fixture is an interior point");
  const std::array<double, 12> before = thermal_block(l);
  ok(before[8] > 0.0 && before[9] > 0.0 && before[10] > 0.0,
     "the fixture has non-zero k_imp_eff, k_rec_eff and functional to lose");

  const double dRd = l.temp_deriv(l.Tleaf_ + 2.0, [&] { return l.R_d_; });
  ok(thermal_block(l) == before, "temp_deriv restores the whole block");
  ok(std::isfinite(dRd), "temp_deriv is finite");

  using CC = phylloptim::Leaf::CostCurve;
  const double term = l.dprofit_energy_balance_term<CC::TF24>(
      l.ci_, l.stom_cond_CO2_, 1.0, 1e-5, 1.0, 0.1, 0.1, 0.1, 1.0, -1e3,
      l.Tleaf_ + 2.0, l.transpiration_ - l.residual_flux_, l.residual_flux_,
      0.0, 0.0, 1.0);
  ok(std::isfinite(term), "dprofit_energy_balance_term is finite");
  ok(thermal_block(l) == before, "dprofit_energy_balance_term restores the whole block");
  // The linked arm evaluates A after the restore, so it must leave the same
  // eleven in place too.
  const double linked = l.dprofit_energy_balance_term<CC::SOX>(
      l.ci_, l.stom_cond_CO2_, 1.0, 1e-5, 1.0, 0.1, 0.1, 0.1, 1.0, -1e3,
      l.Tleaf_ + 2.0, l.transpiration_ - l.residual_flux_, l.residual_flux_,
      0.0, 0.0, 1.0);
  ok(std::isfinite(linked), "the linked dprofit_energy_balance_term is finite");
  ok(thermal_block(l) == before,
     "and the linked arm restores the whole block as well");

  // The repair-flux cost reads k_rec_eff(T), so it reaches dR_d/dT: the
  // derivative moves with that coefficient wherever repair has an activation
  // energy, as the fixture's does.
  phylloptim::Leaf no_flux = l;
  no_flux.c_repair_flux_ = 0.0;
  ok(no_flux.temp_deriv(40.0, [&] { return no_flux.R_d_; }) !=
         l.temp_deriv(40.0, [&] { return l.R_d_; }),
     "the repair-flux cost's temperature response reaches dR_d/dT");
}

// Every derivative path that reads the capacities off the parameter pack, the
// trait gradient, and the closed form refuse a leaf with the layer on; each is
// checked for its message, and each beside a control that runs with the gate
// off. Two paths deliberately do NOT refuse, and are checked to run: the supply
// draw, which reads no capacity, and the marginal's curvature, which differences
// the solve's own dprofit.
void test_thermal_damage_refusals() {
  printf("thermal damage: the derivative paths that cannot carry it refuse\n");
  using CC = phylloptim::Leaf::CostCurve;
  const std::string refused = "refused with use_thermal_damage_ on";
  Drivers d;
  d.PPFD = 1500.0;
  // 35 C, not 25: the gradient check below needs a finite difference clean enough
  // to referee the composite. At 25 C differencing the solve has a floor of ~2e-4
  // relative for inputs that only shift respiration a little (c_repair_flux_ moves
  // profit ~12x less per step than R_d_25, whose own floor there is 1.8e-5 with the
  // layer OFF); at 35 C every column agrees to ~1e-9.
  d.leaf_temp = 35.0;
  const std::vector<double> psi{1.0, 2.0, 3.0}, depth{0.5, 1.0, 1.5};
  phylloptim::Leaf l =
      make_thermal_leaf(d, psi, depth, false, false, set_atls_inputs);
  l.find_root_collar_psi();
  phylloptim::Leaf on =
      make_thermal_leaf(d, psi, depth, false, true, set_atls_inputs);
  on.find_root_collar_psi();
  ok(on.operating_point_kind() ==
         phylloptim::Leaf::OperatingPointKind::Interior,
     "the gate-on fixture is an interior point");
  ok(std::isfinite(profit_row_at_current_point(l, phylloptim::par_vcmax_25)),
     "control: the derivative surface answers with the gate off");

  ok(throws_with([&] { profit_row_at_current_point(on, phylloptim::par_vcmax_25); },
                 "Leaf::operating_collar: " + refused),
     "operating_collar refuses");
  const phylloptim::leaf_pars<double> pars = on.passive_pars();
  const auto on_supply = on.roots_.held_supply();
  const auto draw = on.supply_draw_at<double>(on.opt_root_psi_, on_supply);
  ok(std::isfinite(draw.E_up.value), "supply_draw_at reads no capacity, and runs");
  ok(std::isfinite(on.profit_curvature_at_collar<CC::TF24>()),
     "profit_curvature_at_collar differences the solve's own dprofit, and runs");
  ok(throws_with([&] { on.operating_collar<CC::TF24, double>(draw, pars); },
                 "Leaf::operating_collar: " + refused),
     "operating_collar refuses at double too");
  ok(throws_with([&] { on.billed_outputs_at<CC::TF24, double>(on.opt_root_psi_, draw, pars); },
                 "Leaf::billed_outputs_at: " + refused),
     "billed_outputs_at refuses");
  ok(throws_with([&] { on.vcmax_at<double>(pars); }, "Leaf::vcmax_at: " + refused),
     "vcmax_at refuses");
  ok(throws_with([&] { on.jmax_at<double>(pars); }, "Leaf::jmax_at: " + refused),
     "jmax_at refuses");
  ok(throws_with([&] { on.respiration_at<double>(pars); },
                 "Leaf::respiration_at: " + refused),
     "respiration_at refuses");
  ok(throws_with([&] { on.photo_capacity_at<double>(pars); },
                 "Leaf::vcmax_at: " + refused),
     "photo_capacity_at refuses, through vcmax_at");
  ok(throws_with([&] { on.photo_capacity_at<double>(pars, on.photo_temp_); },
                 "Leaf::vcmax_at: " + refused),
     "and at a scalar leaf temperature too");

  // The trait gradient, both entry points.
  double theta[phylloptim::gradient::n_theta];
  for (int i = 0; i < phylloptim::gradient::n_theta; ++i) {
    theta[i] = pars[std::size_t(i)];
  }
  phylloptim::gradient::Drivers gd;
  gd.root_network = fixture::root_network(
      std::vector<double>(3, 1.0 / 3.0 / d.area_leaf), depth);
  gd.PPFD = d.PPFD;
  gd.psi_soil = psi;
  gd.soil_depth = depth;
  gd.atm_vpd = d.atm_vpd;
  gd.ca = d.ca;
  gd.leaf_temp = d.leaf_temp;
  gd.atm_o2_kpa = d.atm_o2_kpa;
  gd.atm_kpa = d.atm_kpa;
  const std::vector<phylloptim::gradient::Drivers> obs{gd, gd};
  const int par[] = {phylloptim::par_vcmax_25};
  phylloptim::gradient::Settings s;
  s.model = {phylloptim::Leaf::CostCurve::TF24,
             phylloptim::gradient::Route::Collar};
  phylloptim::gradient::Result r;
  // The trait gradient CARRIES the layer: both routes difference the solve's own
  // dprofit and outputs, which read the block the layer writes. So with the gate
  // on it runs, and the IFT composite agrees with differencing the solve.
  // The gate-on leaf's own theta: its damage slots carry the inputs it was given.
  double theta_on[phylloptim::gradient::n_theta];
  {
    const phylloptim::leaf_pars<double> p_on = on.passive_pars();
    for (int i = 0; i < phylloptim::gradient::n_theta; ++i) {
      theta_on[i] = p_on[std::size_t(i)];
    }
  }
  // A trait, a damage state and a cost-carried input, through both routes.
  const int pars_on[] = {phylloptim::par_vcmax_25,
                         phylloptim::par_I,
                         phylloptim::par_c_repair_flux};
  phylloptim::gradient::Result r_on;
  phylloptim::gradient::at(on, theta_on, gd, false, pars_on, 3, s, r_on);
  ok(r_on.status != phylloptim::gradient::Status::Error,
     "gradient::at runs with the gate on");
  {
    phylloptim::gradient::Settings sfd = s;
    sfd.method = phylloptim::gradient::Method::Fd;
    phylloptim::gradient::Result rfd;
    phylloptim::gradient::at(on, theta_on, gd, false, pars_on, 3, sfd, rfd);
    double worst = 0.0;
    for (std::size_t k = 0; k < r_on.grad.size(); ++k) {
      const double scale = std::max(1.0, std::abs(rfd.grad[k]));
      worst = std::max(worst, std::abs(r_on.grad[k] - rfd.grad[k]) / scale);
    }
    printf("    gate on: IFT against differencing the solve, worst %.3e\n", worst);
    ok(worst < 1e-6,
       "with the gate on, the IFT composite agrees with differencing the solve");
    // dA/dI: impaired area does not photosynthesise, so A must FALL with I.
    ok(r_on.grad[1 * phylloptim::gradient::n_outputs + 0] < 0.0,
       "assimilation falls as the impaired fraction rises");
  }
  const std::vector<phylloptim::gradient::Result> gon =
      phylloptim::gradient::batch(on, theta_on, 1, obs, false, par, 1, s);
  ok(gon.size() == 2 && gon[0].status != phylloptim::gradient::Status::Error,
     "gradient::batch runs with the gate on");
  const std::vector<phylloptim::gradient::Result> g =
      phylloptim::gradient::batch(l, theta, 1, obs, false, par, 1, s);
  ok(g.size() == 2 && g[0].status != phylloptim::gradient::Status::Error &&
         g[1].status != phylloptim::gradient::Status::Error,
     "control: the batch runs with the gate off");

  // The closed form, on the one configuration it otherwise accepts.
  phylloptim::Leaf cf = make_thermal_leaf(d, {1.0}, {1.0}, false, false,
                                          set_atls_inputs);
  cf.set_model(CC::TF24, false, true);
  cf.optimise();
  ok(std::isfinite(cf.opt_psi_leaf_), "control: the closed form runs gate off");
  cf.use_thermal_damage_ = true;
  ok(throws_with([&] { cf.optimise(); },
                 "the closed form is refused with use_thermal_damage_ on"),
     "the closed form refuses");
  cf.set_model(CC::TF24, false);
  cf.optimise();
  ok(std::isfinite(cf.opt_psi_leaf_), "and the exact method runs gate on");
}

}  // namespace

void run_thermal_damage() {
  test_thermal_damage_gate_off_is_inert();
  test_thermal_damage_bypasses_the_cache();
  test_thermal_damage_dprofit_is_the_solves_own();
  test_thermal_damage_has_no_corners();
  test_thermal_damage_fd_helpers_restore_the_block();
  test_thermal_damage_refusals();
}
