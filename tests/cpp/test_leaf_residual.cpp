// The residual (cuticular) conductance, its caps, and the smoothness sweep.

#include "leaf_harness.hpp"

namespace residual {

using Kind = phylloptim::Leaf::OperatingPointKind;

using CC = phylloptim::Leaf::CostCurve;

// The manuscript's leaf on the leaf route: one soil potential, kmax 3.14e-5, D 2 kPa,
// PPFD 900, 25 C.
Drivers manuscript_drivers() {
  Drivers d;
  d.PPFD = 900.0;
  d.atm_vpd = 2.0;
  d.leaf_temp = 25.0;
  // kmax = K_s * theta / h = 3.14e-5 at these three
  d.K_s = 1.0; d.theta = 3.14e-5; d.h = 1.0;
  return d;
}

void test_temperature_factor() {
  printf("residual conductance: the temperature factor is one smooth curve\n");
  phylloptim::Leaf l = fixture::multilayer_leaf();
  ok(l.g_min_factor(25.0) == 1.0, "the factor is exactly 1 at 25 C");
  ok(l.g_min == 0.0, "g_min defaults to zero");
  for (double T : {-10.0, 0.0, 25.0, 40.0, 60.0}) {
    ok(l.g_min_at(T) == 0.0, "at g_min = 0 the conductance is exactly zero at " +
                                 std::to_string(T) + " C");
  }
  double prev = 0.0;
  for (double T = 0.0; T <= 60.0; T += 0.5) {
    const double f = l.g_min_factor(T);
    ok(std::isfinite(f) && f > prev, "the factor rises with temperature");
    prev = f;
    // Its slope, against a central difference of itself.
    const double h = 1e-4;
    const double fd = (l.g_min_factor(T + h) - l.g_min_factor(T - h)) / (2 * h);
    near(l.g_min_factor_temp_deriv(T) / fd, 1.0, 1e-6,
         "the slope matches a central difference at " + std::to_string(T));
  }
  // Smooth across the transition: the second difference changes little from
  // one side of T_t to the other, where two slopes meeting at a point would
  // put all of it in one cell.
  const double Tt = l.g_min_T_trans_;
  auto d2 = [&](double T) {
    const double h = 0.25;
    return (std::log(l.g_min_factor(T + h)) - 2 * std::log(l.g_min_factor(T)) +
            std::log(l.g_min_factor(T - h))) / (h * h);
  };
  const double peak = d2(Tt);
  ok(peak > 0.0, "log g_min is convex at the transition");
  ok(std::abs(d2(Tt + 1.0) - d2(Tt - 1.0)) < 0.2 * peak,
     "and its curvature is continuous there, not concentrated in one cell");
  // The documented asymptotes: k_lo below, k_lo + k_hi above.
  const double lo = (std::log(l.g_min_factor(5.0)) - std::log(l.g_min_factor(0.0))) / 5.0;
  const double hi = (std::log(l.g_min_factor(65.0)) - std::log(l.g_min_factor(60.0))) / 5.0;
  near(lo, l.g_min_k_low_, 1e-6, "the low-temperature slope is k_low");
  near(hi, l.g_min_k_low_ + l.g_min_k_high_, 1e-4,
       "the high-temperature slope is k_low + k_high");
  bool threw = false;
  try { l.g_min_T_width_ = 0.0; l.g_min_factor(30.0); } catch (const std::runtime_error &) { threw = true; }
  ok(threw, "a zero width is refused");
}

void test_cap_is_a_smooth_minimum() {
  printf("residual conductance: the cap is a smooth minimum of demand and ceiling\n");
  using L = phylloptim::Leaf;
  const double E = 1e-5, phi = L::residual_cap_fraction;
  ok(L::residual_cap(0.0, E) == 0.0, "no demand, no flux, exactly");
  ok(L::residual_cap(1e-7, 0.0) == 0.0, "no ceiling, no flux");
  near(L::residual_cap(1e-8, E) / 1e-8, 1.0, 1e-12,
       "well inside the ceiling the cap is inert");
  near(L::residual_cap(1e-6, E) / 1e-6, 1.0, 1e-6,
       "at a fifth of the cap it is inert to 1e-6");
  ok(L::residual_cap(1e-3, E) <= phi * E && L::residual_cap(1e-3, E) > 0.999 * phi * E,
     "far past the ceiling it is the ceiling's fraction");
  // Monotone and smooth in the demand, with the slopes it reports. The slopes
  // are checked where they are O(1): far inside the cap dr*/dE_inf is below
  // rounding and far past it dr*/dr is, so a difference measures noise there.
  double prev = 0.0;
  for (double lr = -9.0; lr <= -2.0; lr += 0.05) {
    const double r = std::pow(10.0, lr), h = 1e-5 * r;
    const double v = L::residual_cap(r, E);
    // Monotone to rounding: on the plateau successive values are the cap's
    // fraction to the last bit, in either order.
    ok(v >= prev * (1 - 1e-12) && v <= phi * E * (1 + 1e-12),
       "r* rises with r and stays under the cap");
    prev = v;
    const double q = r / (phi * E);
    if (q < 0.05 || q > 3.0) continue;
    const double fd_r = (L::residual_cap(r + h, E) - L::residual_cap(r - h, E)) / (2 * h);
    near(L::residual_cap_ddemand(r, E) / fd_r, 1.0, 1e-6, "dr*/dr matches a difference");
    if (q < 0.5) continue;  // below this dr*/dE_inf is under the difference's rounding
    const double hE = 1e-5 * E;
    const double fd_E = (L::residual_cap(r, E + hE) - L::residual_cap(r, E - hE)) / (2 * hE);
    near(L::residual_cap_dceiling(r, E) / fd_E, 1.0, 1e-5, "dr*/dE_inf matches a difference");
  }
  // The limit taken past q = 1e4 joins the formula to rounding.
  const double q_lim = phylloptim::residual_cap_q_limit;
  const double r_at = q_lim * phi * E;
  near(L::residual_cap(r_at * (1 - 1e-9), E), L::residual_cap(r_at * (1 + 1e-9), E), 1e-12,
       "the cap is continuous where it switches to its limit");
}

// At g_min = 0 every new term is an exact zero appended: the conductance is the
// old expression bit for bit, the wet bound is the upstream potential, and the
// shut state's fluxes are zero.
void test_reduces_to_nothing_at_zero() {
  printf("residual conductance: at g_min = 0 the term is an exact zero\n");
  Drivers d;
  for (double psi : {0.5, 2.0, 4.0, 6.0}) {
    phylloptim::Leaf l = make(d, {psi}, {1.0}, 0.0);
    l.find_root_collar_psi();
    const std::string at = " at psi_soil " + std::to_string(psi);
    ok(l.residual_flux_ == 0.0, "the residual flux is exactly zero" + at);
    ok(l.g_min_leaf_ == 0.0, "and so is g_min at the leaf temperature" + at);
    ok(l.cuticular_cond_CO2() == 0.0, "and the cuticle's conductance" + at);
    ok(l.stom_cond_CO2_ ==
           l.atm_kpa_ * l.transpiration_ * phylloptim::kg_to_mol_h2o /
               l.vpd_leaf_ / l.H2O_CO2_stom_diff_ratio_,
       "the conductance is the old expression, bit for bit" + at);
    ok(l.leaf_wet_bound(psi) == psi, "the leaf route's wet end is psi_soil" + at);
  }
}

// E = E_s + r* and E_up = E, at every solved kind, with the residual seated.
void test_water_budget_closes() {
  printf("residual conductance: the water budget closes at every kind\n");
  std::map<std::string, int> kinds;
  int rows = 0;
  double worst_up = 0.0;
  for (double g_min : {kGmin, 4 * kGmin}) {
    for (bool eb : {false, true}) {
      for (int layers : {1, 3}) {
        for (double psi0 : {0.5, 2.0, 4.0, 5.5, 6.0, 7.0, 17.0}) {
          for (double ppfd : {0.0, 50.0, 900.0}) {
            Drivers d;
            d.PPFD = ppfd;
            std::vector<double> psi_soil, depth;
            for (int i = 0; i < layers; ++i) {
              psi_soil.push_back(psi0 + 0.3 * i);
              depth.push_back(1.0 * (i + 1));
            }
            phylloptim::Leaf l = make(d, psi_soil, depth, g_min, eb);
            l.find_root_collar_psi();
            const std::string at = " (g_min " + std::to_string(g_min) + ", eb " +
                                   std::to_string(eb) + ", " + std::to_string(layers) +
                                   " layers, psi " + std::to_string(psi0) + ", PPFD " +
                                   std::to_string(ppfd) + ")";
            ++kinds[l.operating_point_kind_name(l.operating_point_kind())];
            ++rows;
            ok(std::isfinite(l.profit_) && std::isfinite(l.ci_),
               "the solve answers with finite profit and ci" + at);
            // Except past the end of the stem curve, which carries nothing.
            ok(l.residual_flux_ > 0.0 ||
                   l.operating_point_kind() == Kind::HydraulicShutdown,
               "the cuticle loses water" + at);
            const double E_s = l.transpiration_ - l.residual_flux_;
            ok(E_s >= -1e-18, "the stomatal flux is non-negative" + at);
            ok(l.stom_cond_CO2_ >= l.cuticular_cond_CO2() * (1 - 1e-12),
               "the conductance is at least the cuticle's" + at);
            // E_up = E: by continuity at an interior point (to the collar
            // tolerance), by the wet bound's own root at a pinned or shut one.
            const double gap = std::abs(l.E_up_ - l.transpiration_) /
                               std::max(l.transpiration_, 1e-12);
            worst_up = std::max(worst_up, gap);
            ok(gap < 1e-7, "the collar carries what the soil delivers" + at);
            if (l.operating_point_kind() == Kind::BoundarySoil ||
                l.operating_point_kind() == Kind::HydraulicShutdown) {
              // The bound is placed to 1e-15 in the collar, so the stomatal
              // flux there is the flux balance's rounding, ~1e-13 of E.
              ok(std::abs(E_s) <= 1e-11 * l.transpiration_,
                 "a leaf at its wet bound moves only the residual flux" + at);
              // In series with the boundary layer where the energy balance
              // gives it one (#122).
              const double cut = l.cuticular_cond_CO2();
              const double cut_series =
                  l.use_energy_balance_
                      ? 1.0 / (1.0 / cut - l.boundary_layer_co2_excess())
                      : cut;
              near(l.stom_cond_CO2_ / cut_series, 1.0, 1e-9,
                   "and conducts CO2 through its cuticle alone" + at);
            }
            // The residual at the reported ci: supply == demand on the total
            // conductance, to the ci root-find's tolerance.
            const double res =
                l.assim_colimited_ * phylloptim::umol_to_mol -
                l.stom_cond_CO2_ * (l.ca_ - l.ci_) / (l.atm_kpa_ * phylloptim::kPa_to_Pa);
            ok(std::abs(res) < 1e-15, "ci balances supply and demand" + at);
            if (ppfd == 0.0) {
              ok(l.ci_ > l.ca_, "in the dark ci sits above ambient" + at);
              near(l.assim_colimited_, -l.R_d_, 1e-9, "and the leaf nets -R_d" + at);
            }
          }
        }
      }
    }
  }
  std::string names;
  for (const auto &[k, n] : kinds) names += k + "(" + std::to_string(n) + ") ";
  printf("    %d rows | kinds: %s| worst |E_up - E|/E %.2e\n", rows, names.c_str(), worst_up);
  ok(kinds.count("interior") > 0 && kinds.count("boundary-soil") > 0 &&
         kinds.count("hydraulic-shutdown") > 0,
     "the grid reached interior, pinned and shut leaves");
  ok(kinds.count("solver-refused") == 0 && kinds.count("non-finite-gradient") == 0,
     "and refused none");
}

// The manuscript's figure-g0 claims, from the solver rather than from splines
// (#182). At rho = 0 the residual flux is water-only, which is the manuscript's
// model: TF24 gives back 0.59-0.85 of it over psi_soil 0.05-4, CF77 gives back
// nothing, A falls by under 1.4% and psi shifts by under 0.013 MPa. With the
// cuticle's CO2 pathway (the default rho) a constant price gives back exactly
// rho * 1.67 of the flux: the cuticle's conductance to CO2 substitutes for the
// stomata's, one for one on the water the stomata would have spent on it.
void test_manuscript_figure_g0() {
  printf("residual conductance: the manuscript's figure-g0 claims, from the solver\n");
  const Drivers d = manuscript_drivers();
  const double lambda = 1.0e5;  // the manuscript's CF77 price
  for (double rho : {0.0, phylloptim::Leaf().g_min_co2_ratio_}) {
    for (const char *m : {"TF24", "CF77"}) {
      double cmin = 1e9, cmax = -1e9, dA_max = 0.0, dpsi_max = 0.0;
      int n = 0;
      for (int i = 0; i < 22; ++i) {
        const double ps = 0.05 + (4.0 - 0.05) * i / 21.0;
        phylloptim::Leaf a = make(d, {ps}, {1.0}, 0.0, false, true, rho);
        phylloptim::Leaf b = make(d, {ps}, {1.0}, kGmin, false, true, rho);
        a.CF77_lambda_ = lambda;
        b.CF77_lambda_ = lambda;
        a.set_model(m, "leaf");
        b.set_model(m, "leaf");
        a.optimise();
        b.optimise();
        if (a.opt_psi_leaf_ > 0.999 * a.psi_crit) continue;  // at the dry bound the comparison is vacuous
        const double Es_a = a.transpiration_ - a.residual_flux_;
        const double Es_b = b.transpiration_ - b.residual_flux_;
        const double comp = (Es_a - Es_b) / b.residual_flux_;
        cmin = std::min(cmin, comp);
        cmax = std::max(cmax, comp);
        dA_max = std::max(dA_max, 100.0 * (a.assim_colimited_ - b.assim_colimited_) / a.assim_colimited_);
        dpsi_max = std::max(dpsi_max, b.opt_psi_leaf_ - a.opt_psi_leaf_);
        ++n;
      }
      printf("    rho %.2f %-4s: gives back %.3f..%.3f of r | A falls <= %.3f%% | psi shifts <= %.4f MPa (%d soils)\n",
             rho, m, cmin, cmax, dA_max, dpsi_max, n);
      // CF77 at the manuscript's price pins to the dry bound over most of the
      // range, as its own figure does; the interior soils are the comparison.
      ok(n >= (std::string(m) == "TF24" ? 20 : 8),
         std::string(m) + ": the sweep is interior over enough of the soil range");
      if (std::string(m) == "TF24") {
        if (rho == 0.0) {
          ok(cmin > 0.55 && cmin < 0.65, "TF24 gives back at least ~0.59 of the residual flux (water-only)");
          ok(cmax > 0.80 && cmax < 0.90, "and at most ~0.85 of it");
          ok(dA_max < 1.45, "A falls by at most ~1.4%");
          ok(dpsi_max < 0.0135, "psi shifts by at most 0.013 MPa");
        }
        ok(cmin > 0.0 && cmax < 1.0, "TF24 gives back part of the flux, never all: eq-g0shift's range");
      } else {
        // A constant price moves the stomatal flux by exactly the conductance
        // the cuticle now supplies to CO2: rho * 1.67 of the residual flux.
        const double want = rho * phylloptim::Leaf().H2O_CO2_stom_diff_ratio_;
        near(cmin, want, 1e-4, "CF77 gives back exactly rho * 1.67 of the residual flux (min)");
        near(cmax, want, 1e-4, "CF77 gives back exactly rho * 1.67 of the residual flux (max)");
        if (rho == 0.0) {
          ok(cmax < 1e-3, "a water-only residual flux moves a constant-price leaf not at all");
        }
      }
    }
  }
}

// eq-g0shift: dE_s/dr = -Theta'' / (|B'| + Theta''), with B = h'(A) dA/dE and
// Theta the cost as a function of the total flux, against a difference of the
// solve in g_min. Water-only (rho = 0), which is the model the identity is
// derived for; the CO2 pathway adds a term the manuscript does not have.
void test_shift_identity() {
  printf("residual conductance: eq-g0shift holds against a difference of the solve\n");
  const Drivers d = manuscript_drivers();
  for (double ps : {0.5, 1.5, 3.0}) {
    // The solve at the midpoint and one step either side in g_min.
    const double g0 = kGmin, dg = 0.5 * kGmin;
    auto solve = [&](double g) {
      phylloptim::Leaf l = make(d, {ps}, {1.0}, g, false, true, 0.0);
      l.set_model("TF24", "leaf");
      l.optimise();
      return l;
    };
    phylloptim::Leaf lo = solve(g0 - dg), mid = solve(g0), hi = solve(g0 + dg);
    const double lhs = ((hi.transpiration_ - hi.residual_flux_) -
                        (lo.transpiration_ - lo.residual_flux_)) /
                       (hi.residual_flux_ - lo.residual_flux_);
    // Theta'' and B' at the midpoint's optimum, from the solver's own leaf
    // route: lambda_TF24 is dC/dE exactly, and A(E) is read off the evaluator
    // at the two neighbouring leaf potentials.
    const double psi = mid.opt_psi_leaf_, h = 2e-3;
    auto E_of = [&](double p) { return mid.transpiration(p, ps); };
    auto A_of = [&](double p) { mid.profit_psi_leaf_for<CC::TF24>(p, ps); return mid.assim_colimited_; };
    const double Em = E_of(psi - h), Ep = E_of(psi + h);
    const double theta2 = (mid.lambda_TF24(psi + h) - mid.lambda_TF24(psi - h)) / (Ep - Em);
    // B(E) = dA/dE, itself by a central difference, at psi +- h/2.
    auto B_at = [&](double p) {
      const double q = h / 2;
      return (A_of(p + q) - A_of(p - q)) / (E_of(p + q) - E_of(p - q));
    };
    const double Bp = (B_at(psi + h) - B_at(psi - h)) / (Ep - Em);
    const double rhs = -theta2 / (std::abs(Bp) + theta2);
    printf("    psi_soil %.1f: dE_s/dr = %.4f (solve)  %.4f (identity)  Theta'' %.3e  B' %.3e\n",
           ps, lhs, rhs, theta2, Bp);
    ok(theta2 > 0.0 && Bp < 0.0, "the cost is convex and the benefit concave at psi_soil " + std::to_string(ps));
    ok(lhs > -1.0 && lhs < 0.0, "the leaf gives back part of the flux, never all");
    near(lhs / rhs, 1.0, 2e-2, "eq-g0shift holds to the second differences' precision at psi_soil " + std::to_string(ps));
  }
}

// The g_min rows on the tape against a difference of the solve, at interior,
// pinned and shut points, with the water rows at the shut one.
void test_rows_match_a_difference() {
  printf("residual conductance: the g_min rows match a difference of the solve\n");
  int compared = 0;
  std::map<std::string, int> seen;
  for (int layers : {1, 3}) {
    for (double psi0 : {1.0, 3.0, 4.5, 6.0, 17.0}) {
      for (double ppfd : {10.0, 900.0}) {
        Drivers d;
        d.PPFD = ppfd;
        std::vector<double> psi_soil, depth;
        for (int i = 0; i < layers; ++i) { psi_soil.push_back(psi0 + 0.3 * i); depth.push_back(i + 1.0); }
        phylloptim::Leaf l = make(d, psi_soil, depth, kGmin);
        l.find_root_collar_psi();
        const Kind kind = l.operating_point_kind();
        if (kind != Kind::Interior && kind != Kind::BoundarySoil && kind != Kind::HydraulicShutdown) continue;
        const std::string at = " (" + std::string(l.operating_point_kind_name(kind)) + ", " +
                               std::to_string(layers) + " layers, psi " + std::to_string(psi0) +
                               ", PPFD " + std::to_string(ppfd) + ")";
        ++seen[l.operating_point_kind_name(kind)];
        const double row = profit_row_at_current_point(l, phylloptim::par_g_min);
        // The difference, swept in step and read at its floor (hazard 15).
        double best = std::numeric_limits<double>::infinity(), fd_best = 0.0;
        for (double rel : {1e-2, 3e-3, 1e-3, 3e-4}) {
          const double h = rel * kGmin;
          phylloptim::Leaf up = make(d, psi_soil, depth, kGmin + h), dn = make(d, psi_soil, depth, kGmin - h);
          up.find_root_collar_psi(); dn.find_root_collar_psi();
          if (up.operating_point_kind() != kind || dn.operating_point_kind() != kind) continue;
          const double fd = (up.profit_ - dn.profit_) / (2 * h);
          const double err = std::abs(fd - row) / std::max(std::abs(row), 1e-3);
          if (err < best) { best = err; fd_best = fd; }
        }
        if (!std::isfinite(best)) continue;
        ok(best < 1e-4, "dprofit/dg_min row matches a difference of the solve" + at + ": " +
                            std::to_string(row) + " vs " + std::to_string(fd_best));
        // Past the end of the stem curve no water moves, so g_min costs nothing.
        ok(row < 0.0 || kind == Kind::HydraulicShutdown,
           "and the cuticle costs the leaf carbon" + at);
        ++compared;
      }
    }
  }
  std::string names;
  for (const auto &[k, n] : seen) names += k + "(" + std::to_string(n) + ") ";
  printf("    %d points compared | kinds: %s\n", compared, names.c_str());
  ok(seen.count("interior") > 0 && seen.count("hydraulic-shutdown") > 0,
     "the sweep reached an interior point and a shut one");
}

// A shut leaf with a residual conductance sits on its wet bound: the collar
// where the soil delivers the cuticle's loss, the stem where it carries it,
// the stomata shut, and every output the evaluator's there (#182).
void test_shut_leaf_sits_on_its_wet_bound() {
  printf("residual conductance: a shut leaf sits on its wet bound and loses the residual flux\n");
  Drivers d;
  for (int layers : {1, 3}) {
    for (double psi0 : {6.0, 7.0}) {
      std::vector<double> psi_soil, depth;
      for (int i = 0; i < layers; ++i) { psi_soil.push_back(psi0 + 0.3 * i); depth.push_back(i + 1.0); }
      phylloptim::Leaf l = make(d, psi_soil, depth, kGmin);
      l.find_root_collar_psi();
      const std::string at = " (" + std::to_string(layers) + " layers, psi " + std::to_string(psi0) + ")";
      // Its stomata stay shut: at the wet bound, where a leaf with a cuticle sits
      // as an ordinary point (boundary-soil), or shut outright.
      ok(l.operating_point_kind() == Kind::BoundarySoil ||
             l.operating_point_kind() == Kind::HydraulicShutdown,
         std::string("the soil keeps the stomata shut, got ") +
             l.operating_point_kind_name(l.operating_point_kind()) + at);
      ok(l.residual_flux_ > 0.0 && l.transpiration_ == l.residual_flux_ + (l.transpiration_ - l.residual_flux_),
         "it still loses the residual flux" + at);
      ok(std::abs(l.transpiration_ - l.residual_flux_) <= 1e-11 * l.transpiration_,
         "with its stomata shut" + at);
      near(l.E_up_ / l.transpiration_, 1.0, 1e-9, "drawn from the soil at the collar" + at);
      ok(l.opt_psi_leaf_ > l.opt_root_psi_, "the stem is drier than the collar: water moves" + at);
      ok(std::isfinite(l.ci_) &&
             std::abs(l.stom_cond_CO2_ / l.cuticular_cond_CO2() - 1.0) < 1e-9,
         "and CO2 moves through the cuticle alone" + at);
      near(l.profit_, l.assim_colimited_ - l.hydraulic_cost_TF24(l.opt_psi_leaf_), 1e-12,
           "profit == A - C at the stem" + at);
    }
  }
  // Continuous into shutdown as the soil dries past the end of the stem curve.
  double prev = NAN, worst = 0.0;
  int shut = 0, open = 0;
  for (int i = 0; i <= 120; ++i) {
    const double ps = 12.0 + 0.035 * i;
    phylloptim::Leaf l = make(d, {ps}, {1.0}, kGmin);
    l.find_root_collar_psi();
    if (l.operating_point_kind() == Kind::HydraulicShutdown) ++shut; else ++open;
    if (!std::isnan(prev)) worst = std::max(worst, std::abs(l.profit_ - prev));
    prev = l.profit_;
  }
  ok(shut > 0 && open > 0, "the drying sweep reaches shutdown");
  ok(worst < 2e-1, "and profit moves by less than 2e-1 per 0.035 MPa across it");
}

// The costs that read the flux charge the TOTAL (#182).
void test_costs_charge_the_total_flux() {
  printf("residual conductance: the flux-reading costs charge the total flux\n");
  Drivers d;
  for (double ps : {1.0, 3.0}) {
    phylloptim::Leaf cf = make(d, {ps}, {1.0}, kGmin, false, true);
    cf.CF77_lambda_ = 1.0e5;
    cf.set_model("CF77", "leaf");
    cf.optimise();
    ok(cf.residual_flux_ > 0.0, "CF77: the residual flux is on");
    ok(cf.hydraulic_cost_ == cf.CF77_lambda_ * cf.transpiration_,
       "CF77 charges lambda times the TOTAL flux, bit for bit");
    near(cf.lambda_emergent_, cf.CF77_lambda_, 1e-12, "and its emergent lambda is its price");

    phylloptim::Leaf fl = make(d, {ps}, {1.0}, kGmin, false, true);
    fl.TF24_floor_lambda_o = 1.5e4;
    fl.set_model("TF24_floor", "leaf");
    fl.optimise();
    ok(fl.shadow_cost() == fl.TF24_floor_lambda_o * fl.transpiration_,
       "TF24_floor's shadow cost is the price times the TOTAL flux");
    near(fl.hydraulic_cost_, fl.hydraulic_cost_TF24_kernel(fl.opt_psi_leaf_) + fl.shadow_cost(), 1e-12,
         "and its cost is TF24's plus that");

    phylloptim::Leaf lc = make(d, {ps}, {1.0}, kGmin, false, true);
    lc.LeastCost_beta = 146.0;
    lc.set_model("LeastCost", "leaf");
    lc.optimise();
    near(lc.leastcost_flux(lc.opt_psi_leaf_, ps),
         lc.transpiration_ * phylloptim::kg_to_mol_h2o + 146.0 * lc.vcmax_ * phylloptim::umol_to_mol,
         1e-12, "LeastCost weighs the TOTAL flux");
  }
}

void test_closed_form_refuses() {
  printf("residual conductance: the closed form refuses it\n");
  Drivers d;
  phylloptim::Leaf l = make(d, {1.0}, {1.0}, kGmin, false, true);
  l.set_model("TF24", "leaf", "closed");
  const std::string msg = message_of([&] { l.optimise(); });
  ok(msg.find("residual conductance") != std::string::npos,
     "the closed form names the residual conductance in its refusal: " + msg);
  phylloptim::Leaf z = make(d, {1.0}, {1.0}, 0.0, false, true);
  z.set_model("TF24", "leaf", "closed");
  ok(message_of([&] { z.optimise(); }).empty(), "and runs at g_min = 0");
}

// The one corner the solve keeps, E_s >= 0 (#183): the value function is C1
// there and the argmax is not. Measured at the light-compensation corner by
// refining the step: profit's one-sided slopes converge, A's do not.
void test_light_compensation_corner() {
  printf("residual conductance: profit is C1 at the light-compensation corner, A is not\n");
  Drivers d;
  auto solve = [&](double ppfd) {
    Drivers di = d;
    di.PPFD = ppfd;
    phylloptim::Leaf l = make(di, {1.0}, {1.0}, kGmin);
    l.find_root_collar_psi();
    return l;
  };
  double lo = 20.0, hi = 35.0;
  for (int i = 0; i < 48; ++i) {
    const double m = 0.5 * (lo + hi);
    if (solve(m).operating_point_kind() == Kind::Interior) hi = m; else lo = m;
  }
  const double qc = 0.5 * (lo + hi);
  ok(solve(qc - 1.0).operating_point_kind() == Kind::BoundarySoil &&
         solve(qc + 1.0).operating_point_kind() == Kind::Interior,
     "the corner separates a wet-pinned leaf from an open one");
  double jump_coarse = NAN, jump_fine = NAN, a_jump_fine = NAN, a_slope = NAN;
  for (double h : {1.0, 0.01}) {
    phylloptim::Leaf a = solve(qc - 2 * h), b = solve(qc - h), c = solve(qc + h), e = solve(qc + 2 * h);
    const double sl = (b.profit_ - a.profit_) / h, sr = (e.profit_ - c.profit_) / h;
    const double al = (b.assim_colimited_ - a.assim_colimited_) / h,
                 ar = (e.assim_colimited_ - c.assim_colimited_) / h;
    if (h == 1.0) jump_coarse = std::abs(sr - sl);
    else { jump_fine = std::abs(sr - sl); a_jump_fine = std::abs(ar - al); a_slope = std::max(std::abs(al), std::abs(ar)); }
    printf("    step %.2f: profit slope %.5f | %.5f   A slope %.5f | %.5f\n", h, sl, sr, al, ar);
  }
  ok(jump_fine < 0.1 * jump_coarse,
     "profit's slope jump shrinks with the step: the value function is C1 at the corner");
  ok(a_jump_fine > 0.5 * a_slope,
     "A's one-sided slopes stay apart: the argmax has a corner there, as #183 allows");
}

// The acceptance sweep of #183, which every PR that touches the solve runs
// (`make -C tests/cpp smoothness`): second differences of the solved outputs
// and of two derivative rows across psi_soil, PPFD, VPD and Tair, with the
// energy balance off and on, 1 and 3 layers, g_min zero and positive, and the
// thermal-damage layer on. A slope jump is read as a fraction of the sweep's
// largest slope, and asserted below 0.05 everywhere except across the one
// corner: a triple of points that straddles a change of kind, which on these
// sweeps is only the stomatal corner E_s >= 0, where a leaf at its wet bound
// meets an open one. How many triples reached the assertion is counted and
// printed (hazard 16).
struct SweepPoint {
  double profit, A, E, psi, lambda, row_vcmax, row_gmin;
  std::string kind;
};

// A slope jump at the sweep's step is either a kink or curvature the step
// cannot resolve, and the two are told apart by refining: a kink's jump is
// independent of the step, curvature's falls with it. Every triple flagged at
// the sweep's step is re-measured at a hundredth of it, and the refined jump
// is what is asserted on, against the sweep's own scale.
// Totals over every call, for the summary line.
long g_sweep_asserted = 0, g_sweep_corner = 0;

void sweep_kinks(const std::vector<double> &x, const std::vector<SweepPoint> &o,
                 double SweepPoint::*f, const char *name, const std::string &where,
                 double limit, const std::function<SweepPoint(double)> &at) {
  std::vector<double> d1;
  for (size_t i = 1; i < x.size(); ++i) d1.push_back(((o[i].*f) - (o[i - 1].*f)) / (x[i] - x[i - 1]));
  double mx = 0.0;
  for (double v : d1) if (std::isfinite(v)) mx = std::max(mx, std::abs(v));
  double worst = 0.0, worst_refined = 0.0;
  size_t at_i = 0, at_r = 0;
  int excluded = 0, refined = 0, asserted = 0;
  for (size_t i = 1; i < d1.size(); ++i) {
    // the triple (i-1, i, i+1) of points
    const std::string &k0 = o[i - 1].kind, &k1 = o[i].kind, &k2 = o[i + 1].kind;
    const bool corner = k0 != k1 || k1 != k2;
    if (corner) { ++excluded; continue; }
    const double k = std::abs(d1[i] - d1[i - 1]) / mx;
    if (!std::isfinite(k)) continue;
    ++asserted;
    if (k > worst) { worst = k; at_i = i; }
    if (k <= limit) continue;
    // Refine: the same jump at a hundredth of the step, about the same point.
    const double h = (x[i] - x[i - 1]) / 100.0;
    const SweepPoint a = at(x[i] - h), c = at(x[i] + h);
    if (a.kind != k1 || c.kind != k1) continue;  // the refinement crossed a corner
    const double s1 = ((o[i].*f) - (a.*f)) / h, s2 = ((c.*f) - (o[i].*f)) / h;
    const double kr = std::abs(s2 - s1) / mx;
    ++refined;
    if (std::isfinite(kr) && kr > worst_refined) { worst_refined = kr; at_r = i; }
  }
  // A discontinuity is not a kink, and the refinement above cannot see one: it
  // is centred on a sweep point, and a step between two points straddles it.
  // So every interval whose change dwarfs both neighbours' is bisected. A jump
  // keeps its size as the interval shrinks; a steep smooth stretch does not.
  // Measured against the output's range over the sweep, away from corners.
  double lo_v = std::numeric_limits<double>::infinity(), hi_v = -lo_v;
  for (const SweepPoint &p : o) {
    const double v = p.*f;
    if (std::isfinite(v)) { lo_v = std::min(lo_v, v); hi_v = std::max(hi_v, v); }
  }
  const double range = hi_v - lo_v;
  double worst_jump = 0.0, jump_at = x.front();
  int bisected = 0;
  for (size_t i = 0; i + 1 < x.size(); ++i) {
    if (o[i].kind != o[i + 1].kind) continue;  // a corner
    const double fa0 = o[i].*f, fb0 = o[i + 1].*f;
    const double df = std::abs(fb0 - fa0);
    const double prev = i > 0 ? std::abs(fa0 - (o[i - 1].*f)) : 0.0;
    const double next = i + 2 < x.size() ? std::abs((o[i + 2].*f) - fb0) : 0.0;
    if (!std::isfinite(df) || !(range > 0.0) || df <= 1e-6 * range ||
        df <= 4.0 * std::max(prev, next)) {
      continue;
    }
    ++bisected;
    double a = x[i], b = x[i + 1], fa = fa0, fb = fb0;
    bool crossed = false;
    for (int k = 0; k < 30; ++k) {
      const double m = 0.5 * (a + b);
      const SweepPoint pm = at(m);
      if (pm.kind != o[i].kind) { crossed = true; break; }
      const double fm = pm.*f;
      if (std::abs(fm - fa) > std::abs(fb - fm)) { b = m; fb = fm; } else { a = m; fa = fm; }
    }
    if (crossed) continue;  // a corner inside the interval
    const double jump = std::abs(fb - fa) / range;
    if (jump > worst_jump) { worst_jump = jump; jump_at = a; }
  }
  printf("    %-28s %-8s %d triples asserted; worst at the step %.3e of scale (x = %.3f, %s); %d refined, worst %.3e (x = %.3f); %d corner triples excluded; %d bisected, worst jump %.3e of range (x = %.6f)\n",
         where.c_str(), name, asserted, worst, x[at_i], o[at_i].kind.c_str(), refined, worst_refined,
         x[at_r], excluded, bisected, worst_jump, jump_at);
  g_sweep_asserted += asserted;
  g_sweep_corner += excluded;
  ok(asserted > 0, std::string(name) + " reached the assertion on at least one triple, " + where);
  ok(worst_refined < limit, std::string(name) + " has no slope jump above " + std::to_string(limit) +
                                " of scale away from a corner once the step is refined, " + where);
  ok(worst_jump < 1e-6, std::string(name) + " has no jump away from a corner, " + where);
}

void test_smoothness_sweep() {
  printf("residual conductance: the #183 second-difference sweep\n");
  struct Config { bool eb; int layers; bool damage; double g_min; };
  const Config configs[] = {{false, 1, false, kGmin}, {false, 3, false, kGmin},
                            {true, 1, false, kGmin},  {true, 3, false, kGmin},
                            {false, 1, true, kGmin},  {true, 1, true, kGmin},
                            {false, 1, false, 0.0},   {true, 1, false, 0.0}};
  g_sweep_asserted = g_sweep_corner = 0;
  std::map<std::string, long> kinds;
  for (const Config &c : configs) {
    auto solve = [&](double psi0, double ppfd, double vpd, double tair) {
      Drivers d;
      d.PPFD = ppfd; d.atm_vpd = vpd; d.leaf_temp = tair;
      std::vector<double> psi_soil, depth;
      for (int i = 0; i < c.layers; ++i) { psi_soil.push_back(psi0 + 0.25 * i); depth.push_back(i + 1.0); }
      phylloptim::Leaf l = make(d, psi_soil, depth, c.g_min, c.eb, false,
                                std::numeric_limits<double>::quiet_NaN(),
                                c.damage ? +[](phylloptim::Leaf &x) { x.use_thermal_damage_ = true; set_atls_inputs(x); }
                                         : static_cast<void (*)(phylloptim::Leaf &)>(nullptr));
      l.find_root_collar_psi();
      SweepPoint p{l.profit_, l.assim_colimited_, l.transpiration_, l.opt_psi_leaf_,
                   l.lambda_emergent_, NAN, NAN, l.operating_point_kind_name(l.operating_point_kind())};
      ++kinds[p.kind];
      // The derivative rows, off the thermal-damage layer and the energy balance,
      // which the surface refuses, and with a residual conductance to differentiate.
      if (!c.damage && !c.eb && c.g_min > 0.0 &&
          (l.operating_point_kind() == Kind::Interior || l.operating_point_kind() == Kind::BoundarySoil ||
                        l.operating_point_kind() == Kind::HydraulicShutdown)) {
        p.row_vcmax = profit_row_at_current_point(l, phylloptim::par_vcmax_25);
        p.row_gmin = profit_row_at_current_point(l, phylloptim::par_g_min);
      }
      return p;
    };
    const std::string where = std::string(c.eb ? "EB on, " : "EB off, ") + std::to_string(c.layers) +
                              (c.damage ? " layer, damage on" : " layers") +
                              (c.g_min > 0.0 ? "" : ", g_min 0");
    auto run = [&](const char *driver, double from, double to, double step,
                   std::function<SweepPoint(double)> at) {
      std::vector<double> x;
      std::vector<SweepPoint> o;
      for (double v = from; v <= to + 1e-9; v += step) { x.push_back(v); o.push_back(at(v)); }
      const std::string w = where + ", " + driver;
      sweep_kinks(x, o, &SweepPoint::profit, "profit", w, 0.05, at);
      sweep_kinks(x, o, &SweepPoint::A, "A", w, 0.05, at);
      sweep_kinks(x, o, &SweepPoint::E, "E", w, 0.05, at);
      sweep_kinks(x, o, &SweepPoint::psi, "psi_leaf", w, 0.05, at);
      sweep_kinks(x, o, &SweepPoint::lambda, "lambda", w, 0.05, at);
      if (!c.damage && !c.eb && c.g_min > 0.0) {
        sweep_kinks(x, o, &SweepPoint::row_vcmax, "d/dvcmax", w, 0.05, at);
        sweep_kinks(x, o, &SweepPoint::row_gmin, "d/dg_min", w, 0.05, at);
      }
    };
    run("psi_soil", 0.3, 6.4, 0.01, [&](double v) { return solve(v, 900.0, 2.0, 25.0); });
    run("PPFD", 0.0, 300.0, 1.0, [&](double v) { return solve(1.0, v, 2.0, 25.0); });
    // The deficit stays below saturation at the air temperature, past which it
    // is refused: to 3.15 kPa at 25 C (e_s 3.17), and to 6 kPa at 40 C (7.38),
    // which also crosses the emissivity's soft bound with the balance on.
    run("VPD", 0.05, 3.15, 0.02, [&](double v) { return solve(1.0, 900.0, v, 25.0); });
    run("VPD at 40 C", 0.05, 6.0, 0.02, [&](double v) { return solve(1.0, 900.0, v, 40.0); });
    run("Tair", 5.0, 48.0, 0.25, [&](double v) { return solve(1.0, 900.0, 0.5, v); });
  }
  printf("    in all: %ld triples asserted, %ld corner triples excluded; kinds solved:",
         g_sweep_asserted, g_sweep_corner);
  for (const auto &kv : kinds) printf(" %s %ld", kv.first.c_str(), kv.second);
  printf("\n");
}

// The derivative surface follows the energy balance on every curve: ProfitMax,
// the last to, answers with it on, and the supply draw reads no temperature.
void test_surface_refuses_the_energy_balance() {
  printf("the derivative surface follows the energy balance on every curve\n");
  using CC = phylloptim::Leaf::CostCurve;
  Drivers d;
  d.PPFD = 1500.0;
  const std::vector<double> psi{1.0, 1.25}, depth{1.0, 2.0};
  phylloptim::Leaf on = make(d, psi, depth, kGmin, true);
  on.find_root_collar_psi();
  ok(on.operating_point_kind() == Kind::Interior,
     "the gate-on fixture is an interior point");
  ok(std::isfinite(profit_row_at_current_point(on, phylloptim::par_vcmax_25)),
     "TF24 answers with the energy balance on");
  phylloptim::Leaf pm = make(d, psi, depth, kGmin, true);
  pm.set_model(CC::ProfitMax, true);
  pm.optimise();
  const phylloptim::leaf_pars<double> pars = pm.passive_pars();
  const auto pm_supply = pm.roots_.held_supply();
  const auto draw = pm.supply_draw_at<double>(pm.opt_root_psi_, pm_supply);
  ok(std::isfinite(draw.E_up.value), "supply_draw_at reads no temperature, and runs");
  const auto out = pm.billed_outputs_at<CC::ProfitMax, double>(pm.opt_root_psi_, draw, pars);
  ok(out.profit == pm.profit_, "and ProfitMax's profit is the solve's with it on");
}

// #184's first requirement from plant: a leaf forced shut -- by a dry soil, the
// dark or heat -- is an ordinary state, so it comes back as a value and never as
// an exception. plant reaches it on ~4% of TF24f leaf solves at 1 m/yr of rain
// and 50% at 0.04 m/yr. Refusals are for drivers outside the model's range, and
// none is offered here: each air temperature takes a deficit below saturation.
// Swept over soil potentials from wet to plant's 1000 MPa ceiling, with and
// without a residual conductance, the energy balance off and on, on the collar
// route at one and three layers and on every leaf route.
void test_a_shut_leaf_is_a_value_not_an_exception() {
  printf("residual conductance: a dry, dark or hot leaf is a value, never an exception (#184)\n");
  const double soils[] = {0.05, 0.5, 3.0, 6.0, 10.0, 30.0, 100.0, 1000.0};
  struct Air { double tair, vpd; };
  const Air airs[] = {{10.0, 0.8}, {30.0, 2.0}, {45.0, 3.0}};
  const double lights[] = {0.0, 1500.0};
  const double gmins[] = {0.0, kGmin, 2e-3};
  int solved = 0, threw = 0, non_finite = 0, shut = 0;
  std::string first_failure;
  auto record = [&](const phylloptim::Leaf &l, const std::string &where) {
    const bool finite = std::isfinite(l.profit_) && std::isfinite(l.assim_colimited_) &&
                        std::isfinite(l.transpiration_) && std::isfinite(l.opt_psi_leaf_) &&
                        std::isfinite(l.residual_flux_) && !std::isnan(l.ci_);
    if (!finite) {
      ++non_finite;
      if (first_failure.empty()) first_failure = "non-finite output at " + where;
    }
    ++solved;
  };
  for (double g : gmins) {
    for (bool eb : {false, true}) {
      for (double ppfd : lights) {
        for (const Air &a : airs) {
          Drivers d;
          d.PPFD = ppfd; d.atm_vpd = a.vpd; d.leaf_temp = a.tair;
          for (double psi : soils) {
            char buf[160];
            snprintf(buf, sizeof buf, "g_min %g EB %d PPFD %g Tair %g psi_soil %g", g, eb, ppfd, a.tair, psi);
            for (int layers : {1, 3}) {
              std::vector<double> ps, dp;
              for (int i = 0; i < layers; ++i) { ps.push_back(psi + 0.25 * i); dp.push_back(i + 1.0); }
              const std::string where = std::string(buf) + ", collar, " + std::to_string(layers) + " layers";
              try {
                phylloptim::Leaf l = make(d, ps, dp, g, eb);
                l.find_root_collar_psi();
                if (l.operating_point_kind() == Kind::HydraulicShutdown) ++shut;
                record(l, where);
              } catch (const std::exception &e) {
                ++threw;
                if (first_failure.empty()) first_failure = where + ": " + e.what();
              }
            }
            for (int curve = 0; curve < phylloptim::Leaf::n_cost_curves; ++curve) {
              const auto K = static_cast<phylloptim::Leaf::CostCurve>(curve);
              const std::string where = std::string(buf) + ", leaf route " +
                                        phylloptim::Leaf::curve_name(curve);
              try {
                phylloptim::Leaf l = make(d, {psi}, {1.0}, g, eb, true);
                l.CF77_lambda_ = 1.5e5;
                l.TF24_floor_lambda_o = kFloorLambdaO;
                l.LeastCost_beta = 146.0;
                l.set_model(K, false);
                l.optimise();
                record(l, where);
              } catch (const std::exception &e) {
                ++threw;
                if (first_failure.empty()) first_failure = where + ": " + e.what();
              }
            }
          }
        }
      }
    }
  }
  printf("    %d solves (%d collar solves shut), %d threw, %d with a non-finite output%s%s\n",
         solved + threw, shut, threw, non_finite, first_failure.empty() ? "" : "; first: ",
         first_failure.c_str());
  ok(threw == 0, "no dry, dark or hot leaf throws, on any route");
  ok(non_finite == 0, "every such leaf reports finite outputs");
  ok(shut > 0, "the sweep reaches the shut state");
}

// #184's second requirement: where the soil cannot supply what the cuticle
// loses, the leaf is smooth in psi_soil -- soft caps, not a min(). The sweep
// runs from where the caps are inert to where they bind hard, which the
// realised fraction r*/r reports, at a conductance large enough to get there,
// and on past the end of the root curve to 17 MPa.
//
// Two ceilings bind in turn: the stem's (what it can carry from the collar),
// then the roots' (what they can deliver from this soil, which falls to zero
// at the end of the root curve, where G is G(inf) to rounding; #85). Past it no
// layer conducts and the shut leaf sits at the wettest layer's gravity balance. The roots' ceiling
// shrinks as U (U / (U + r))^2 so the leaf reaches that balance with the
// value AND the slope it has on the far side (residual_root_ceiling_at says
// why); the cap is asserted on both, directly, as well as swept through.
void test_the_capped_regime_is_smooth_in_the_soil() {
  printf("residual conductance: smooth in psi_soil where the soil cannot supply the cuticle (#184)\n");
  for (bool eb : {false, true}) {
    for (int layers : {1, 3}) {
      for (double g : {5e-4, 5e-3}) {
        auto solve = [&](double psi0) {
          Drivers d;
          d.PPFD = 900.0; d.atm_vpd = 2.0; d.leaf_temp = 25.0;
          std::vector<double> ps, dp;
          for (int i = 0; i < layers; ++i) { ps.push_back(psi0 + 0.25 * i); dp.push_back(i + 1.0); }
          phylloptim::Leaf l = make(d, ps, dp, g, eb);
          l.find_root_collar_psi();
          return l;
        };
        auto at = [&](double psi0) {
          const phylloptim::Leaf l = solve(psi0);
          return SweepPoint{l.profit_, l.assim_colimited_, l.transpiration_, l.opt_psi_leaf_,
                            l.lambda_emergent_, NAN, NAN,
                            l.operating_point_kind_name(l.operating_point_kind())};
        };
        char buf[64];
        snprintf(buf, sizeof buf, ", g_min %g", g);
        const std::string where = std::string(eb ? "EB on, " : "EB off, ") + std::to_string(layers) +
                                  " layers" + buf;
        std::vector<double> x;
        std::vector<SweepPoint> o;
        double frac_first = NAN, frac_min = 1.0, frac_last = NAN;
        for (double v = 4.0; v <= 17.0 + 1e-9; v += 0.01) {
          x.push_back(v);
          o.push_back(at(v));
          const phylloptim::Leaf l = solve(v);
          const double frac = l.residual_flux_ / l.residual_demand_at(l.Tleaf_);
          if (std::isnan(frac_first)) frac_first = frac;
          frac_min = std::min(frac_min, frac);
          frac_last = frac;
        }
        printf("    %s: realised r*/r runs %.4f -> %.4f\n", where.c_str(), frac_first, frac_last);
        ok(frac_first > 0.99 && frac_last == 0.0,
           "the sweep runs from an inert cap to a soil that supplies no residual flux, " + where);
        sweep_kinks(x, o, &SweepPoint::profit, "profit", where, 0.05, at);
        sweep_kinks(x, o, &SweepPoint::A, "A", where, 0.05, at);
        sweep_kinks(x, o, &SweepPoint::E, "E", where, 0.05, at);
        // psi_leaf to the end of the stem curve: with the tail carried as
        // G(inf) - G (#190) a shut leaf's leaf potential is resolved to there,
        // where before ~12 MPa it was defined only to G(inf)'s rounding. Within
        // a gravity head of the end itself the shut collar is held AT the end,
        // where what the stem has left steps from e^-40 of G(inf) to zero,
        // until the wettest layer's balance passes it: a corner the curve's
        // end sets, at fluxes of ~1e-50.
        {
          const phylloptim::Leaf l0 = solve(4.0);
          const double end = phylloptim::vulnerability_psi_end(l0.stem_b, l0.stem_c);
          std::vector<double> xs;
          std::vector<SweepPoint> os;
          for (std::size_t i = 0; i < x.size(); ++i) {
            if (x[i] < end - 0.1) { xs.push_back(x[i]); os.push_back(o[i]); }
          }
          sweep_kinks(xs, os, &SweepPoint::psi, "psi_leaf", where, 0.05, at);
        }

        // The end of the root curve: where the wettest layer stops conducting,
        // bisected.
        double open = 12.0, closed = 16.0;
        for (int k = 0; k < 60; ++k) {
          const double m = 0.5 * (open + closed);
          bool conducts = false;
          const phylloptim::Leaf l = solve(m);
          l.roots_.wet_floor_layer(l.roots_.psi_soil_, conducts);
          (conducts ? open : closed) = m;
        }
        const phylloptim::Leaf wet = solve(open), dry = solve(closed);
        const double step_leaf = std::abs(dry.opt_psi_leaf_ - wet.opt_psi_leaf_);
        // One-sided slopes a few uptake bands (1e-8) away, at a step where the
        // difference is the slope: measured converged from 1e-5 down.
        const double h = 1e-5;
        auto slope = [&](double from, double sign, double SweepPoint::*f) {
          return sign * ((at(from + sign * 2 * h).*f) - (at(from + sign * h).*f)) / h;
        };
        const double s_wet = slope(open, -1.0, &SweepPoint::psi), s_dry = slope(closed, 1.0, &SweepPoint::psi);
        const double p_wet = slope(open, -1.0, &SweepPoint::profit), p_dry = slope(closed, 1.0, &SweepPoint::profit);
        printf("    %s: at the root curve's end (psi_soil %.6f) psi_leaf steps %.2e MPa; slopes psi_leaf %.5f | %.5f, profit %.5f | %.5f\n",
               where.c_str(), open, step_leaf, s_wet, s_dry, p_wet, p_dry);
        ok(step_leaf < 1e-7, "the shut leaf is continuous through the root curve's end, " + where);
        ok(std::abs(s_wet - s_dry) < 1e-3 * std::max(1.0, std::abs(s_dry)) &&
               std::abs(p_wet - p_dry) < 1e-3 * std::max(0.05, std::abs(p_dry)),
           "and its slope in psi_soil is continuous there too, " + where);
      }
    }
  }
}

// The rows where the roots' ceiling binds -- and, for contrast, at an open
// leaf in a wet soil: profit's and the residual flux's, in g_min, in the root
// curve that sets the ceiling and in the stem curve that carries the flux,
// each against a difference of the solve. A row the ceiling's rows were missing from would
// pass everywhere the ceiling is inert, which is everywhere else in this file.
static double g_perturbed_root_P50 = NAN, g_perturbed_stem_P50 = NAN;

void test_rows_where_the_roots_cap_binds() {
  printf("residual conductance: rows where the roots' ceiling binds, against a difference of the solve (#184)\n");
  auto seat = [](phylloptim::Leaf &l) {
    l.set_traits(l.vcmax_25, l.stem_c,
                 std::isnan(g_perturbed_stem_P50) ? l.stem_P50 : g_perturbed_stem_P50,
                 l.roots_.root_c,
                 std::isnan(g_perturbed_root_P50) ? l.roots_.root_P50 : g_perturbed_root_P50,
                 l.TF24_beta2, l.jmax_25, l.phi_J, l.curv_fact_elec_trans, l.curv_fact_colim,
                 l.TF24_cost_scale, l.R_d_25, l.JS22_gamma, l.CMax_a, l.CMax_b, l.g_min);
  };
  int compared = 0;
  double worst = 0.0, worst_r = 0.0;
  std::string worst_r_at;
  for (double psi0 : {2.0, 7.2, 7.35, 7.385}) {
    for (int layers : {1, 3}) {
      Drivers d;
      d.PPFD = 900.0; d.atm_vpd = 2.0;
      std::vector<double> ps, dp;
      for (int i = 0; i < layers; ++i) { ps.push_back(psi0 + 0.25 * i); dp.push_back(i + 1.0); }
      auto solved = [&](double g) {
        phylloptim::Leaf l = make(d, ps, dp, g, false, false,
                                  std::numeric_limits<double>::quiet_NaN(), seat);
        l.find_root_collar_psi();
        return l;
      };
      g_perturbed_root_P50 = g_perturbed_stem_P50 = NAN;
      const phylloptim::Leaf base = solved(5e-3);
      const double r = base.residual_demand_at(base.Tleaf_);
      const double frac = base.residual_flux_ / r;
      const double root_share = base.residual_root_ceiling_;
      struct Par { int index; const char *name; double value; double *knob; bool is_gmin; };
      const Par ps_[] = {{phylloptim::par_g_min, "g_min", 5e-3, nullptr, true},
                         {phylloptim::par_root_P50, "root_P50", base.roots_.root_P50, &g_perturbed_root_P50, false},
                         {phylloptim::par_stem_P50, "stem_P50", base.stem_P50, &g_perturbed_stem_P50, false}};
      for (const Par &p : ps_) {
        phylloptim::Leaf l = solved(5e-3);
        double residual_row = NAN;
        const double row = profit_row_at_current_point(l, p.index, &residual_row);
        double best = std::numeric_limits<double>::infinity(), fd_best = NAN;
        double best_r = std::numeric_limits<double>::infinity(), fd_r_best = NAN;
        // Swept over five decades and read at the floor (hazard 15): near the
        // roots' cap profit is strongly curved in the root curve, and a step of
        // 1e-4 of root_P50 still carries 5e-4 of truncation where 1e-6 is exact
        // to 5e-8.
        for (double rel : {1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7}) {
          const double h = rel * p.value;
          double up, dn, r_up, r_dn;
          if (p.is_gmin) {
            const phylloptim::Leaf a = solved(p.value + h), b = solved(p.value - h);
            up = a.profit_; dn = b.profit_; r_up = a.residual_flux_; r_dn = b.residual_flux_;
          } else {
            *p.knob = p.value + h; const phylloptim::Leaf a = solved(5e-3);
            *p.knob = p.value - h; const phylloptim::Leaf b = solved(5e-3);
            *p.knob = NAN;
            up = a.profit_; dn = b.profit_; r_up = a.residual_flux_; r_dn = b.residual_flux_;
          }
          // The residual flux's row against its own difference, with a floor
          // scaled to the flux itself (kg m^-2 s^-1, so no fixed absolute floor
          // serves). A row below 1e-3 of the flux per unit trait is one whose
          // components nearly cancel, and past the root table's last knot they
          // carry the table's fit error (#147, ~1e-4 of each, #184 section 2):
          // measured 5e-5 of a 3.4e-12 row at a 2.9e-8 flux, flat in the step
          // from 3e-5 to 3e-6, so not noise and not this flux's own algebra. A
          // floor of 1e-6 failed on Linux at a 7e-17 row; 1e-4 left 2x margin.
          const double fd_r = (r_up - r_dn) / (2 * h);
          const double err_r = std::abs(fd_r - residual_row) /
                               std::max(std::abs(residual_row), 1e-3 * l.residual_flux_ / p.value);
          if (err_r < best_r) { best_r = err_r; fd_r_best = fd_r; }
          const double fd = (up - dn) / (2 * h);
          // Relative, with a floor of 1e-3 profit per unit trait: a row of a
          // millionth of the profit's scale sits on the solve's noise floor
          // divided by the step, so there it is held to 1e-7 absolute.
          const double err = std::abs(fd - row) / std::max(std::abs(row), 1e-3);
          if (err < best) { best = err; fd_best = fd; }
        }
        worst = std::max(worst, best);
        ++compared;
        char buf[200];
        snprintf(buf, sizeof buf, " (psi_soil %.3f, %d layers, r*/r %.3f, roots' ceiling %.3e): %.6e vs %.6e",
                 psi0, layers, frac, root_share, row, fd_best);
        ok(best < 1e-4, std::string("dprofit/d") + p.name + " matches a difference of the solve" + buf);
        snprintf(buf, sizeof buf, " (psi_soil %.3f, %d layers): %.6e vs %.6e", psi0, layers,
                 residual_row, fd_r_best);
        ok(best_r < 1e-4, std::string("the residual flux's row in ") + p.name +
                              " matches a difference of the solve" + buf);
        if (best_r > worst_r) {
          worst_r = best_r;
          worst_r_at = std::string(p.name) + buf;
        }
      }
    }
  }
  printf("    %d points compared, worst relative error at the difference's floor %.3e (profit), %.3e (residual flux, d/d%s)\n",
         compared, worst, worst_r, worst_r_at.c_str());
}

}  // namespace residual

void run_residual() {
  residual::test_smoothness_sweep();
  residual::test_temperature_factor();
  residual::test_cap_is_a_smooth_minimum();
  residual::test_reduces_to_nothing_at_zero();
  residual::test_water_budget_closes();
  residual::test_manuscript_figure_g0();
  residual::test_shift_identity();
  residual::test_rows_match_a_difference();
  residual::test_shut_leaf_sits_on_its_wet_bound();
  residual::test_costs_charge_the_total_flux();
  residual::test_closed_form_refuses();
  residual::test_light_compensation_corner();
  residual::test_surface_refuses_the_energy_balance();
  residual::test_a_shut_leaf_is_a_value_not_an_exception();
  residual::test_the_capped_regime_is_smooth_in_the_soil();
  residual::test_rows_where_the_roots_cap_binds();
}
