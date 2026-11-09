// The derivative surface: the rows handed to a consumer, against differences of the solve and against the solve itself.

#include "leaf_harness.hpp"

namespace {

// R_d's TEMPERATURE RESPONSE (#41).
//
// ⚠️ THE GOLDEN FILE IS NEARLY BLIND TO THIS, and that is the reason this test
// exists. Every reference value in the model is DEFINED at 25 C, so a change to any
// response curve is inert there by construction; the golden grid carries one hot
// block for exactly this reason, and this test is what pins the response itself.
// The pack-reading kernels against the member-reading ones they are twins of,
// at double, where the two must be THE SAME NUMBERS and not merely close.
//
// This is the only check on the reverse-mode surface's VALUES. The transpose
// identity reads that surface in two directions and is satisfied exactly by a
// kernel evaluated at the wrong temperature or in the wrong light, because both
// directions are then wrong together.
//
// ⚠️ THE ENERGY-BALANCE ARM IS THE POINT, not a second case for completeness.
// Off that path the block's temperature and `leaf_temp_` are one number and a
// kernel reading either passes; on it they run up to 10 K apart, and reading
// `leaf_temp_` moved assimilation by 100%, electron transport by 84% and dA/dci
// by 8.3%. The block members are compared FIRST so a failure names the cause
// rather than the three consequences below it.
void test_pack_kernels_are_the_models_own() {
  printf("the pack-reading kernels are the model's own numbers\n");
  using AD = xad::fwd<double>::active_type;
  int compared = 0, agreed = 0, refused = 0;
  double widest_gap = 0.0;
  // Identical numbers, and identical refusals. At a leaf 50 C hot the
  // colimitation discriminant leaves dA/dci NaN, and two NaNs from two spellings
  // of one function are the agreement being asserted rather than an exception to
  // it -- `near` cannot say that, since NaN compares false against everything.
  // Counted, so a fixture that started refusing everywhere cannot pass by it.
  auto same = [&](double got, double want, const std::string &what) {
    if (std::isnan(got) && std::isnan(want)) {
      ++refused;
      ++checks;
      return;
    }
    ++agreed;
    near(got, want, 0.0, what);
  };
  for (bool gate : {false, true}) {
    for (double psi0 : {0.5, 1.0, 2.0, 3.0}) {
      for (double ppfd : {100.0, 1500.0}) {
        for (double temp : {25.0, 40.0}) {
          Drivers d;
          d.PPFD = ppfd;
          d.leaf_temp = temp;
          // The thermal cost scales jmax, so it is on at the hot rows: the
          // scalar-temperature capacity must carry its factor's slope too.
          phylloptim::Leaf l =
              make_pm_leaf(d, {psi0, psi0 + 0.35, psi0 + 0.7}, {1.0, 2.0, 3.0},
                           gate, 2.0, 0.05, temp > 30.0);
          l.find_root_collar_psi();
          const double ci = l.ci_;
          if (!std::isfinite(ci)) continue;
          ++compared;
          widest_gap = std::max(widest_gap, std::abs(l.Tleaf_ - l.leaf_temp_));
          const phylloptim::leaf_pars<double> p = l.passive_pars();
          const std::string at = " at psi_soil=" + std::to_string(psi0) +
                                 " ppfd=" + std::to_string(ppfd) + " Tair=" +
                                 std::to_string(temp) +
                                 (gate ? " eb=1" : " eb=0");
          same(l.vcmax_at<double>(p), l.vcmax_, "vcmax" + at);
          same(l.jmax_at<double>(p), l.jmax_, "jmax" + at);
          same(l.respiration_at<double>(p), l.R_d_, "respiration" + at);
          same(l.electron_transport_at<double>(p), l.electron_transport(),
               "electron transport" + at);
          const phylloptim::Leaf::PhotoCapacity<double> cap =
              l.photo_capacity_at<double>(p);
          same(l.assim_colimited_kernel<double>(ci, cap), l.assim_colimited(ci),
               "colimited assimilation" + at);
          AD ci_ad = ci;
          xad::derivative(ci_ad) = 1.0;
          same(l.assim_slope_at<double>(ci, cap),
               xad::derivative(l.assim_colimited_kernel(ci_ad)),
               "dA/dci" + at);
          // The cost's slope walks its own four arguments rather than the pack,
          // so it needs the same referee the assimilation's has.
          AD sigma_ad = l.opt_psi_leaf_;
          xad::derivative(sigma_ad) = 1.0;
          same(l.cost_slope_at<double>(l.opt_psi_leaf_, p),
               xad::derivative(l.hydraulic_cost_TF24_kernel(sigma_ad)),
               "dC/dsigma" + at);
          // At the block's own temperature on a scalar, all six capacities are
          // the members bit for bit, Gamma* and Km included; and with that
          // temperature moving, their slopes are the model's own central
          // difference in temperature, to its truncation error.
          const phylloptim::Leaf::PhotoCapacity<double> at_T =
              l.photo_capacity_at<double>(p, l.photo_temp_);
          same(at_T.vcmax, l.vcmax_, "vcmax at a scalar temperature" + at);
          same(at_T.transport, l.electron_transport_,
               "electron transport at a scalar temperature" + at);
          same(at_T.respiration, l.R_d_, "respiration at a scalar temperature" + at);
          same(at_T.gamma, l.gamma_, "Gamma* at a scalar temperature" + at);
          same(at_T.km, l.km_, "Km at a scalar temperature" + at);
          {
            using TT = odelia::ode::tangent_scalar<double>;
            phylloptim::leaf_pars<TT> pt;
            for (std::size_t i = 0; i < pt.size(); ++i) pt[i] = TT(p[i]);
            TT temp_t(l.photo_temp_);
            odelia::ode::set_direction(temp_t, 1.0);
            const phylloptim::Leaf::PhotoCapacity<TT> c =
                l.photo_capacity_at<TT>(pt, temp_t);
            const auto slope_ok = [&](const TT& got, auto read, const char* what) {
              const double want = l.temp_deriv(l.photo_temp_, read);
              near(odelia::ode::derivative_along(got), want,
                   1e-6 * std::abs(want) + 1e-12,
                   std::string("d") + what + "/dT on a scalar temperature" + at);
            };
            slope_ok(c.vcmax, [&] { return l.vcmax_; }, "vcmax");
            slope_ok(c.transport, [&] { return l.electron_transport_; }, "J");
            slope_ok(c.respiration, [&] { return l.R_d_; }, "R_d");
            slope_ok(c.gamma, [&] { return l.gamma_; }, "Gamma*");
            slope_ok(c.km, [&] { return l.km_; }, "Km");
            if (l.use_thermal_cost_) {
              near(odelia::ode::derivative_along(l.thermal_cost_on<TT>(temp_t)),
                   l.thermal_cost_temp_deriv(l.photo_temp_),
                   1e-12 * std::abs(l.thermal_cost_temp_deriv(l.photo_temp_)) + 1e-15,
                   "the thermal cost's slope on a scalar temperature" + at);
              phylloptim::Leaf bad = l;
              bad.Tcrit_ = bad.T50_;
              ok(throws_with([&] { bad.thermal_cost_on<TT>(temp_t); }, "T50 > Tcrit"),
                 "and it refuses T50 <= Tcrit as thermal_cost_at does" + at);
            }
          }
        }
      }
    }
  }
  ok(compared >= 20, "the sweep reached at least 20 operating points");
  ok(agreed > 4 * refused,
     "most comparisons were of numbers, not of refusals: " +
         std::to_string(agreed) + " numbers, " + std::to_string(refused) +
         " refusals");
  // The claim is only evidence where the two temperatures actually differ, so
  // a fixture that stopped separating them would have stopped checking.
  ok(widest_gap > 1.0,
     "the sweep reached a leaf more than 1 K off air temperature, gap=" +
         std::to_string(widest_gap));
}

// Light reaches carbon with a ROW, and the row is the right size.
//
// ⚠️ AN EXACT ZERO IS THE FAILURE THIS EXISTS FOR, not a small number. The pack
// slot is the only route light has into the electron transport; read off the
// `PPFD_` member instead, every trait that moves the light a cohort stands in --
// the extinction coefficient, the leaf area a stem carries -- reaches
// assimilation with no row and the gradient reads zero rather than wrong. The
// transpose identity is satisfied exactly by a zero row on both sides, so it
// cannot see this either.
//
// The difference is taken at a FIXED ci, so it differences the kernel and not the
// solve, and the row is required to FALL with light: at 100 umol the leaf is
// light-limited and at 1500 it is rubisco-limited, so a row that did not saturate
// would be reading light somewhere it does not belong.
void test_light_reaches_carbon_with_a_row() {
  printf("light reaches carbon with a row\n");
  using AD = xad::fwd<double>::active_type;
  double previous = std::numeric_limits<double>::infinity();
  int compared = 0;
  for (double ppfd : {100.0, 400.0, 900.0, 1500.0}) {
    Drivers d;
    d.PPFD = ppfd;
    phylloptim::Leaf l = make_leaf(d, {1.0}, {1.0});
    l.find_root_collar_psi();
    const double ci = l.ci_;
    if (!std::isfinite(ci)) continue;
    ++compared;
    const std::string at = " at ppfd=" + std::to_string(ppfd);

    phylloptim::leaf_pars<AD> p;
    const phylloptim::leaf_pars<double> seated = l.passive_pars();
    for (std::size_t i = 0; i < p.size(); ++i) p[i] = AD(seated[i]);
    xad::derivative(p[phylloptim::par_PPFD]) = 1.0;
    const double row = xad::derivative(
        l.assim_colimited_kernel<AD>(AD(ci), l.photo_capacity_at<AD>(p)));

    const double h = ppfd * 1e-6;
    Drivers up_d = d, dn_d = d;
    up_d.PPFD = ppfd + h;
    dn_d.PPFD = ppfd - h;
    phylloptim::Leaf up = make_leaf(up_d, {1.0}, {1.0});
    phylloptim::Leaf dn = make_leaf(dn_d, {1.0}, {1.0});
    const double fd = (up.assim_colimited(ci) - dn.assim_colimited(ci)) / (2 * h);

    ok(row != 0.0, "dA/dPPFD is not an exact zero" + at);
    ok(row > 0.0, "more light cannot lower assimilation" + at);
    near(row, fd, 1e-7, "dA/dPPFD against a central difference" + at);
    ok(row < previous, "the light response saturates" + at);
    previous = row;
  }
  ok(compared == 4, "every light level solved");
}

// Every operating point the model answers at hands over FINITE rows.
//
// ⚠️ A NON-FINITE ROW IS NOT REPORTED ANYWHERE DOWNSTREAM, which is why this is
// asserted here. plant records its rows against a constant 1.0, so the row it
// validates is always finite and the sentinel rides in on the tape behind it;
// its census gradient then polls for a DECLARED refusal, finds none, and hands
// back NaN as a number the sweep computed. Measured before this check existed:
// three of five census drivers returned NaN in every trait column of three
// metrics, with nothing refused.
//
// The two routes are separate claims. `duptake_dpsi` returns NaN deliberately
// where its analytic branch is not valid, and a single-layer wet bound -- where
// a shut leaf's collar sits at g_min = 0 -- is such a point; the
// coordinates additionally run the implicit function theorem, which has no
// residual to invert where a leaf with no conductance nets zero at no ci.
void test_every_answered_point_hands_over_finite_rows() {
  printf("every answered operating point hands over finite rows\n");
  using T = odelia::ode::tangent_scalar<double>;
  std::map<std::string, int> reached;
  int draws = 0, uptakes = 0, coords = 0;
  for (int layers : {1, 3}) {
    for (double psi0 : {0.25, 1.0, 2.0, 4.0, 6.0, 17.0}) {
      for (double ppfd : {1.0, 100.0, 900.0}) {
        Drivers d;
        d.PPFD = ppfd;
        std::vector<double> psi_soil, soil_depth;
        for (int i = 0; i < layers; ++i) {
          psi_soil.push_back(psi0 + 0.35 * i);
          soil_depth.push_back(1.0 * (i + 1));
        }
        phylloptim::Leaf l = make_leaf(d, psi_soil, soil_depth);
        l.find_root_collar_psi();
        const auto kind = l.operating_point_kind();
        const std::string name = l.operating_point_kind_name(kind);
        if (kind == phylloptim::Leaf::OperatingPointKind::Unsolved ||
            kind == phylloptim::Leaf::OperatingPointKind::SolverRefused ||
            kind == phylloptim::Leaf::OperatingPointKind::NonFiniteGradient) {
          continue;
        }
        ++reached[name];
        const std::string at = " at " + name + " psi_soil=" +
                               std::to_string(psi0) + " ppfd=" +
                               std::to_string(ppfd) + " layers=" +
                               std::to_string(layers);

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
        ++draws;
        ok(std::isfinite(odelia::util::to_passive(draw.E_up.slope)),
           "the draw's collar slope is finite" + at);
        for (std::size_t j = 0; j < draw.duptake_dp.size(); ++j) {
          ok(std::isfinite(draw.duptake_dp[j]),
             "layer " + std::to_string(j) + "'s supplied collar slope is finite" + at);
        }

        const double curvature =
            kind == phylloptim::Leaf::OperatingPointKind::Interior
                ? l.profit_curvature_at_collar<phylloptim::Leaf::CostCurve::TF24>()
                : std::numeric_limits<double>::quiet_NaN();
        if (kind == phylloptim::Leaf::OperatingPointKind::Interior &&
            !std::isfinite(curvature)) {
          continue;
        }
        const T collar =
            l.operating_collar<phylloptim::Leaf::CostCurve::TF24, T>(draw, pars, curvature);
        const auto out =
            l.billed_outputs_at<phylloptim::Leaf::CostCurve::TF24, T>(collar, draw, pars);
        ++uptakes;
        ok(std::isfinite(odelia::ode::derivative_along(out.profit)),
           "profit's row is finite" + at);
        for (std::size_t j = 0; j < out.uptake.size(); ++j) {
          ok(std::isfinite(odelia::ode::derivative_along(out.uptake[j])),
             "layer " + std::to_string(j) + "'s uptake row is finite" + at);
        }
        // The coordinates are refused by name at a shutdown; billed_outputs_at answers
        // there instead, which the branch above just exercised.
        if (kind != phylloptim::Leaf::OperatingPointKind::HydraulicShutdown) {
          auto cap = l.photo_capacity_at<T>(pars);
          const auto co = l.collar_coords_at<T>(l.opt_psi_leaf_, l.ci_, collar,
                                                draw, pars, cap, true);
          ++coords;
          ok(std::isfinite(odelia::ode::derivative_along(co.sigma.value)) &&
                 std::isfinite(odelia::util::to_passive(co.sigma.slope)),
             "sigma's row and collar slope are finite" + at);
          ok(std::isfinite(odelia::ode::derivative_along(co.ci.value)) &&
                 std::isfinite(odelia::util::to_passive(co.ci.slope)),
             "ci's row and collar slope are finite" + at);
          // Omitting the channel is what profit asks for at an interior point,
          // and it must move no number: the step it multiplies is exactly zero
          // in value, so both coordinates are the residuals' own bit for bit.
          // Asserted here because the value is where an omission would show and
          // the row is not -- a channel that went missing leaves every number
          // finite.
          const auto held = l.collar_coords_at<T>(l.opt_psi_leaf_, l.ci_, collar,
                                                  draw, pars, cap, false);
          ok(odelia::util::to_passive(held.sigma.value) ==
                 odelia::util::to_passive(co.sigma.value),
             "sigma is the same number with the channel omitted" + at);
          ok(odelia::util::to_passive(held.ci.value) ==
                 odelia::util::to_passive(co.ci.value),
             "ci is the same number with the channel omitted" + at);
          ok(std::isfinite(odelia::ode::derivative_along(held.sigma.value)) &&
                 std::isfinite(odelia::ode::derivative_along(held.ci.value)),
             "both coordinates keep a finite row without the channel" + at);
        }
      }
    }
  }
  std::string names;
  for (const auto &[name, n] : reached) {
    names += (names.empty() ? "" : ", ") + name + "(" + std::to_string(n) + ")";
  }
  printf("    %d draws, %d output sets, %d coordinate pairs | kinds: %s\n", draws,
         uptakes, coords, names.c_str());
  // ⚠️ THE CLAIM IS ONLY EVIDENCE ABOUT THE KINDS IT REACHED, and the wet-pinned
  // dark leaf is the one that carried the sentinel: a fixture that stopped
  // reaching it would pass this while checking nothing.
  ok(reached.count("boundary-soil") > 0, "the sweep reached a wet-pinned point");
  ok(reached.count("hydraulic-shutdown") > 0, "the sweep reached a shutdown");
  ok(reached.count("interior") > 0, "the sweep reached an interior point");
}

// TF24_floor's rows at a nonzero price, against re-solving at a perturbed price.
//
// ⚠️ THE TWO ROWS FAIL DIFFERENTLY, which is why both are asserted. The profit row
// is the envelope theorem, so it is right only if the solve's collar is stationary
// in the profit it reports; the collar row is the implicit function theorem on
// dprofit_dcollar_at, so it is right only if dprofit_dcollar_at is the FOC the solve drove to
// zero. The floor's price reads the collar directly, so a solve or a marginal
// missing that term fails one row or the other -- and a difference taken through
// the same solve cannot see it, because it re-solves the same wrong condition.
void test_tf24_floor_rows_at_a_price() {
  printf("TF24_floor's collar and profit rows at lambda_o > 0 match a re-solve\n");
  using T = odelia::ode::tangent_scalar<double>;
  using K = phylloptim::Leaf::CostCurve;
  int compared = 0;
  for (int layers : {1, 3}) {
    for (double price : {1e3, 1e4}) {
      Drivers d;
      d.PPFD = 1500.0;
      std::vector<double> psi_soil, soil_depth;
      for (int i = 0; i < layers; ++i) {
        psi_soil.push_back(1.0 + 0.35 * i);
        soil_depth.push_back(1.0 * (i + 1));
      }
      auto solved = [&](double lam) {
        phylloptim::Leaf l = make_leaf(d, psi_soil, soil_depth);
        l.TF24_floor_lambda_o = lam;
        l.set_model(K::TF24_floor, true);
        l.optimise();
        return l;
      };
      phylloptim::Leaf l = solved(price);
      const std::string at = " at lambda_o=" + std::to_string(price) +
                             " layers=" + std::to_string(layers);
      if (l.operating_point_kind() != phylloptim::Leaf::OperatingPointKind::Interior) {
        continue;
      }

      phylloptim::leaf_pars<T> pars;
      const phylloptim::leaf_pars<double> seated = l.passive_pars();
      for (std::size_t i = 0; i < pars.size(); ++i) pars[i] = T(seated[i]);
      odelia::ode::set_direction(pars[phylloptim::par_TF24_floor_lambda_o], 1.0);
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
      const double curvature = l.profit_curvature_at_collar<K::TF24_floor>();
      const T collar = l.operating_collar<K::TF24_floor, T>(draw, pars, curvature);
      const auto out = l.billed_outputs_at<K::TF24_floor, T>(collar, draw, pars);

      const double h = 1e-4 * price;
      const phylloptim::Leaf up = solved(price + h), dn = solved(price - h);
      const double dcollar = (up.opt_root_psi_ - dn.opt_root_psi_) / (2 * h);
      const double dprofit = (up.profit_ - dn.profit_) / (2 * h);
      // Scaled by the row's own size: both are small numbers per unit price.
      auto rel = [](double got, double want) {
        return std::abs(got - want) / std::max(std::abs(want), 1e-300);
      };
      ok(rel(odelia::ode::derivative_along(collar), dcollar) < 1e-4,
         "the collar row matches a re-solve" + at + ": " +
             std::to_string(odelia::ode::derivative_along(collar)) + " vs " +
             std::to_string(dcollar));
      ok(rel(odelia::ode::derivative_along(out.profit), dprofit) < 1e-4,
         "the profit row matches a re-solve" + at + ": " +
             std::to_string(odelia::ode::derivative_along(out.profit)) + " vs " +
             std::to_string(dprofit));
      ++compared;
    }
  }
  printf("    %d interior points compared\n", compared);
  ok(compared == 4, "every layer count and price reached an interior point");
}

// The supplied rows against a difference of the SOLVE, for every parameter in
// the pack, at every kind operating_collar places by a condition, on every
// curve the surface serves. A difference re-solves, so it re-derives the bound
// or the stationarity condition from scratch and shares nothing with the rows
// but the model -- which is what makes it a referee for rows the surface
// SUPPLIES. The fixtures were chosen on TF24's kinds; another curve may land a
// fixture elsewhere, so each is refereed at the kind it lands on, and coverage
// is asserted per curve: every kind reached, every parameter compared at each.
// A parameter whose step moves the point to another kind is skipped and
// counted, because a difference across a kink is not a derivative.
//
// A parameter the curve's solve does not read on this path (another curve's
// price or level, the single path's resistance, the damage layer, which the
// surface refuses) must give an EXACT zero both as a row and as a difference: a
// surface that starts reading one fails here by name. Every other row is held
// to the tolerance below, and every read parameter must be compared somewhere,
// so a wrong row in any of them fails.
template <phylloptim::Leaf::CostCurve K, bool EB = false>
void referee_every_parameter_at_every_kind(int* compared_out) {
  using T = odelia::ode::tangent_scalar<double>;
  using Kind = phylloptim::Leaf::OperatingPointKind;
  using CC = phylloptim::Leaf::CostCurve;
  namespace pl = phylloptim;
  const std::string curve =
      pl::Leaf::curve_name(static_cast<int>(K)) + (EB ? " with the energy balance" : "");
  // With the energy balance on, as make_pm_leaf seats it: daylight NIR fixed at
  // the fixture's light, so PPFD's row is the light's alone.
  const auto seat_balance = [&](pl::Leaf& leaf, double ppfd) {
    leaf.use_energy_balance_ = EB;
    if (EB) {
      leaf.wind_speed_ = 2.0;
      leaf.d_ = 0.05;
      leaf.NIR_int_ = ppfd / pl::umol_par_per_joule;
      // ProfitMax's thermal cost reads leaf temperature, so it is on there.
      leaf.use_thermal_cost_ = K == CC::ProfitMax;
    }
  };
  struct Fixture {
    std::vector<double> psi_soil;
    double ppfd, vpd;
    Kind want;
    double g_min = 0.0;
  };
  // Golden-grid drivers: layers 0.25 MPa apart; the pinned rows test-gradient.R
  // names, a leaf in the dark (a net source sitting on its wet bound, with a
  // cuticle so the bound is a regular point: at g_min = 0 the solve steps 1e-6
  // of the bracket inside it, and a difference of the solve then measures that
  // step), a soil dry enough that the leaf transpires all its stem can carry (a
  // dry pin at the end of the stem curve), a soil past the end of the curve
  // (shut), and an interior control.
  std::vector<Fixture> fixtures{
      {{1.0, 1.25}, 1500.0, 2.0, Kind::Interior},
      {{4.0, 4.25, 4.5, 4.75, 5.0}, 900.0, 2.0, Kind::BoundarySoil},
      {{3.0, 3.25, 3.5}, 900.0, 3.1, Kind::BoundarySoil},  // 3.1 < e_s(25 C)
      {{6.0}, 500.0, 0.5, Kind::BoundaryDry},
      {{1.0}, 1.0, 2.0, Kind::BoundarySoil, 5e-4},
      {{16.0, 16.25, 16.5}, 900.0, 2.0, Kind::HydraulicShutdown},
      // The same kinds with a cuticle, where g_min is a regular parameter.
      {{1.0, 1.25}, 1500.0, 2.0, Kind::Interior, 5e-4},
      {{6.0}, 500.0, 0.5, Kind::BoundaryDry, 5e-4},
      {{16.0, 16.25, 16.5}, 900.0, 2.0, Kind::HydraulicShutdown, 5e-4},
      // Where the other curves reach the kinds TF24's drivers miss: a dim leaf
      // on a wet soil (CF77's interior), a bright one on a moist soil
      // (LeastCost's and CF77's dry pins), and a dark one on a wet soil (the
      // wet bound of the curves that price the drop or scale A).
      {{0.5}, 100.0, 0.5, Kind::Interior},
      {{0.5}, 100.0, 0.5, Kind::Interior, 5e-4},
      {{0.5}, 500.0, 2.0, Kind::BoundaryDry},
      {{0.3}, 1.0, 0.5, Kind::BoundarySoil, 5e-4},
      // A cuticle large enough that the stem's ceiling caps it on a dry soil,
      // where the residual flux's slopes move with every stem trait.
      {{5.0}, 900.0, 2.0, Kind::Interior, 5e-3},
      {{6.5, 6.75, 7.0}, 900.0, 2.0, Kind::Interior, 5e-3},
      {{7.2}, 900.0, 2.0, Kind::BoundarySoil, 2e-2}};
  // With the energy balance on, the cooling a wide-open leaf earns moves the dry
  // pin: dim light on a dry soil reaches it.
  if (EB) {
    fixtures.push_back({{6.0}, 100.0, 0.5, Kind::BoundaryDry});
    fixtures.push_back({{6.0}, 100.0, 0.5, Kind::BoundaryDry, 5e-4});
  }
  // TF24_floor at a small price, so its fixtures stay on TF24's kinds.
  const double price = 1e3;
  // A slot some curve owns is read only by its owner, and TF24's two by
  // TF24_floor too, whose hydraulic half is TF24's; the single path's
  // resistance and the damage layer, which the surface refuses, by none.
  const auto owns = [](const auto& parameters, int k) {
    for (const auto& q : parameters) if (q.slot == k) return true;
    return false;
  };
  const auto unread = [&](int k) {
    if (k == pl::par_resistance || (k >= pl::par_A && k < pl::par_PPFD)) return true;
    // NIR reaches the solve only through the energy balance.
    if (k == pl::par_NIR) return !EB;
    bool owned = false;
    for (int c = 0; c < pl::Leaf::n_cost_curves; ++c) {
      pl::Leaf::with_curve(static_cast<CC>(c), [&](auto tag) {
        owned = owned || owns(pl::Curve<decltype(tag)::value>::parameters, k);
      });
    }
    if (!owned) return false;
    const bool reads = owns(pl::Curve<K>::parameters, k) ||
                       (K == CC::TF24_floor && owns(pl::Curve<CC::TF24>::parameters, k));
    return !reads;
  };
  const auto par_name = [](int k) {
    return k == pl::par_PPFD  ? std::string("PPFD")
           : k == pl::par_NIR ? std::string("NIR")
                              : pl::gradient::par_names()[std::size_t(k)];
  };
  std::map<std::string, std::vector<int>> per_par_kind;
  double worst_rel = 0.0, worst_abs = 0.0;
  std::string worst_rel_at, worst_abs_at, skipped_at;
  int compared = 0, skipped = 0, zero = 0, at_zero_g_min = 0;
  // ProfitMax's normaliser with the balance on, where |A| can peak inside the
  // stream: points refereed with the peak there and at psi_crit.
  int peak_inside = 0, peak_at_crit = 0;
  for (const Fixture& fx : fixtures) {
    std::vector<double> depth, root;
    for (std::size_t i = 0; i < fx.psi_soil.size(); ++i) {
      depth.push_back(1.0 * double(i + 1));
      root.push_back(1.0 / 0.05 / double(fx.psi_soil.size()));
    }
    const pl::RootNetwork network = fixture::root_network(root, depth);
    // Every slot the pack holds is written to the leaf, the damage layer's
    // included although its gate is off, so a difference sees whatever the
    // solve reads.
    const auto solve = [&](const pl::leaf_pars<double>& p) {
      pl::Leaf l = fixture::multilayer_leaf();
      l.set_traits(p[pl::par_vcmax_25], p[pl::par_stem_c], p[pl::par_stem_P50],
                   p[pl::par_root_c], p[pl::par_root_P50], p[pl::par_TF24_beta2],
                   p[pl::par_jmax_25], p[pl::par_phi_J],
                   p[pl::par_curv_fact_elec_trans], p[pl::par_curv_fact_colim],
                   p[pl::par_TF24_cost_scale], p[pl::par_R_d_25],
                   p[pl::par_JS22_gamma], p[pl::par_CMax_a], p[pl::par_CMax_b],
                   p[pl::par_g_min]);
      l.CF77_lambda_ = p[pl::par_CF77_lambda];
      l.TF24_floor_lambda_o = p[pl::par_TF24_floor_lambda_o];
      l.LeastCost_beta = p[pl::par_LeastCost_beta];
      l.A_ = p[pl::par_A];
      l.I_ = p[pl::par_I];
      l.D_ = p[pl::par_D];
      l.topt_offset_ = p[pl::par_topt_offset];
      l.dTopt_max_ = p[pl::par_dTopt_max];
      l.K_A_ = p[pl::par_K_A];
      l.k_i_ = p[pl::par_k_i];
      l.k_i_baseline_ = p[pl::par_k_i_baseline];
      l.k_rec_ = p[pl::par_k_rec];
      l.k_rec_Ea_ = p[pl::par_k_rec_Ea];
      l.c_acclim_maint_ = p[pl::par_c_acclim_maint];
      l.c_repair_maint_ = p[pl::par_c_repair_maint];
      l.c_protect_maint_ = p[pl::par_c_protect_maint];
      l.c_repair_flux_ = p[pl::par_c_repair_flux];
      seat_balance(l, fx.ppfd);
      if (EB) l.NIR_int_ = p[pl::par_NIR];
      l.set_physiology(network, p[pl::par_PPFD], fx.psi_soil, depth,
                       p[pl::par_kmax], fx.vpd, 40.0, 25.0, 21.0, 101.3);
      l.set_model(K, true);
      l.optimise();
      return l;
    };
    pl::Leaf base = fixture::multilayer_leaf();
    seat_balance(base, fx.ppfd);
    base.set_physiology(network, fx.ppfd, fx.psi_soil, depth, 0.000157 / 5.0,
                        fx.vpd, 40.0, 25.0, 21.0, 101.3);
    pl::leaf_pars<double> p0 = base.passive_pars();
    p0[pl::par_g_min] = fx.g_min;
    p0[pl::par_TF24_floor_lambda_o] = price;
    // The single path's slot, NA on this one; any value, since nothing reads it.
    p0[pl::par_resistance] = 1e3;
    // NIR is unset off the balance, where nothing reads it; any finite value.
    if (!EB) p0[pl::par_NIR] = 0.0;
    pl::Leaf l;
    try {
      l = solve(p0);
    } catch (const pl::util::infeasible_error&) {
      skipped_at += " fixture:infeasible";
      continue;
    }
    const Kind kind = l.operating_point_kind();
    const std::string name = l.operating_point_kind_name(kind);
    if (kind != Kind::Interior && kind != Kind::BoundarySoil &&
        kind != Kind::BoundaryDry && kind != Kind::HydraulicShutdown) {
      skipped_at += " fixture:" + name;
      continue;
    }
    if (!fixture::surface_has_objective(l)) {
      skipped_at += " fixture:no-objective";
      continue;
    }
    const double curvature = kind == Kind::Interior
                                 ? l.profit_curvature_at_collar<K>()
                                 : std::numeric_limits<double>::quiet_NaN();
    // operating_collar refuses an interior point too flat to divide by.
    if (kind == Kind::Interior && !(curvature < -pl::Leaf::interior_curvature_floor)) {
      skipped_at += " fixture:flat-interior";
      continue;
    }
    // The solve places psi_leaf by inverting the stem's table, whose slope in E
    // departs from the closed form's 1/(kmax f) along the curve's tail (1e-11 at
    // 3 MPa, 1e-5 at 8, 5e-4 past 10): a difference of psi_leaf carries that and
    // the row does not. Measured at the point, it widens psi_leaf's tolerance only.
    // Both measured on a copy: the inverses move the leaf, and the forward
    // difference reads the base leaf itself.
    pl::Leaf probe = l;
    double inverse_err = 0.0;
    try {
      const double E = probe.transpiration(probe.opt_psi_leaf_, probe.opt_root_psi_);
      const double dE = 1e-6 * E;
      const double fd_inv =
          (probe.transpiration_to_psi_leaf(E + dE, probe.opt_root_psi_) -
           probe.transpiration_to_psi_leaf(E - dE, probe.opt_root_psi_)) / (2.0 * dE);
      const double closed = 1.0 / (probe.leaf_specific_conductance_max_ *
                                   probe.stem_curve_integral_deriv(probe.opt_psi_leaf_));
      inverse_err = std::abs(fd_inv - closed) / std::abs(closed);
      // A shut leaf's stem sits at its collar, where there is no slope to compare.
      if (!std::isfinite(inverse_err)) inverse_err = 0.0;
    } catch (const std::exception&) {
      // At the end of the stem curve there is no inverse to step across.
    }
    // psi_leaf is the collar's function, so the solve's noise in it is the
    // collar's times dpsi_leaf/dcollar, which on the curve's tail runs to 1e5.
    double psi_per_collar = 0.0;
    try {
      const double dp = 1e-6;
      psi_per_collar = std::abs(
          (probe.find_psi_leaf_from_psi_root(probe.opt_root_psi_ + dp, probe.supply_psi_soil()) -
           probe.find_psi_leaf_from_psi_root(probe.opt_root_psi_ - dp, probe.supply_psi_soil())) /
          (2.0 * dp));
      if (!std::isfinite(psi_per_collar)) psi_per_collar = 0.0;
    } catch (const std::exception&) {
    }
    if (K == CC::ProfitMax && EB) {
      (l.profitmax_peak_psi_ == l.psi_crit ? peak_at_crit : peak_inside) += 1;
    }
    for (int k = 0; k < pl::n_leaf_pars; ++k) {
      // At g_min = 0 the solve steps 1e-6 of the bracket inside the wet bound,
      // and at any g_min > 0 it sits on it (hazard 18), so a forward difference
      // carries that step over h (-4.3e-6 MPa / h in the collar here). The row
      // is the slope from above; on a bound it is refereed at the fixtures with
      // a cuticle. An interior point has no such step.
      if (k == pl::par_g_min && fx.g_min == 0.0 && kind != Kind::Interior) {
        ++at_zero_g_min;
        continue;
      }
      pl::leaf_pars<T> tp;
      for (std::size_t i = 0; i < tp.size(); ++i) tp[i] = T(p0[i]);
      odelia::ode::set_direction(tp[std::size_t(k)], 1.0);
      std::vector<T> soil, r_h, r_v;
      for (int i = 0; i < l.supply_n_layers(); ++i) {
        soil.push_back(T(l.roots_.psi_soil_[std::size_t(i)]));
        r_h.push_back(T(l.roots_.network_.r_R_H_min[std::size_t(i)]));
        r_v.push_back(T(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
      }
      const pl::SupplyInputs<T> supply{soil, r_h, r_v, tp[pl::par_root_P50],
                                       tp[pl::par_root_c]};
      const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), supply);
      const T collar = l.operating_collar<K, T>(draw, tp, curvature);
      const auto out = l.billed_outputs_at<K, T>(collar, draw, tp);
      const auto phys = l.physical_outputs_at<K, T>(collar, draw, tp);
      const auto row_of = [](const T& x) { return odelia::ode::derivative_along(x); };
      // The two variants share their core: what both return is the same number.
      bool shared = row_of(phys.profit) == row_of(out.profit) &&
                    row_of(phys.residual_flux) == row_of(out.residual_flux);
      for (std::size_t j = 0; j < out.uptake.size(); ++j) {
        shared = shared && row_of(phys.uptake[j]) == row_of(out.uptake[j]);
      }
      ok(shared, curve + ": the billed and physical variants share profit, "
                 "uptake and the residual flux at " + name + "/" + par_name(k));
      std::vector<double> rows{row_of(collar),       row_of(out.profit),
                               row_of(out.carbon),   row_of(phys.carbon),
                               row_of(phys.A),       row_of(phys.E),
                               row_of(phys.psi_leaf), row_of(phys.ci),
                               row_of(phys.gc),      row_of(out.residual_flux)};
      for (const T& u : out.uptake) rows.push_back(row_of(u));

      // RELATIVE 1e-4, read off a step sweep rather than assumed: every row's
      // error-versus-step V has its floor at or above 1e-5, and at the
      // gradient's own 1e-6 kmax's profit row at a wet bound reads the
      // difference's noise (2e-4 against a floor of 9e-7). At 1e-4 every row
      // larger than 1e-6 agrees to 1e-5. A parameter at zero takes a forward
      // step, since some (CMax_b among them) may not go negative.
      const double v = p0[std::size_t(k)];
      const bool forward = v == 0.0;
      const double h = forward ? 1e-7 : std::abs(v) * 1e-4;
      pl::leaf_pars<double> p_up = p0, p_down = p0;
      p_up[std::size_t(k)] += h;
      if (!forward) p_down[std::size_t(k)] -= h;
      pl::Leaf up, down;
      try {
        up = solve(p_up);
        down = forward ? l : solve(p_down);
      } catch (const pl::util::infeasible_error&) {
        ++skipped;
        skipped_at += " " + name + "/" + par_name(k) + "(infeasible)";
        continue;
      }
      const double span = forward ? h : 2.0 * h;
      if (up.operating_point_kind() != kind ||
          down.operating_point_kind() != kind) {
        ++skipped;
        skipped_at += " " + name + "/" + par_name(k);
        continue;
      }
      // Each output's difference, and its size at the two ends, which sets the
      // noise floor below.
      std::vector<double> fd, mag;
      const auto diff = [&](double a, double b) {
        fd.push_back((a - b) / span);
        mag.push_back(std::max(std::abs(a), std::abs(b)));
      };
      diff(up.opt_root_psi_, down.opt_root_psi_);
      diff(up.profit_, down.profit_);
      diff(up.carbon(), down.carbon());
      diff(up.carbon(), down.carbon());
      diff(up.assim_colimited_, down.assim_colimited_);
      diff(up.transpiration_, down.transpiration_);
      diff(up.opt_psi_leaf_, down.opt_psi_leaf_);
      diff(up.ci_, down.ci_);
      diff(up.stom_cond_CO2_, down.stom_cond_CO2_);
      diff(up.residual_flux_, down.residual_flux_);
      for (std::size_t j = 0; j < out.uptake.size(); ++j) {
        diff(up.soil_consumption_[j], down.soil_consumption_[j]);
      }
      static const char* const output_name[] = {
          "collar", "profit", "carbon (billed)", "carbon (physical)", "A", "E",
          "psi_leaf", "ci", "gc", "residual_flux"};
      constexpr std::size_t n_named = std::size(output_name);
      for (std::size_t j = 0; j < rows.size(); ++j) {
        const std::string at =
            " (" + curve + ") at " + name +
            (fx.g_min > 0.0 ? " with a cuticle" : "") + " d" +
            (j < n_named ? std::string(output_name[j])
                         : "uptake[" + std::to_string(j - n_named) + "]") +
            "/d" + par_name(k);
        if (unread(k)) {
          // The row is exactly zero, and so is the difference: the solve reads
          // nothing of the parameter, not even through rounding.
          ok(rows[j] == 0.0 && fd[j] == 0.0,
             "a parameter the solve does not read has a zero row and a zero "
             "difference" + at + ": row " + pl::util::format_double(rows[j]) +
                 ", difference " + pl::util::format_double(fd[j]));
          ++zero;
          continue;
        }
        const double err = std::abs(rows[j] - fd[j]);
        const double scale = std::max(std::abs(rows[j]), std::abs(fd[j]));
        // The absolute floor is on the row times the parameter's own scale,
        // so a parameter of order 1e3 (PPFD, the price) is not let off rows of
        // order 1e-9 per unit. The second floor is the solve's own noise in the
        // output over the step (hazard 15): 1e-10 of its size, ten times the
        // worst measured, ci in the dark, where A'(ci) is near zero and a step
        // in kmax or stem_P50 moved ci by 1e-11 of itself.
        const double unit = std::max(1.0, std::abs(v));
        const double noise =
            1e-10 * (j == 6 ? std::max(mag[j], psi_per_collar * mag[0]) : mag[j]) / h;
        const double floor = 1e-7 / unit + noise;
        if (err > worst_abs) { worst_abs = err; worst_abs_at = at; }
        // Worst relative among the rows the relative term governs.
        if (1e-4 * scale > floor && err / scale > worst_rel) {
          worst_rel = err / scale;
          worst_rel_at = at;
        }
        ++compared;
        auto& seen = per_par_kind[name];
        seen.resize(pl::n_leaf_pars, 0);
        ++seen[std::size_t(k)];
        const double rel_tol = 1e-4 + (j == 6 ? 2.0 * inverse_err : 0.0);
        ok(err <= rel_tol * scale + floor,
           "the supplied row matches a difference of the solve" + at +
               ": row " + pl::util::format_double(rows[j]) + ", difference " +
               pl::util::format_double(fd[j]));
      }
    }
  }
  std::string kinds;
  for (const auto& [kind, seen] : per_par_kind) {
    int n = 0;
    for (int c : seen) n += c;
    kinds += " " + kind + " " + std::to_string(n);
  }
  printf("    %s: %d rows compared (%s ), %d unread rows exactly zero,\n"
         "    g_min passed over at %d fixtures without a cuticle,\n"
         "    %d parameters skipped across a kind change:%s\n"
         "    worst relative %.3g%s\n    worst absolute %.3g%s\n",
         curve.c_str(), compared, kinds.c_str(), zero, at_zero_g_min, skipped,
         skipped_at.c_str(), worst_rel, worst_rel_at.c_str(), worst_abs,
         worst_abs_at.c_str());
  if (K == CC::ProfitMax && EB) {
    printf("    ProfitMax's |A| peak: %d points inside the stream, %d at psi_crit\n",
           peak_inside, peak_at_crit);
    ok(peak_inside > 0, curve + ": rows were compared where |A| peaks inside the stream");
    ok(peak_at_crit > 0, curve + ": and where it peaks at psi_crit");
  }
  // A referee that skipped a parameter at a kind would pass on what is left.
  for (Kind want : fixture::surface_kinds(K)) {
    ok(per_par_kind.count(pl::Leaf::operating_point_kind_name(want)) == 1,
       curve + ": rows were compared at " + pl::Leaf::operating_point_kind_name(want));
  }
  for (const auto& [kind, seen] : per_par_kind) {
    for (int k = 0; k < pl::n_leaf_pars; ++k) {
      if (!unread(k)) {
        ok(seen[std::size_t(k)] > 0,
           curve + ": " + par_name(k) + " was compared at " + kind);
      }
    }
  }
  ok(skipped <= 3, curve + ": at most three parameters were skipped across a "
                   "kind change");
  *compared_out += compared;
}

void test_supplied_rows_match_a_resolve_at_every_kind() {
  printf("supplied rows against a difference of the solve, every parameter, "
         "at every kind\n");
  using CC = phylloptim::Leaf::CostCurve;
  int compared = 0;
  referee_every_parameter_at_every_kind<CC::TF24>(&compared);
  referee_every_parameter_at_every_kind<CC::TF24_floor>(&compared);
  referee_every_parameter_at_every_kind<CC::CF77>(&compared);
  referee_every_parameter_at_every_kind<CC::JS22>(&compared);
  referee_every_parameter_at_every_kind<CC::CMax>(&compared);
  referee_every_parameter_at_every_kind<CC::SOX>(&compared);
  referee_every_parameter_at_every_kind<CC::JW26>(&compared);
  referee_every_parameter_at_every_kind<CC::LeastCost>(&compared);
  referee_every_parameter_at_every_kind<CC::ProfitMax>(&compared);
  referee_every_parameter_at_every_kind<CC::TF24, true>(&compared);
  referee_every_parameter_at_every_kind<CC::TF24_floor, true>(&compared);
  referee_every_parameter_at_every_kind<CC::CF77, true>(&compared);
  referee_every_parameter_at_every_kind<CC::JS22, true>(&compared);
  referee_every_parameter_at_every_kind<CC::CMax, true>(&compared);
  referee_every_parameter_at_every_kind<CC::SOX, true>(&compared);
  referee_every_parameter_at_every_kind<CC::JW26, true>(&compared);
  referee_every_parameter_at_every_kind<CC::LeastCost, true>(&compared);
  referee_every_parameter_at_every_kind<CC::ProfitMax, true>(&compared);
  ok(compared >= 3300, "the referee compared at least 3300 rows");
}

// A shut leaf whose roots outlast its stem, so set_shut_state seats the stem at
// the collar. The surface used to place the stem by inverting T1 against the
// draw there, which divided the draw's root rows by the stem's conductivity at
// the collar (~1e-20 past the end of the stem curve, ~1e-18 just short of it)
// and billed the root traits up to 9e-3 of profit. The leaf is shut and stays
// shut, so a difference of the solve gives 0 and those rows must be zero to
// rounding; every other row is held to the referee's tolerance against the
// same difference.
//
// Only fixtures whose seating is decided exactly: past the stem's end the
// residual flux is exactly zero, and at g_min = 0 it always is. Short of the
// end at g_min > 0 the seating follows the sign of an uptake of ~1e-22, which
// is rounding and need not agree across platforms.
void test_shut_leaf_seated_at_its_collar_bills_the_roots_nothing() {
  printf("a shut leaf seated at its collar bills the root traits nothing\n");
  using T = odelia::ode::tangent_scalar<double>;
  using Kind = phylloptim::Leaf::OperatingPointKind;
  using CC = phylloptim::Leaf::CostCurve;
  namespace pl = phylloptim;
  struct Fixture {
    double psi_top, g_min;
    bool past_the_end;
  };
  const std::vector<Fixture> fixtures{
      {16.0, 5e-4, true}, {16.0, 0.0, true}, {15.1, 0.0, false}};
  const std::vector<int> pars{
      pl::par_vcmax_25, pl::par_stem_c, pl::par_stem_P50, pl::par_root_c,
      pl::par_root_P50, pl::par_TF24_beta2, pl::par_jmax_25, pl::par_phi_J,
      pl::par_curv_fact_elec_trans, pl::par_curv_fact_colim,
      pl::par_TF24_cost_scale, pl::par_R_d_25, pl::par_kmax, pl::par_g_min};
  int compared = 0, skipped = 0;
  double worst_root_profit = 0.0;
  for (const Fixture& fx : fixtures) {
    const std::vector<double> psi_soil{fx.psi_top, fx.psi_top + 0.25,
                                       fx.psi_top + 0.5};
    std::vector<double> depth, root;
    for (std::size_t i = 0; i < psi_soil.size(); ++i) {
      depth.push_back(1.0 * double(i + 1));
      root.push_back(1.0 / 0.05 / double(psi_soil.size()));
    }
    const pl::RootNetwork network = fixture::root_network(root, depth);
    const auto solve = [&](const pl::leaf_pars<double>& p) {
      pl::Leaf l = fixture::multilayer_leaf();
      l.set_traits(p[pl::par_vcmax_25], p[pl::par_stem_c], p[pl::par_stem_P50],
                   p[pl::par_root_c], p[pl::par_root_P50], p[pl::par_TF24_beta2],
                   p[pl::par_jmax_25], p[pl::par_phi_J],
                   p[pl::par_curv_fact_elec_trans], p[pl::par_curv_fact_colim],
                   p[pl::par_TF24_cost_scale], p[pl::par_R_d_25],
                   p[pl::par_JS22_gamma], p[pl::par_CMax_a], p[pl::par_CMax_b],
                   p[pl::par_g_min]);
      l.set_physiology(network, 900.0, psi_soil, depth, p[pl::par_kmax], 2.0,
                       40.0, 25.0, 21.0, 101.3);
      l.find_root_collar_psi();
      return l;
    };
    pl::Leaf base = fixture::multilayer_leaf();
    base.set_physiology(network, 900.0, psi_soil, depth, 0.000157 / 5.0, 2.0,
                        40.0, 25.0, 21.0, 101.3);
    pl::leaf_pars<double> p0 = base.passive_pars();
    p0[pl::par_g_min] = fx.g_min;
    p0[pl::par_root_P50] = 6.0;
    const pl::Leaf l = solve(p0);
    const std::string at = " at soil " + pl::util::format_double(fx.psi_top) +
                           ", g_min " + pl::util::format_double(fx.g_min);
    ok(l.operating_point_kind() == Kind::HydraulicShutdown, "the leaf is shut" + at);
    ok(l.opt_psi_leaf_ == l.opt_root_psi_, "the stem is seated at the collar" + at);
    ok((l.opt_root_psi_ > pl::vulnerability_psi_end(l.stem_b, l.stem_c)) ==
           fx.past_the_end,
       std::string("the collar is ") + (fx.past_the_end ? "past" : "short of") +
           " the end of the stem curve" + at);
    if (l.operating_point_kind() != Kind::HydraulicShutdown) continue;
    for (int k : pars) {
      // A difference across g_min = 0 is not a derivative.
      if (k == pl::par_g_min && fx.g_min == 0.0) continue;
      pl::Leaf lk = l;
      pl::leaf_pars<T> tp;
      for (std::size_t i = 0; i < tp.size(); ++i) tp[i] = T(p0[i]);
      odelia::ode::set_direction(tp[std::size_t(k)], 1.0);
      std::vector<T> soil, r_h, r_v;
      for (int i = 0; i < lk.supply_n_layers(); ++i) {
        soil.push_back(T(lk.roots_.psi_soil_[std::size_t(i)]));
        r_h.push_back(T(lk.roots_.network_.r_R_H_min[std::size_t(i)]));
        r_v.push_back(T(lk.roots_.network_.r_R_V_sum[std::size_t(i)]));
      }
      const pl::SupplyInputs<T> supply{soil, r_h, r_v, tp[pl::par_root_P50],
                                       tp[pl::par_root_c]};
      const auto draw = lk.supply_draw_at<T>(T(lk.opt_root_psi_), supply);
      const T collar = lk.operating_collar<CC::TF24, T>(
          draw, tp, std::numeric_limits<double>::quiet_NaN());
      const auto out = lk.billed_outputs_at<CC::TF24, T>(collar, draw, tp);
      const double h = std::abs(p0[std::size_t(k)]) * 1e-4;
      pl::leaf_pars<double> p_up = p0, p_down = p0;
      p_up[std::size_t(k)] += h;
      p_down[std::size_t(k)] -= h;
      const pl::Leaf up = solve(p_up), down = solve(p_down);
      const std::string name = pl::gradient::par_names()[std::size_t(k)];
      if (up.operating_point_kind() != Kind::HydraulicShutdown ||
          down.operating_point_kind() != Kind::HydraulicShutdown) {
        ++skipped;
        continue;
      }
      const double rows[2] = {odelia::ode::derivative_along(collar),
                              odelia::ode::derivative_along(out.profit)};
      const double fd[2] = {(up.opt_root_psi_ - down.opt_root_psi_) / (2.0 * h),
                            (up.profit_ - down.profit_) / (2.0 * h)};
      for (int j = 0; j < 2; ++j) {
        const std::string what =
            at + " d" + (j == 0 ? "collar" : "profit") + "/d" + name;
        const double scale = std::max(std::abs(rows[j]), std::abs(fd[j]));
        ++compared;
        ok(std::abs(rows[j] - fd[j]) <= 1e-4 * scale + 1e-7,
           "the row matches a difference of the solve" + what + ": row " +
               pl::util::format_double(rows[j]) + ", difference " +
               pl::util::format_double(fd[j]));
      }
      if (k == pl::par_root_P50 || k == pl::par_root_c) {
        worst_root_profit = std::max(worst_root_profit, std::abs(rows[1]));
        ok(std::abs(fd[1]) <= 1e-12,
           "a difference of the solve leaves profit where it is" + at + " d/d" + name);
        ok(std::abs(rows[1]) <= 1e-12,
           "the profit row is zero to rounding" + at + " d/d" + name + ": " +
               pl::util::format_double(rows[1]));
      }
    }
  }
  printf("    %d rows compared, %d parameters skipped across a kind change, "
         "worst root |dprofit| %.3g\n",
         compared, skipped, worst_root_profit);
  ok(compared >= 70, "the test compared at least 70 rows");
}

// Three of the surface's preconditions fail with every number finite, so each is
// checked by name: the curve the solve was seated on, the supply path the draw
// walks, and the draw's E_up at a shutdown, which supply_draw_at returns before
// writing. The control for each is the same call on the matching leaf.
void test_surface_refuses_a_curve_or_supply_it_was_not_solved_on() {
  printf("the derivative surface refuses another curve or supply path\n");
  using CC = phylloptim::Leaf::CostCurve;
  using Kind = phylloptim::Leaf::OperatingPointKind;
  Drivers d;
  d.PPFD = 1500.0;
  phylloptim::Leaf l = make_leaf(d, {1.0, 1.5}, {1.0, 2.0});
  l.find_root_collar_psi();
  ok(l.operating_point_kind() == Kind::Interior, "the fixture is an interior point");
  const phylloptim::leaf_pars<double> pars = l.passive_pars();
  const auto held_supply = l.roots_.held_supply();
  const auto draw = l.supply_draw_at<double>(l.opt_root_psi_, held_supply);
  ok(l.operating_collar<CC::TF24, double>(draw, pars) == l.opt_root_psi_,
     "control: operating_collar answers on the seated curve");
  ok(std::isfinite(l.billed_outputs_at<CC::TF24, double>(l.opt_root_psi_, draw, pars).profit),
     "control: billed_outputs_at answers on the seated curve");
  const std::string wrong =
      "asked for the TF24_floor curve, and the leaf is seated on TF24";
  ok(throws_with([&] { l.operating_collar<CC::TF24_floor, double>(draw, pars); },
                 "Leaf::operating_collar: " + wrong),
     "operating_collar refuses a curve the leaf is not seated on");
  ok(throws_with([&] {
       l.billed_outputs_at<CC::TF24_floor, double>(l.opt_root_psi_, draw, pars);
     }, "Leaf::billed_outputs_at: " + wrong),
     "billed_outputs_at refuses a curve the leaf is not seated on");

  phylloptim::Leaf single = fixture::multilayer_leaf();
  single.setup_transpiration(100);
  single.setup_root_vulnerability(100);
  single.set_supply_single();
  std::vector<double> psi_soil{1.0}, depth{1.0};
  single.set_physiology(fixture::series_resistance(1.0e3), d.PPFD, psi_soil, depth,
                        d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                        d.atm_o2_kpa, d.atm_kpa);
  single.find_root_collar_psi();
  ok(std::isfinite(single.profit_), "control: the single-potential leaf solves");
  const auto single_supply = single.roots_.held_supply();
  ok(throws_with([&] {
       single.supply_draw_at<double>(single.opt_root_psi_, single_supply);
     }, "Leaf::supply_draw_at: the derivative surface reads the multi-layer "
        "supply, and this leaf is on the single path"),
     "supply_draw_at refuses the single-potential path");

  phylloptim::Leaf shut = make_leaf(Drivers{}, {17.0}, {1.0});
  shut.find_root_collar_psi();
  ok(shut.operating_point_kind() == Kind::HydraulicShutdown,
     "the dry fixture is a hydraulic shutdown");
  const auto shut_supply = shut.roots_.held_supply();
  const auto dry = shut.supply_draw_at<double>(shut.opt_root_psi_, shut_supply);
  ok(dry.E_up.value == 0.0 && dry.E_up.slope == 0.0,
     "a shutdown draw's E_up is zero, not indeterminate");
}

// A consumer that records each solve and replays it later (plant does, on every
// step of a reverse sweep) must be handed the rows the solve itself would hand
// over. A pinned solve places its collar at a POLISHED bound and corrects its
// rows back to that bound's root, and neither is in the replay's inputs; a
// replay that lost them gave root-trait rows of the wrong sign at one-layer wet
// bounds. Asserted bit for bit, across interior and both pinned kinds.
void test_a_replay_hands_over_the_rows_its_solve_did() {
  printf("a replayed operating point hands over its solve's rows\n");
  std::map<std::string, int> seen;
  int compared = 0, differ = 0;
  for (int layers : {1, 2}) {
    for (double psi0 : {2.5, 4.0, 5.0, 5.5}) {
      for (double ppfd : {100.0, 900.0}) {
        for (double temp : {25.0, 40.0}) {
          Drivers d;
          d.PPFD = ppfd;
          d.leaf_temp = temp;
          std::vector<double> psi_soil, soil_depth;
          for (int i = 0; i < layers; ++i) {
            psi_soil.push_back(psi0 + 0.35 * i);
            soil_depth.push_back(1.0 * (i + 1));
          }
          phylloptim::Leaf fresh = make_leaf(d, psi_soil, soil_depth);
          fresh.find_root_collar_psi();
          const auto kind = fresh.operating_point_kind();
          if (kind != phylloptim::Leaf::OperatingPointKind::Interior &&
              kind != phylloptim::Leaf::OperatingPointKind::BoundarySoil &&
              kind != phylloptim::Leaf::OperatingPointKind::BoundaryDry) {
            continue;
          }
          ++seen[phylloptim::Leaf::operating_point_kind_name(kind)];
          const double collar = fresh.opt_root_psi_;
          phylloptim::Leaf replay = make_leaf(d, psi_soil, soil_depth);
          replay.replay_operating_point(collar, kind);
          if (replay.opt_root_psi_ != collar) ++differ;
          for (int par : {phylloptim::par_vcmax_25, phylloptim::par_stem_P50,
                          phylloptim::par_root_c, phylloptim::par_root_P50}) {
            ++compared;
            if (profit_row_at_current_point(replay, par) !=
                profit_row_at_current_point(fresh, par)) {
              ++differ;
            }
          }
        }
      }
    }
  }
  printf("       (%d rows; interior %d, boundary-soil %d, boundary-dry %d points)\n",
         compared, seen["interior"], seen["boundary-soil"], seen["boundary-dry"]);
  ok(seen["interior"] > 0 && seen["boundary-soil"] > 0,
     "the grid reaches an interior point and a wet bound (test is not vacuous)");
  ok(differ == 0, "every replayed collar and row equals its solve's, bit for bit");
}

// The dry pin and the interior answer meet continuously. The collar solve
// decides a dry pin on the marginal AT the end of the collar interval, solved to
// collar precision; decided a step inside an end placed only to the bracket's
// tolerance, a stationary point between that step and the end was reported at
// the end, past its maximum, and the solve jumped by 2.4e-4 MPa in the collar
// and 6e-7 in profit where the decision flipped. A stand of such leaves carries
// thousands of those flips, and a difference of a whole run integrated them as
// a bias that no step size removed, 0.37 against 0.22 on plant's seasonal stand.
//
// The first fixture is a leaf of that stand: dry topsoil holding most of the
// roots over a wet subsoil, where the collar objective's maximum sits within
// micro-MPa of the dry end for a whole range of the top layer's potential. The
// top layer is swept finely through the point the old decision flipped at, and
// every step of profit and of the collar is held to the steps beside it. The
// sweep is asserted to run along the end: a step of 1e-5 MPa past the collar
// is outside the interval on most of its points. The second fixture is a true
// pin, where the marginal is not falling at the collar it reports.
void test_a_dry_pin_meets_the_interior_answer_continuously() {
  printf("the dry pin and the interior answer meet continuously\n");
  using Kind = phylloptim::Leaf::OperatingPointKind;
  using CC = phylloptim::Leaf::CostCurve;
  namespace pl = phylloptim;
  // One of plant's seasonal-stand leaves, as it solved it: TF24_floor at
  // lambda_o = 0, five soil layers 0.3 m thick, roots in all five.
  const std::vector<double> soil0{4.2195504249414908, 0.092694683070275272,
                                  0.060412939335358067, 0.048956447693756115,
                                  0.042802900662827763};
  const std::vector<double> depth{0.3, 0.6, 0.9, 1.2, 1.5};
  pl::RootNetwork net;
  net.r_R_H_min = {189.35400538213139, 3668.3966498423301, 9406.0554611473381,
                   23250.799912514827, 91864.343073921438};
  net.r_R_V_sum = {942.31463854872447, 19197.982672469967, 66006.941026179644,
                   181713.86294375348, 638874.06435868028};
  const double kmax = 0.00017675922327204239, ppfd = 676.27236337813406;
  auto solve = [&](double top) {
    pl::Leaf l = fixture::multilayer_leaf();
    l.set_model(CC::TF24_floor, true);
    l.TF24_floor_lambda_o = 0.0;
    // plant's TF24 traits, in set_traits' order.
    l.set_traits(96.0, 2.04, 2.8887256653360196, 2.680147, 3.4, 1.5, 157.44,
                 0.35294117647058826, 0.7, 0.99, 7.5, 1.44, 0.4717887793,
                 0.6385596982, 0.0, 0.0);
    std::vector<double> soil = soil0;
    soil[0] = top;
    l.set_physiology(net, ppfd, soil, depth, kmax, 1.0, 40.0, 25.0, 21.0, 101.3);
    l.find_root_collar_psi_for<CC::TF24_floor>();
    return l;
  };
  // Through the point the old decision flipped at (4.21997), finely enough
  // that its jump of 5e-6 MPa in the collar stands against the sweep's own
  // steps of 1e-7: 1001 solves, a few milliseconds.
  const double lo = 4.21987, hi = 4.22007, step = 2e-7;
  std::vector<double> profit, collar;
  int along_the_end = 0, n = 0, pinned = 0;
  for (double top = lo; top <= hi + 0.5 * step; top += step) {
    pl::Leaf l = solve(top);
    ok(l.operating_point_kind() == Kind::Interior ||
           l.operating_point_kind() == Kind::BoundaryDry,
       "the fixture solves to an interior point or a dry pin");
    pinned += l.operating_point_kind() == Kind::BoundaryDry;
    profit.push_back(l.profit_);
    collar.push_back(l.opt_root_psi_);
    bool read = false;
    l.dprofit_at_collar_psi<CC::TF24_floor>(l.opt_root_psi_ + 1e-5, &read);
    along_the_end += !read;
    ++n;
  }
  // The kinds are printed, not asserted: which side of the decision these
  // points fall on is the solve's to say, and the test is about the join.
  printf("    %d points, %d within 1e-5 MPa of the dry end, %d pinned\n", n,
         along_the_end, pinned);
  ok(along_the_end > n / 2, "the sweep runs along the dry end (test is not vacuous)");
  // Each step against the next: a jump is a step out of all proportion to the
  // one beside it. The old decision's jump was 5e-6 MPa in the collar against
  // steps of ~1e-7 here.
  double worst_profit = 0.0, worst_collar = 0.0;
  for (std::size_t i = 1; i + 1 < profit.size(); ++i) {
    const double dp = std::abs(profit[i] - profit[i - 1]);
    const double dp_next = std::abs(profit[i + 1] - profit[i]);
    const double dc = std::abs(collar[i] - collar[i - 1]);
    const double dc_next = std::abs(collar[i + 1] - collar[i]);
    worst_profit = std::max(worst_profit, dp / std::max(dp_next, 1e-300));
    worst_collar = std::max(worst_collar, dc / std::max(dc_next, 1e-300));
  }
  printf("    largest step against the next: profit %.3g, collar %.3g\n",
         worst_profit, worst_collar);
  ok(worst_profit < 3.0, "profit is continuous along the dry end");
  ok(worst_collar < 3.0, "the collar is continuous along the dry end");

  // A true pin: one soil layer at 6 MPa, where profit still rises at the end.
  {
    pl::Leaf l = make_leaf([] { Drivers d; d.PPFD = 500.0; d.atm_vpd = 0.5; return d; }(),
                           {6.0}, {1.0});
    l.find_root_collar_psi();
    ok(l.operating_point_kind() == Kind::BoundaryDry, "the dry fixture pins");
    // A pin past its maximum, the old decision's failure, has a falling
    // marginal of order one at the collar it reports, and a step inside it
    // EARNS by that marginal times the step. The step is 1e-7 MPa, so such a
    // pin earns ~1e-7, well clear of what the solve resolves; a true pin loses
    // the rise over that step. The marginal AT the end is not asserted: the
    // profit's maximum can sit a hair inside the end, where the marginal has
    // already turned (it falls from zero to order one across the last
    // nano-MPa), and the pin is still the answer the collar resolves.
    pl::Leaf inside = make_leaf([] { Drivers d; d.PPFD = 500.0; d.atm_vpd = 0.5; return d; }(),
                                {6.0}, {1.0});
    const double p_in = inside.evaluate_root_collar_psi_for<CC::TF24>(l.opt_root_psi_ - 1e-7);
    ok(p_in < l.profit_,
       "a step inside a dry pin loses: " + std::to_string(p_in) + " against " +
           std::to_string(l.profit_));
  }
}

// The coordinates operating_collar hands to billed_outputs_at must be the ones billed_outputs_at
// would have solved for itself, so handing them over moves no number and no row.
// Asserted bit for bit on the tangent, where both routes run the same arithmetic
// on the same values; the adjoint, where the rows are summed in another order,
// is test_transpose's.
struct HandoffTally {
  int interior = 0, other = 0, rows = 0, differ = 0, misplaced = 0;
};

// The curves the surface serves beyond TF24 and TF24_floor, for the tests that
// treat those two separately.
template <class F>
void for_each_further_surface_curve(F&& f) {
  using CC = phylloptim::Leaf::CostCurve;
  f(std::integral_constant<CC, CC::CF77>{});
  f(std::integral_constant<CC, CC::JS22>{});
  f(std::integral_constant<CC, CC::CMax>{});
  f(std::integral_constant<CC, CC::SOX>{});
  f(std::integral_constant<CC, CC::JW26>{});
  f(std::integral_constant<CC, CC::LeastCost>{});
  f(std::integral_constant<CC, CC::ProfitMax>{});
}

template <phylloptim::Leaf::CostCurve K>
void compare_held_handoff(phylloptim::Leaf &l, const std::string &at,
                          HandoffTally &tally) {
  using T = odelia::ode::tangent_scalar<double>;
  using Kind = phylloptim::Leaf::OperatingPointKind;
  const Kind kind = l.operating_point_kind();
  const double curvature = kind == Kind::Interior
                               ? l.profit_curvature_at_collar<K>()
                               : std::numeric_limits<double>::quiet_NaN();
  if (kind == Kind::Interior &&
      !(curvature < -phylloptim::Leaf::interior_curvature_floor)) {
    return;
  }
  if (!fixture::surface_has_objective(l)) return;
  ++(kind == Kind::Interior ? tally.interior : tally.other);
  // Traits every curve reads, the shallowest layer's soil potential (-1), and
  // the curve's own parameters, with TF24's for TF24_floor, whose hydraulic half
  // they are.
  std::vector<int> directions{phylloptim::par_vcmax_25, phylloptim::par_stem_P50,
                              phylloptim::par_root_P50, phylloptim::par_kmax, -1};
  for (const auto& q : phylloptim::Curve<K>::parameters) directions.push_back(q.slot);
  if (K == phylloptim::Leaf::CostCurve::TF24_floor) {
    for (const auto& q : phylloptim::Curve<phylloptim::Leaf::CostCurve::TF24>::parameters) {
      directions.push_back(q.slot);
    }
  }
  for (int direction : directions) {
    phylloptim::leaf_pars<T> pars;
    const phylloptim::leaf_pars<double> seated = l.passive_pars();
    for (std::size_t i = 0; i < pars.size(); ++i) pars[i] = T(seated[i]);
    std::vector<T> soil, r_h, r_v;
    for (int i = 0; i < l.supply_n_layers(); ++i) {
      soil.push_back(T(l.roots_.psi_soil_[std::size_t(i)]));
      r_h.push_back(T(l.roots_.network_.r_R_H_min[std::size_t(i)]));
      r_v.push_back(T(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
    }
    odelia::ode::set_direction(direction < 0 ? soil[0] : pars[std::size_t(direction)],
                               1.0);
    const phylloptim::SupplyInputs<T> supply{soil, r_h, r_v,
                                         pars[phylloptim::par_root_P50],
                                         pars[phylloptim::par_root_c]};
    const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), supply);
    const T solved_collar = l.operating_collar<K, T>(draw, pars, curvature);
    const auto solved = l.billed_outputs_at<K, T>(solved_collar, draw, pars);
    const auto solved_p = l.physical_outputs_at<K, T>(solved_collar, draw, pars);
    phylloptim::Leaf::HeldCoords<T> held;
    const T handed_collar = l.operating_collar<K, T>(draw, pars, curvature, &held);
    const auto handed = l.billed_outputs_at<K, T>(handed_collar, draw, pars, &held);
    const auto handed_p =
        l.physical_outputs_at<K, T>(handed_collar, draw, pars, &held);
    if (held.placed() != (kind == Kind::Interior)) {
      ++tally.misplaced;
      printf("    placed=%d at %s\n", int(held.placed()), at.c_str());
    }
    auto same = [](const T &a, const T &b) {
      return odelia::util::to_passive(a) == odelia::util::to_passive(b) &&
             odelia::ode::derivative_along(a) == odelia::ode::derivative_along(b);
    };
    bool all_same = same(solved_collar, handed_collar) &&
                    same(solved.profit, handed.profit) &&
                    same(solved.carbon, handed.carbon) &&
                    same(solved_p.carbon, handed_p.carbon) &&
                    same(solved_p.A, handed_p.A) && same(solved_p.E, handed_p.E) &&
                    same(solved_p.psi_leaf, handed_p.psi_leaf) &&
                    same(solved_p.ci, handed_p.ci) && same(solved_p.gc, handed_p.gc) &&
                    solved.uptake.size() == handed.uptake.size();
    for (std::size_t j = 0; all_same && j < solved.uptake.size(); ++j) {
      all_same = same(solved.uptake[j], handed.uptake[j]);
    }
    ++tally.rows;
    if (!all_same) {
      ++tally.differ;
      printf("    differs along %d at %s: profit row %.17g vs %.17g\n", direction,
             at.c_str(), odelia::ode::derivative_along(solved.profit),
             odelia::ode::derivative_along(handed.profit));
    }
  }
}

void test_held_coordinates_move_no_row() {
  printf("the coordinates operating_collar hands to billed_outputs_at move no number and no row\n");
  using K = phylloptim::Leaf::CostCurve;
  HandoffTally tally;
  for (int layers : {1, 3}) {
    for (double psi0 : {0.5, 1.0, 2.0, 4.0, 6.0}) {
      for (double ppfd : {100.0, 900.0}) {
        for (double temp : {25.0, 40.0}) {
          Drivers d;
          d.PPFD = ppfd;
          d.leaf_temp = temp;
          std::vector<double> psi_soil, soil_depth;
          for (int i = 0; i < layers; ++i) {
            psi_soil.push_back(psi0 + 0.35 * i);
            soil_depth.push_back(1.0 * (i + 1));
          }
          const std::string at = "psi_soil=" + std::to_string(psi0) +
                                 " ppfd=" + std::to_string(ppfd) +
                                 " T=" + std::to_string(temp) +
                                 " layers=" + std::to_string(layers);
          phylloptim::Leaf tf = make_leaf(d, psi_soil, soil_depth);
          tf.find_root_collar_psi();
          compare_held_handoff<K::TF24>(tf, "TF24 " + at, tally);
          // The floor's price reads the collar inside profit, so it is the curve
          // whose profit could tell a held collar from a live one.
          phylloptim::Leaf fl = make_leaf(d, psi_soil, soil_depth);
          fl.TF24_floor_lambda_o = kPriceForFOC;
          fl.set_model(K::TF24_floor, true);
          fl.optimise();
          compare_held_handoff<K::TF24_floor>(fl, "TF24_floor " + at, tally);
          for_each_further_surface_curve([&](auto tag) {
            constexpr K curve = decltype(tag)::value;
            phylloptim::Leaf c = make_leaf(d, psi_soil, soil_depth);
            c.set_model(curve, true);
            try {
              c.optimise();
            } catch (const phylloptim::util::infeasible_error&) {
              return;  // ProfitMax with no stream: see outputs_agree_with_the_solve_on
            }
            compare_held_handoff<curve>(
                c, phylloptim::Leaf::curve_name(static_cast<int>(curve)) + " " + at,
                tally);
          });
        }
      }
    }
  }
  printf("    %d rows | %d interior points, %d at other kinds\n", tally.rows,
         tally.interior, tally.other);
  ok(tally.interior > 0 && tally.other > 0,
     "the sweep reaches interior points and other kinds (test is not vacuous)");
  ok(tally.misplaced == 0, "operating_collar places coordinates at interior points only");
  ok(tally.differ == 0,
     "every collar, profit and uptake, value and row, is the same with the "
     "coordinates handed over");
}

// Coordinates are handed over only where they mean something, and refused,
// by name, where reading them would answer for another point.
void test_held_coordinates_are_refused_off_their_point() {
  printf("held coordinates are refused away from the point that placed them\n");
  using T = odelia::ode::tangent_scalar<double>;
  using CC = phylloptim::Leaf::CostCurve;
  using Kind = phylloptim::Leaf::OperatingPointKind;
  // pars and supply at the leaf's seated state, carrying no direction.
  struct Seat {
    phylloptim::leaf_pars<T> pars;
    std::vector<T> soil, r_h, r_v;
    phylloptim::SupplyInputs<T> supply() const {
      return {soil, r_h, r_v, pars[phylloptim::par_root_P50],
              pars[phylloptim::par_root_c]};
    }
  };
  auto seat = [](const phylloptim::Leaf &l) {
    Seat s;
    const phylloptim::leaf_pars<double> p = l.passive_pars();
    for (std::size_t i = 0; i < p.size(); ++i) s.pars[i] = T(p[i]);
    for (int i = 0; i < l.supply_n_layers(); ++i) {
      s.soil.push_back(T(l.roots_.psi_soil_[std::size_t(i)]));
      s.r_h.push_back(T(l.roots_.network_.r_R_H_min[std::size_t(i)]));
      s.r_v.push_back(T(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
    }
    return s;
  };
  Drivers d;
  phylloptim::Leaf in = make_leaf(d, {1.0}, {1.0});
  in.find_root_collar_psi();
  // Shade death sits on the wet bound, and dim light is what reaches it.
  Drivers dim;
  dim.PPFD = 1.0;
  phylloptim::Leaf bound = make_leaf(dim, {1.0}, {1.0});
  bound.find_root_collar_psi();
  ok(in.operating_point_kind() == Kind::Interior, "the first fixture is interior");
  ok(bound.operating_point_kind() != Kind::Interior &&
         bound.operating_point_kind() != Kind::HydraulicShutdown,
     std::string("the second fixture is at a bound, got ") +
         bound.operating_point_kind_name(bound.operating_point_kind()));

  const Seat si = seat(in);
  const auto si_supply = si.supply();
  const auto draw = in.supply_draw_at<T>(T(in.opt_root_psi_), si_supply);
  const double curvature = in.profit_curvature_at_collar<CC::TF24>();
  phylloptim::Leaf::HeldCoords<T> held;
  const T collar = in.operating_collar<CC::TF24, T>(draw, si.pars, curvature, &held);
  ok(held.placed() && held.at == in.opt_root_psi_,
     "the interior arm places the coordinates at its own collar");

  // Read at another collar on the same leaf: the draw is taken there, so only
  // the coordinates are out of place.
  const double off = in.opt_root_psi_ + 0.01;
  const auto draw_off = in.supply_draw_at<T>(T(off), si_supply);
  ok(throws_with([&] { in.billed_outputs_at<CC::TF24, T>(T(off), draw_off, si.pars, &held); },
                 "Leaf::billed_outputs_at: the held coordinates were placed at a collar of"),
     "billed_outputs_at refuses coordinates placed at another collar");

  // Read at a point of another kind.
  const Seat sb = seat(bound);
  const auto sb_supply = sb.supply();
  const auto draw_b = bound.supply_draw_at<T>(T(bound.opt_root_psi_), sb_supply);
  ok(throws_with([&] {
       bound.billed_outputs_at<CC::TF24, T>(T(bound.opt_root_psi_), draw_b, sb.pars, &held);
     },
                 "were placed at an interior point and are being read at one that is"),
     "billed_outputs_at refuses coordinates placed at an interior point at another kind");

  // operating_collar unplaces what it is handed before deciding, so coordinates left
  // from the interior point cannot survive a call at the bound.
  bound.operating_collar<CC::TF24, T>(draw_b, sb.pars,
                                      std::numeric_limits<double>::quiet_NaN(), &held);
  ok(!held.placed(), "operating_collar at a bound leaves the coordinates unplaced");
  ok(std::isfinite(odelia::util::to_passive(
         bound.billed_outputs_at<CC::TF24, T>(T(bound.opt_root_psi_), draw_b, sb.pars,
                                              &held).profit)),
     "billed_outputs_at solves its own with unplaced coordinates");

  // At double the interior arm records nothing, so there is nothing to hand over.
  phylloptim::leaf_pars<double> pd = in.passive_pars();
  std::vector<double> soil_d(in.roots_.psi_soil_), rh_d(in.roots_.network_.r_R_H_min),
      rv_d(in.roots_.network_.r_R_V_sum);
  const phylloptim::SupplyInputs<double> supply_d{soil_d, rh_d, rv_d,
                                              pd[phylloptim::par_root_P50],
                                              pd[phylloptim::par_root_c]};
  const auto draw_d = in.supply_draw_at<double>(in.opt_root_psi_, supply_d);
  phylloptim::Leaf::HeldCoords<double> held_d;
  const double collar_d = in.operating_collar<CC::TF24, double>(draw_d, pd, curvature, &held_d);
  ok(!held_d.placed(), "operating_collar at double leaves the coordinates unplaced");
  ok(in.billed_outputs_at<CC::TF24, double>(collar_d, draw_d, pd, &held_d).profit ==
         in.profit_,
     "billed_outputs_at at double with unplaced coordinates is the solve's profit");
  (void)collar;
}

// What plant reads, against what the solve itself left.
//
// billed_outputs_at and dprofit_dcollar_at ASSEMBLE the leaf's answers from kernels, where the
// solve reaches the same numbers by a different route entirely -- through
// set_leaf_states_rates and the members. So agreement is a statement about the
// assembly rather than a tautology, and it is asserted BIT-FOR-BIT because
// nothing here is reassociated: the assembly is meant to be the same arithmetic.
//
// ⚠️ EVERY KIND, because the branches differ. An interior point reads the held
// collar (the envelope omission), a bound reads its own condition, and the two
// shut leaf is its own case: it sits on the wet bound with the stem at the
// collar, where the stomatal flux is zero but the rows are the bound's own. The
// census below is
// part of the result -- a fixture reaching only one of them has narrowed the
// check without narrowing what it claims.
// The surface's outputs at S = double against the numbers the solve reported:
// every value the surface hands over is the solve's, so its rows are rows of
// the right number. Bit for bit, except carbon on a curve whose price reads only
// the flux (see below) and gc (likewise).
template <phylloptim::Leaf::CostCurve K, bool EB = false>
void outputs_agree_with_the_solve_on(double price) {
  namespace pl = phylloptim;
  using Kind = pl::Leaf::OperatingPointKind;
  const std::string curve =
      pl::Leaf::curve_name(static_cast<int>(K)) + (EB ? " with the energy balance" : "");
  const std::vector<std::vector<double>> soils = {
      {2.0}, {1.0, 2.0, 3.0}, {0.5, 1.0, 1.5, 2.0, 2.5}};
  const std::vector<std::vector<double>> depths = {
      {1.0}, {0.5, 1.0, 1.5}, {0.4, 0.8, 1.2, 1.6, 2.0}};
  int n_interior = 0, n_bound = 0, n_shut = 0, n_marg = 0;
  // With the balance on, dim light as well: where LeastCost leaves psi_crit.
  std::vector<double> lights{Drivers().PPFD, 30.0, 1.0};
  if (EB) lights.push_back(100.0);
  for (double ppfd : lights) for (double T : {25.0, 40.0}) {
    for (double psi0 : {0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 17.0}) {
      for (std::size_t k = 0; k < soils.size(); ++k) {
        Drivers d;
        d.leaf_temp = T;
        d.PPFD = ppfd;
        std::vector<double> ps;
        bool bad = false;
        for (double v : soils[k]) {
          ps.push_back(v + psi0 - 1.0);
          if (ps.back() <= 0.0) bad = true;
        }
        if (bad) continue;
        // ProfitMax's thermal cost on with the balance: it reads leaf temperature.
        pl::Leaf l = EB ? make_pm_leaf(d, ps, depths[k], true, 2.0, 0.05,
                                       K == pl::Leaf::CostCurve::ProfitMax)
                        : make_leaf(d, ps, depths[k]);
        l.CF77_lambda_ = price;
        l.TF24_floor_lambda_o = price;
        l.set_model(K, true);
        // A solve the supply refuses as infeasible leaves nothing to compare.
        try {
          l.optimise();
        } catch (const pl::util::infeasible_error&) {
          continue;
        }
        const auto kind = l.operating_point_kind();
        if (kind == Kind::SolverRefused || kind == Kind::NonFiniteGradient ||
            kind == Kind::Unsolved) {
          continue;
        }
        if (!fixture::surface_has_objective(l)) continue;
        if (kind == Kind::Interior) ++n_interior;
        else if (kind == Kind::HydraulicShutdown) ++n_shut;
        else ++n_bound;
        const std::string at = " (" + curve + ", " +
                               l.operating_point_kind_name(kind) + ")";

        const auto pars = l.passive_pars();
        const double p0 = l.opt_root_psi_;
        const auto held_supply = l.roots_.held_supply();
        const auto draw = l.supply_draw_at<double>(p0, held_supply);
        const auto got = l.billed_outputs_at<K, double>(p0, draw, pars);
        const auto phys = l.physical_outputs_at<K, double>(p0, draw, pars);
        ok(got.profit == l.profit_, "billed_outputs_at's profit IS the solve's" + at);
        ok(phys.profit == l.profit_, "physical_outputs_at's profit IS the solve's" + at);
        ok(got.uptake.size() == l.soil_consumption_.size(),
           "one uptake per layer" + at);
        for (std::size_t j = 0; j < got.uptake.size(); ++j) {
          ok(got.uptake[j] == l.soil_consumption_[j],
             "billed_outputs_at's layer draw IS the solve's" + at);
        }
        if constexpr (pl::Curve<K>::price_reads_flux_only) {
          // The solve's carbon is profit plus the price (Leaf::carbon); the
          // physical outputs' is A less the carbon part, equal to rounding.
          near(phys.carbon, l.carbon(), 1e-12, "the physical carbon is the solve's" + at);
          near(got.carbon, l.carbon(), 1e-10,
               "the billed carbon is the solve's, to the collar's tolerance" + at);
        } else {
          ok(phys.carbon == l.carbon(), "the physical carbon IS the solve's" + at);
          ok(got.carbon == l.carbon(), "the billed carbon IS the solve's" + at);
        }
        ok(phys.A == l.assim_colimited_, "A IS the solve's" + at);
        ok(phys.E == l.transpiration_, "E IS the solve's" + at);
        ok(phys.psi_leaf == l.opt_psi_leaf_, "psi_leaf IS the solve's" + at);
        ok(phys.ci == l.ci_ || (std::isnan(phys.ci) && std::isnan(l.ci_)),
           "ci IS the solve's" + at);
        // The marginal's association of the conductance, which the reported
        // one forms in another order: equal to rounding, not bit for bit.
        near(phys.gc, l.stom_cond_CO2_, 1e-14, "gc is the solve's" + at);
        if (kind == Kind::Interior) {
          bool feasible = false;
          const double up = l.dprofit_at_collar_psi<K>(p0, &feasible);
          // dprofit_at_collar_psi drives the model to that collar; put the point
          // back before reading anything else off this leaf, because every path
          // out of the solve writes its own rates.
          l.replay_operating_point(p0, kind);
          if (feasible && std::isfinite(up)) {
            ++n_marg;
            const auto held2 = l.roots_.held_supply();
            const auto d2 = l.supply_draw_at<double>(p0, held2);
            ok(l.dprofit_dcollar_at<K, double>(p0, d2, pars,
                                               l.normalisers_at<K, double>(d2, pars)) == up,
               "dprofit_dcollar_at IS upstream's dprofit_at_collar_psi" + at);
          }
        }
      }
    }
  }
  printf("    %s: %d interior (%d refereed against dprofit), %d at a bound, %d shut\n",
         curve.c_str(), n_interior, n_marg, n_bound, n_shut);
  const auto reach = fixture::surface_kinds(K);
  const bool shuts = std::find(reach.begin(), reach.end(),
                               Kind::HydraulicShutdown) != reach.end();
  ok(n_interior > 0 && n_bound > 0 && (n_shut > 0 || !shuts),
     curve + ": the sweep reached an interior point, a bound and, where the "
             "curve has one, a shut point");
  ok(n_marg > 0, curve + ": the sweep refereed the marginal somewhere");
}

// Each curve's scalar forms against its double forms, at S = double: the surface
// bills what the solve optimised, factor for factor. Bit for bit: each scalar
// form reads the pack's primitives, which equal the members the double form
// reads (passive_pars), in the double form's operand order.
void test_scalar_forms_are_the_double_forms() {
  printf("each curve's scalar forms are its double forms at S = double\n");
  namespace pl = phylloptim;
  using CC = pl::Leaf::CostCurve;
  Drivers d;
  pl::Leaf base = make_leaf(d, {1.0, 1.25, 1.5}, {0.5, 1.0, 1.5});
  base.optimise();
  const pl::leaf_pars<double> p = base.passive_pars();
  // ProfitMax's double forms read the normalisers its solve seeds.
  pl::Leaf seeded = base;
  seeded.set_model(CC::ProfitMax, true);
  seeded.optimise();
  std::map<std::string, int> compared;
  for (int c = 0; c < pl::Leaf::n_cost_curves; ++c) {
    pl::Leaf::with_curve(static_cast<CC>(c), [&](auto tag) {
      constexpr CC K = decltype(tag)::value;
      using C = pl::Curve<K>;
      if constexpr (pl::SurfaceForms<C>) {
        const std::string name = pl::Leaf::curve_name(c);
        for (double collar : {0.5, 1.0, 2.0, 3.0}) {
          for (double drop : {0.05, 0.5, 1.0, 2.0}) {
            const double sigma = collar + drop;
            if (C::cost_domain_ends_at_psi_crit && sigma >= base.psi_crit) continue;
            pl::Leaf l = C::prepares_normalisers ? seeded : base;
            const std::string at = " (" + name + ", collar " + std::to_string(collar) +
                                   ", leaf " + std::to_string(sigma) + ")";
            if constexpr (C::shape == pl::ObjectiveShape::Difference) {
              ok(C::template scalar_cost<double>(l, sigma, collar, p) ==
                     C::cost(l, sigma, collar),
                 "the scalar cost is the cost" + at);
            }
            if constexpr (C::shape == pl::ObjectiveShape::Normalised) {
              const pl::Leaf::Normalisers<double> n{
                  l.profitmax_k_soil_, l.profitmax_k_span_, l.profitmax_A_max_};
              ok(C::template scalar_hydraulic_cost_norm<double>(l, sigma, p, n) ==
                     C::hydraulic_cost_norm(l, sigma),
                 "the scalar normalised cost is the normalised cost" + at);
              ok(C::template scalar_cost_slope<double>(l, sigma, collar, p, n) ==
                     C::cost_slope(l, sigma, collar),
                 "the scalar slope is the slope" + at);
            } else {
              ok(C::template scalar_cost_slope<double>(l, sigma, collar, p) ==
                     C::cost_slope(l, sigma, collar),
                 "the scalar slope is the slope" + at);
            }
            if constexpr (C::cost_reads_upstream) {
              ok(C::template scalar_cost_slope_upstream<double>(l, sigma, collar, p) ==
                     C::cost_slope_upstream(l, sigma, collar),
                 "the scalar upstream slope is the upstream slope" + at);
            }
            ok(C::template scalar_carbon_cost<double>(l, sigma, p) ==
                   C::carbon_cost(l, sigma),
               "the scalar carbon part is the carbon part" + at);
            if constexpr (C::price_reads_flux_only) {
              const double E = l.transpiration(sigma, collar);
              ok(C::template scalar_price<double>(l, E, p) == C::price(l, E),
                 "the scalar price is the price" + at);
            }
            ++compared[name];
          }
        }
      }
    });
  }
  std::string seen;
  for (const auto& [name, n] : compared) seen += " " + name + " " + std::to_string(n);
  printf("    points compared per curve:%s\n", seen.c_str());
  // Every curve has them, and every curve was compared on the whole grid.
  for (int c = 0; c < pl::Leaf::n_cost_curves; ++c) {
    const std::string name = pl::Leaf::curve_name(c);
    ok(compared[name] > 0, "the scalar forms were compared on " + name);
  }
}

// TF24_floor at lambda_o = 0 IS TF24 on the surface too, bit for bit: every
// value and every row of every output, at every kind the sweep reaches, in every
// parameter but the price, whose row is the flux on TF24_floor and absent on
// TF24. It holds because each scalar form adds the floor's term to its parent's
// expression, and a term of zero adds an exact zero (hazards.md, "One solver").
void test_surface_reduces_tf24_floor_to_tf24() {
  printf("the surface on TF24_floor at a zero price is TF24's, bit for bit\n");
  namespace pl = phylloptim;
  using CC = pl::Leaf::CostCurve;
  using T = odelia::ode::tangent_scalar<double>;
  using Kind = pl::Leaf::OperatingPointKind;
  std::map<std::string, int> seen;
  int compared = 0, differ = 0;
  // A grid over light and soil, and the referee's dry pin, which the grid's
  // drivers do not reach; each with and without a cuticle, set before the
  // physiology seats it.
  std::vector<pl::Leaf> bases;
  for (double g_min : {0.0, 5e-4}) {
    for (double ppfd : {1.0, 900.0, 1500.0}) {
      for (double psi0 : {1.0, 2.0, 4.0, 6.0, 17.0}) {
        Drivers d;
        pl::Leaf l = fixture::multilayer_leaf();
        l.setup_transpiration(100);
        l.setup_root_vulnerability(100);
        l.g_min = g_min;
        l.set_physiology(fixture::root_network({1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
                                               {0.5, 1.0, 1.5}),
                         ppfd, {psi0, psi0 + 0.25, psi0 + 0.5}, {0.5, 1.0, 1.5},
                         d.K_s * d.theta / d.h, d.atm_vpd, d.ca, d.leaf_temp,
                         d.atm_o2_kpa, d.atm_kpa);
        bases.push_back(l);
      }
    }
    pl::Leaf dry = fixture::multilayer_leaf();
    dry.g_min = g_min;
    dry.set_physiology(fixture::root_network({1.0 / 0.05}, {1.0}), 500.0, {6.0},
                       {1.0}, 0.000157 / 5.0, 0.5, 40.0, 25.0, 21.0, 101.3);
    bases.push_back(dry);
  }
  for (const pl::Leaf& base : bases) {
    pl::Leaf a = base, b = base;
    a.set_model(CC::TF24, true);
    b.TF24_floor_lambda_o = 0.0;
    b.set_model(CC::TF24_floor, true);
    a.optimise();
    b.optimise();
    const Kind kind = a.operating_point_kind();
    if (kind != Kind::Interior && kind != Kind::BoundarySoil &&
        kind != Kind::BoundaryDry && kind != Kind::HydraulicShutdown) {
      continue;
    }
    ok(b.operating_point_kind() == kind && b.opt_root_psi_ == a.opt_root_psi_ &&
           b.profit_ == a.profit_,
       "the solves reduce, so the surfaces read the same point");
    const double ca = kind == Kind::Interior
                          ? a.profit_curvature_at_collar<CC::TF24>()
                          : std::numeric_limits<double>::quiet_NaN();
    const double cb = kind == Kind::Interior
                          ? b.profit_curvature_at_collar<CC::TF24_floor>()
                          : std::numeric_limits<double>::quiet_NaN();
    if (kind == Kind::Interior && !(ca < -pl::Leaf::interior_curvature_floor)) {
      continue;
    }
    ok(ca == cb || (std::isnan(ca) && std::isnan(cb)),
       "the curvature reduces bit for bit");
    ++seen[a.operating_point_kind_name(kind)];
    for (int k = 0; k < pl::n_leaf_pars; ++k) {
      if (k == pl::par_TF24_floor_lambda_o) continue;
      const auto rows = [&]<CC K>(pl::Leaf& l, double curvature) {
        pl::leaf_pars<T> tp;
        const auto p = l.passive_pars();
        for (std::size_t i = 0; i < tp.size(); ++i) tp[i] = T(p[i]);
        odelia::ode::set_direction(tp[std::size_t(k)], 1.0);
        std::vector<T> soil, r_h, r_v;
        for (int i = 0; i < l.supply_n_layers(); ++i) {
          soil.push_back(T(l.roots_.psi_soil_[std::size_t(i)]));
          r_h.push_back(T(l.roots_.network_.r_R_H_min[std::size_t(i)]));
          r_v.push_back(T(l.roots_.network_.r_R_V_sum[std::size_t(i)]));
        }
        const pl::SupplyInputs<T> supply{soil, r_h, r_v, tp[pl::par_root_P50],
                                         tp[pl::par_root_c]};
        const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), supply);
        pl::Leaf::HeldCoords<T> held;
        const T collar = l.operating_collar<K, T>(draw, tp, curvature, &held);
        const auto o = l.billed_outputs_at<K, T>(collar, draw, tp, &held);
        const auto q = l.physical_outputs_at<K, T>(collar, draw, tp, &held);
        std::vector<double> v{collar.value(), o.profit.value(), o.carbon.value(),
                              q.carbon.value(), q.A.value(), q.E.value(),
                              q.psi_leaf.value(), q.ci.value(), q.gc.value(),
                              o.residual_flux.value()};
        for (const T* x : {&collar, &o.profit, &o.carbon, &q.carbon, &q.A, &q.E,
                           &q.psi_leaf, &q.ci, &q.gc, &o.residual_flux}) {
          v.push_back(odelia::ode::derivative_along(*x));
        }
        for (const T& u : o.uptake) {
          v.push_back(u.value());
          v.push_back(odelia::ode::derivative_along(u));
        }
        return v;
      };
      const auto ra = rows.template operator()<CC::TF24>(a, ca);
      const auto rb = rows.template operator()<CC::TF24_floor>(b, cb);
      bool same = ra.size() == rb.size();
      for (std::size_t j = 0; same && j < ra.size(); ++j) {
        same = ra[j] == rb[j] || (std::isnan(ra[j]) && std::isnan(rb[j]));
      }
      if (!same) {
        ++differ;
        ok(false, std::string("the surface reduces at ") + a.operating_point_kind_name(kind) +
                      " in " + pl::gradient::par_names()[std::size_t(k)]);
      }
      ++compared;
    }
  }
  std::string kinds;
  for (const auto& [kind, n] : seen) kinds += " " + kind + " " + std::to_string(n);
  printf("    %d parameter directions compared at%s, %d differ\n", compared,
         kinds.c_str(), differ);
  ok(seen.size() == 4, "the reduction was compared at all four kinds");
}

// A ProfitMax draw keeps its own copy of the supply, so its normalisers read
// the same numbers after the caller's supply is gone; a draw on any other
// curve keeps none.
void test_a_profitmax_draw_outlives_its_supply() {
  printf("a ProfitMax draw keeps the supply its normalisers read\n");
  namespace pl = phylloptim;
  using CC = pl::Leaf::CostCurve;
  Drivers d;
  pl::Leaf l = make_leaf(d, {1.0, 1.25, 1.5}, {0.5, 1.0, 1.5});
  l.set_model(CC::ProfitMax, true);
  l.optimise();
  const auto pars = l.passive_pars();
  const auto drawn_in_a_scope = [&] {
    const auto supply = l.roots_.held_supply();
    return l.supply_draw_at<double>(l.opt_root_psi_, supply);
  };
  const auto draw = drawn_in_a_scope();
  const auto held = l.roots_.held_supply();
  const auto live = l.supply_draw_at<double>(l.opt_root_psi_, held);
  const auto a = l.normalisers_at<CC::ProfitMax, double>(draw, pars);
  const auto b = l.normalisers_at<CC::ProfitMax, double>(live, pars);
  ok(draw.inputs != nullptr && a.k_soil == b.k_soil && a.k_span == b.k_span &&
         a.A_max == b.A_max,
     "the normalisers read the draw's own copy of the supply");
  pl::Leaf t = make_leaf(d, {1.0, 1.25, 1.5}, {0.5, 1.0, 1.5});
  t.optimise();
  const auto t_supply = t.roots_.held_supply();
  ok(t.supply_draw_at<double>(t.opt_root_psi_, t_supply).inputs == nullptr,
     "a curve without normalisers keeps no copy");
}

void test_outputs_agree_with_the_solve() {
  printf("the surface's outputs and marginal against the solve's own numbers\n");
  using K = phylloptim::Leaf::CostCurve;
  outputs_agree_with_the_solve_on<K::TF24>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::TF24_floor>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::CF77>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::JS22>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::CMax>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::SOX>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::JW26>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::LeastCost>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::ProfitMax>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::TF24, true>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::TF24_floor, true>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::CF77, true>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::JS22, true>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::CMax, true>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::SOX, true>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::JW26, true>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::LeastCost, true>(kPriceForFOC);
  outputs_agree_with_the_solve_on<K::ProfitMax, true>(kPriceForFOC);
}

}  // namespace

void run_surface() {
  test_scalar_forms_are_the_double_forms();
  test_a_profitmax_draw_outlives_its_supply();
  test_outputs_agree_with_the_solve();
  test_surface_reduces_tf24_floor_to_tf24();
  test_pack_kernels_are_the_models_own();
  test_light_reaches_carbon_with_a_row();
  test_every_answered_point_hands_over_finite_rows();
  test_a_replay_hands_over_the_rows_its_solve_did();
  test_a_dry_pin_meets_the_interior_answer_continuously();
  test_held_coordinates_move_no_row();
  test_held_coordinates_are_refused_off_their_point();
  test_supplied_rows_match_a_resolve_at_every_kind();
  test_shut_leaf_seated_at_its_collar_bills_the_roots_nothing();
  test_surface_refuses_a_curve_or_supply_it_was_not_solved_on();
  test_tf24_floor_rows_at_a_price();
}
