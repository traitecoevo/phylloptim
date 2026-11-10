// -*-c++-*-
#ifndef PHYLLOPTIM_GRADIENT_HPP_
#define PHYLLOPTIM_GRADIENT_HPP_

// Trait gradients at a solved or prescribed operating point, from the
// derivative surface: one forward pass carrying a direction per slot of theta,
// and one on the route's decision variable, through supply_draw_at,
// operating_collar and physical_outputs_at. Every row is the implicit function
// theorem's at the converged point, on either route and supply path, with the
// energy balance and the thermal-damage layer. Every trait gradient the package
// gives is this one: leaf_gradient() is `batch` with one observation. From R,
// dispatch rather than the model dominates a gradient's cost, so this crosses
// the boundary once per likelihood evaluation, not once per parameter per
// observation. It returns dY/dtheta for the model parameters; the likelihood and
// the caller's parameterisation Jacobian stay with the caller.
//
// tests/testthat/test-gradient-referee.R checks the result against a
// difference of the solve for every curve on both routes; gradient_golden.tsv
// records it.

#include <phylloptim/leaf_model.hpp>
#include <phylloptim/roots.hpp>
#include <phylloptim/util.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <string>
#include <vector>

namespace phylloptim {
namespace gradient {

// How long `theta` is: one entry per par_names(), the slots a calibration
// proposes, which come first in the enumeration and stop at par_PPFD.
// ⚠️ Not n_leaf_pars, which also counts par_PPFD (hazard 12).
inline constexpr int n_theta = par_PPFD;

// The proposed slots' names: each curve's parameters from its definition
// (`curve_parameters`), the slots no curve owns here. par_PPFD is seated from
// the drivers and never fitted, so it has none. Append, never reorder (hazard
// 13); test-gradient-batch.R compares R's `.gradient_par_names()` with this list.
struct SlotName {
  int slot;
  const char* name;
};
inline constexpr SlotName slots_no_curve_owns[] = {
    {par_vcmax_25, "vcmax_25"},       {par_stem_c, "stem_c"},
    {par_stem_P50, "stem_P50"},       {par_root_c, "root_c"},
    {par_root_P50, "root_P50"},       {par_jmax_25, "jmax_25"},
    {par_phi_J, "phi_J"},             {par_curv_fact_elec_trans, "curv_fact_elec_trans"},
    {par_curv_fact_colim, "curv_fact_colim"}, {par_R_d_25, "R_d_25"},
    {par_g_min, "g_min"},             {par_kmax, "leaf_specific_conductance_max"},
    {par_resistance, "resistance"},   {par_A, "A_"},
    {par_I, "I_"},                    {par_D, "D_"},
    {par_topt_offset, "topt_offset_"}, {par_dTopt_max, "dTopt_max_"},
    {par_K_A, "K_A_"},                {par_k_i, "k_i_"},
    {par_k_i_baseline, "k_i_baseline_"}, {par_k_rec, "k_rec_"},
    {par_k_rec_Ea, "k_rec_Ea_"},      {par_c_acclim_maint, "c_acclim_maint_"},
    {par_c_repair_maint, "c_repair_maint_"},
    {par_c_protect_maint, "c_protect_maint_"},
    {par_c_repair_flux, "c_repair_flux_"}};

constexpr std::array<const char*, n_theta> par_names_by_slot() {
  std::array<const char*, n_theta> out{};
  for (const SlotName& s : slots_no_curve_owns) out[s.slot] = s.name;
  for (const CurveParameterEntry& e : curve_parameters) {
    out[e.parameter.slot] = e.parameter.field;
  }
  return out;
}

// Every slot is named exactly once, by a curve or by the list above.
constexpr bool every_slot_named_once() {
  std::array<int, n_theta> seen{};
  for (const SlotName& s : slots_no_curve_owns) ++seen[s.slot];
  for (const CurveParameterEntry& e : curve_parameters) ++seen[e.parameter.slot];
  for (int n : seen) {
    if (n != 1) return false;
  }
  return true;
}
static_assert(every_slot_named_once(),
              "a slot before par_PPFD is unnamed or named twice: a curve's slot "
              "goes in its Curve<K>::parameters, any other in slots_no_curve_owns");

inline const std::vector<std::string>& par_names() {
  static constexpr auto by_slot = par_names_by_slot();
  static const std::vector<std::string> names(by_slot.begin(), by_slot.end());
  return names;
}

// --- the five differentiated outputs ------------------------------------------
//
// A, gc, psi_leaf, collar and profit, the order R names them in (it reads
// `gradient_output_names()`): append, never reorder. `collar` is psi* itself, so dcollar/dtheta is
// dpsi*/dtheta. The first four are what a gas-exchange calibration observes;
// `profit` is the carbon plant's mass budget consumes.
inline constexpr int n_outputs = 5;
inline constexpr int out_psi_leaf = 2;
inline constexpr int out_collar = 3;
inline constexpr int out_profit = 4;

inline const std::vector<std::string>& output_names() {
  static const std::vector<std::string> names{"A", "gc", "psi_leaf", "collar",
                                              "profit"};
  return names;
}

// --- which model's optimum is being differentiated ---------------------------
//
// A model is a (curve, route) pair, and every curve runs on both routes. The
// collar route maximises over the root-collar potential; a leaf route over
// psi_leaf with the upstream potential pinned at psi_soil, which is a different
// model, not a restriction of the collar one (hazard 11). The default is the
// Leaf's: TF24_floor on the leaf route.
enum class Route { Collar, Leaf };

struct Model {
  Leaf::CostCurve curve = Leaf::CostCurve::TF24_floor;
  Route route = Route::Leaf;
  bool collar() const { return route == Route::Collar; }
};

// Seat the model once, at the top of a gradient. Seated on the prescribed path
// too, which never calls route_solve, because the entry points below read the
// seated curve and would otherwise differentiate the leaf's last model.
inline void route_seat(Leaf& l, const Model& m) {
  l.set_model(m.curve, m.collar());
}

inline void route_solve(Leaf& l, const Model& m) {
  if (m.collar()) {
    l.find_root_collar_psi();
  } else {
    l.optimise();
  }
}

inline double route_psi(const Leaf& l, const Model& m) {
  return m.collar() ? l.opt_root_psi_ : l.opt_psi_leaf_;
}

inline void route_evaluate(Leaf& l, const Model& m, double psi) {
  if (m.collar()) {
    l.evaluate_root_collar_psi(psi);
  } else {
    l.evaluate_psi_leaf_at(psi);
  }
}

// The same, false where the leaf refuses `psi` as outside its feasible interval.
// Only that refusal is caught: anything else is a failure of the evaluation and
// propagates.
inline bool route_evaluate_feasible(Leaf& l, const Model& m, double psi) {
  try {
    route_evaluate(l, m, psi);
  } catch (const util::infeasible_error& e) {
    const std::string token = util::infeasible_token("prescribed_potential");
    if (std::string(e.what()).compare(0, token.size(), token) != 0) {
      throw;
    }
    return false;
  }
  return true;
}

// Read straight off the members, which operating_point_values() copies too.
inline void outputs(const Leaf& l, double* y) {
  y[0] = l.assim_colimited_;
  y[1] = l.stom_cond_CO2_;
  y[2] = l.opt_psi_leaf_;
  y[3] = l.opt_root_psi_;
  y[4] = l.profit_;
}

// --- one observation's drivers ------------------------------------------------
//
// Everything set_physiology takes except kmax and resistance, which come out of
// `theta`. Held C++-side and resolved once, because handing a RootNetwork
// across the R boundary costs more than the gradient it carries.
struct Drivers {
  RootNetwork root_network;
  double PPFD = 0.0;
  std::vector<double> psi_soil;
  std::vector<double> soil_depth;
  double atm_vpd = 0.0;
  double ca = 0.0;
  double leaf_temp = 0.0;
  double atm_o2_kpa = 0.0;
  double atm_kpa = 0.0;
};

// What the operating point is, not whether the call worked, except `Error`, the
// batch's per-row failure (see `batch`). `Prescribed` and `Refused` end the
// caller-supplied-psi path: the collar was imposed, and for `Refused` it fell
// outside the feasible interval. Append only: a consumer holding the integer
// breaks quietly if a value moves.
//
// Read off the solve's `OperatingPointKind` (status_of), which records which
// condition placed the point: an interior optimum is interior, a bound of
// either kind or a shut leaf is pinned, and a point the solve could not place,
// ProfitMax with no objective, and an interior point flatter than
// Leaf::interior_curvature_floor are no-gradient.
enum class Status { Interior, Pinned, NoGradient, Error, Prescribed, Refused };

// No `default:`, so a new Status is a compiler diagnostic rather than silently
// labelled "error".
inline std::string status_name(Status s) {
  switch (s) {
    case Status::Interior:   return "interior";
    case Status::Pinned:     return "pinned";
    case Status::NoGradient: return "no-gradient";
    case Status::Prescribed: return "prescribed";
    case Status::Refused:    return "refused";
    case Status::Error:      return "error";
  }
  return "error";  // unreachable; silences -Wreturn-type on a bad cast
}

// No `default:`, so a new kind must say what status it is.
inline Status status_of(Leaf::OperatingPointKind kind) {
  switch (kind) {
    case Leaf::OperatingPointKind::Interior:          return Status::Interior;
    case Leaf::OperatingPointKind::BoundarySoil:
    case Leaf::OperatingPointKind::BoundaryDry:
    case Leaf::OperatingPointKind::HydraulicShutdown: return Status::Pinned;
    case Leaf::OperatingPointKind::Prescribed:        return Status::Prescribed;
    case Leaf::OperatingPointKind::Unsolved:
    case Leaf::OperatingPointKind::SolverRefused:
    case Leaf::OperatingPointKind::NonFiniteGradient: return Status::NoGradient;
  }
  return Status::NoGradient;  // unreachable
}

struct Settings {
  // Which model's optimum to differentiate: the Leaf's default, TF24_floor on
  // the leaf route, unless set.
  Model model;
};

// A collar potential the caller imposes in place of the solved one, and how it
// responds to each parameter. `dpsi_dtheta` is npars long in `pars` order, not
// n_theta; null means zero (the partial at fixed collar).
struct Prescribed {
  double psi = 0.0;
  const double* dpsi_dtheta = nullptr;
};

struct Result {
  // The five outputs the gradient is taken at.
  double value[n_outputs];
  // npars * n_outputs, parameter-major: d(output j)/d(pars[k]) at
  // [k * n_outputs + j].
  std::vector<double> grad;
  Status status = Status::Error;
  // The curvature of profit in the decision variable, at an interior or a
  // prescribed point.
  double H = util::na_value;
  // d2profit/dpsi dtheta (npars long) and dY/dpsi at fixed traits: kept, being
  // also the coefficients of a caller's own sensitivity ODE for psi, which
  // cannot be recovered from the rows. M is -H dpsi*/dtheta at an interior
  // point, the marginal's own rows at a prescribed one, and NA elsewhere.
  std::vector<double> M;
  double dY_dpsi[n_outputs];
  // The decision variable the outputs were evaluated at -- psi* when solved,
  // the requested one when prescribed (refused or not).
  double psi = util::na_value;
  std::string message;

  void reset(std::size_t npars) {
    for (int j = 0; j < n_outputs; ++j) {
      value[j] = util::na_value;
      dY_dpsi[j] = util::na_value;
    }
    grad.assign(npars * n_outputs, util::na_value);
    M.assign(npars, util::na_value);
    status = Status::Error;
    H = util::na_value;
    psi = util::na_value;
    message.clear();
  }
};

// --- pushing a parameter vector back onto the leaf ----------------------------
//
// ⚠️ Order is load-bearing (hazard 10): set_traits() returns the leaf to its
// just-constructed state, so the drivers are re-supplied after it, and
// set_physiology() re-derives vcmax_/jmax_/R_d_ behind the temperature cache
// set_traits() cleared. The other order runs the gradient at the first vcmax
// the object saw.
inline void apply(Leaf& l, const double* theta, const Drivers& d, bool single) {
  l.set_traits(theta[par_vcmax_25], theta[par_stem_c], theta[par_stem_P50],
               theta[par_root_c], theta[par_root_P50],
               theta[par_TF24_beta2], theta[par_jmax_25],
               theta[par_phi_J], theta[par_curv_fact_elec_trans],
               theta[par_curv_fact_colim], theta[par_TF24_cost_scale],
               theta[par_R_d_25], theta[par_JS22_gamma],
               theta[par_CMax_a], theta[par_CMax_b], theta[par_g_min]);
  // Not traits and nothing is derived from them, so written directly. Only a
  // set slot is written: on a model that does not own a slot it is NA, and
  // writing it would overwrite the leaf's aligned default.
  if (!std::isnan(theta[par_CF77_lambda])) {
    l.CF77_lambda_ = theta[par_CF77_lambda];
  }
  if (!std::isnan(theta[par_TF24_floor_lambda_o])) {
    l.TF24_floor_lambda_o = theta[par_TF24_floor_lambda_o];
  }
  if (!std::isnan(theta[par_LeastCost_beta])) {
    l.LeastCost_beta = theta[par_LeastCost_beta];
  }
  // The thermal-damage layer's inputs, ONLY with its gate on -- off, the slots
  // hold whatever a caller zero-filled and are never read -- and BEFORE the
  // drivers, because set_physiology is what seats the block that reads them.
  if (l.use_thermal_damage_) {
    for (const DamageSlot& ds : damage_slots) {
      l.*ds.member = theta[ds.slot];
    }
  }
  if (single) {
    // R's `series_resistance()`: one series resistance in `r_R_V_sum`, built
    // from theta, since `resistance` is a differentiated parameter.
    RootNetwork net;
    net.r_R_V_sum.assign(1, theta[par_resistance]);
    l.set_physiology(net, d.PPFD, d.psi_soil, d.soil_depth, theta[par_kmax],
                     d.atm_vpd, d.ca, d.leaf_temp, d.atm_o2_kpa, d.atm_kpa);
  } else {
    l.set_physiology(d.root_network, d.PPFD, d.psi_soil, d.soil_depth,
                     theta[par_kmax], d.atm_vpd, d.ca, d.leaf_temp,
                     d.atm_o2_kpa, d.atm_kpa);
  }
}

// --- the pass -----------------------------------------------------------------
//
// One forward pass whose lane k carries d/d(theta[k]) and whose last lane
// carries d/dpsi at fixed traits, psi the route's decision variable. Every slot
// gets a lane whatever `pars` asks for: one instantiation of the surface, at a
// cost that is flat in the parameter count.
inline constexpr int n_lanes = n_theta + 1;
inline constexpr int lane_decision = n_theta;
using Lanes = xad::fwd<double, n_lanes>::active_type;

// The rows at the seated point, which the solve (or a prescribed evaluation)
// has placed; `out.status` is set and is not NoGradient on entry.
template <Leaf::CostCurve K>
inline void surface_rows(Leaf& l, const int* pars, std::size_t npars,
                         bool single, const Prescribed* prescribed, Result& out) {
  const bool interior = out.status == Status::Interior;
  double H = std::numeric_limits<double>::quiet_NaN();
  if (interior) {
    H = l.profit_curvature_at_collar<K>();
    out.H = H;
    // operating_collar refuses it too: nearer zero every row is amplification.
    if (!(H < -Leaf::interior_curvature_floor)) {
      out.status = Status::NoGradient;
      return;
    }
  } else if (prescribed != nullptr) {
    // No row reads it: it is the other coefficient of a tracking caller's
    // ds/dt = k(M + H s), and may be any sign off the optimum.
    out.H = l.profit_curvature_at_collar<K>();
  }
  // The pack at the seated values, lane k on slot k. A slot nothing reads may
  // be unset (another curve's price, NIR off the energy balance).
  const leaf_pars<double> seated = l.passive_pars();
  leaf_pars<Lanes> p;
  for (std::size_t i = 0; i < p.size(); ++i) {
    p[i] = Lanes(std::isfinite(seated[i]) ? seated[i] : 0.0);
  }
  for (int k = 0; k < n_theta; ++k) {
    xad::derivative(p[std::size_t(k)])[std::size_t(k)] = 1.0;
  }
  // The soil and the root network are drivers; on the single path the series
  // resistance is the pack's, so its slot carries a row. Each keeps its own
  // length: a network is sized to its deepest rooted layer, which may lie
  // above the bottom of the soil profile.
  const SupplyInputs<double> held = l.held_supply();
  std::vector<Lanes> soil, r_h, r_v;
  for (double psi : held.psi_soil) soil.push_back(Lanes(psi));
  for (std::size_t i = 0; i < held.r_R_H_min.size(); ++i) {
    r_h.push_back(Lanes(held.r_R_H_min[i]));
    r_v.push_back(single ? p[par_resistance] : Lanes(held.r_R_V_sum[i]));
  }
  const SupplyInputs<Lanes> supply{soil, r_h, r_v, p[par_root_P50], p[par_root_c]};
  const double decision = l.decision_psi();
  const auto draw = l.supply_draw_at<Lanes>(Lanes(decision), supply);
  Leaf::HeldCoords<Lanes> placed;
  Lanes x = prescribed != nullptr ? Lanes(decision)
                                  : l.operating_collar<K, Lanes>(draw, p, H, &placed);
  // The decision variable's own lane, zero in value: every output's last lane
  // is then dY/dpsi at fixed traits.
  Lanes x_lane(decision);
  xad::derivative(x_lane)[std::size_t(lane_decision)] = 1.0;
  x = x + (x_lane - Lanes(decision));
  const auto o = l.physical_outputs_at<K, Lanes>(x, draw, p,
                                                 prescribed != nullptr ? nullptr : &placed);
  // A, gc, psi_leaf, collar, profit; a leaf route solves for no collar.
  const bool collar = l.model_route() == "collar";
  const Lanes* y[n_outputs] = {&o.A, &o.gc, &o.psi_leaf, collar ? &x : nullptr,
                               &o.profit};
  const auto lane = [](const Lanes* v, int k) {
    return v == nullptr ? std::numeric_limits<double>::quiet_NaN()
                        : xad::derivative(*v)[std::size_t(k)];
  };
  for (int j = 0; j < n_outputs; ++j) {
    out.dY_dpsi[j] = lane(y[j], lane_decision);
  }
  // At a prescribed point the decision variable moves only as the caller says.
  Lanes marginal;
  if (prescribed != nullptr) {
    const Leaf::Normalisers<Lanes> norms = l.normalisers_at<K, Lanes>(draw, p);
    marginal = l.dprofit_dcollar_at<K, Lanes>(Lanes(decision), draw, p, norms);
  }
  for (std::size_t k = 0; k < npars; ++k) {
    const int slot = pars[k];
    const double dpsi = prescribed == nullptr || prescribed->dpsi_dtheta == nullptr
                            ? 0.0
                            : prescribed->dpsi_dtheta[k];
    for (int j = 0; j < n_outputs; ++j) {
      out.grad[k * n_outputs + std::size_t(j)] =
          lane(y[j], slot) + (prescribed != nullptr ? out.dY_dpsi[j] * dpsi : 0.0);
    }
    if (interior) {
      out.M[k] = -H * xad::derivative(x)[std::size_t(slot)];
    } else if (prescribed != nullptr) {
      out.M[k] = xad::derivative(marginal)[std::size_t(slot)];
    }
  }
}

// --- one observation ---------------------------------------------------------
//
// Fills `out` progressively and throws on a failure of the surface, so a partly
// determined row keeps its value and status when `batch` catches.
inline void at(Leaf& l, const double* theta, const Drivers& d, bool single,
               const int* pars, std::size_t npars, const Settings& s,
               Result& out, const Prescribed* prescribed = nullptr) {
  out.reset(npars);
  apply(l, theta, d, single);
  // Seated here, not only by `batch`, because at() is a public entry point.
  route_seat(l, s.model);
  // The default solves; `prescribed` imposes the decision variable. A
  // prescribed psi outside the feasible interval is reported, not thrown: a fit
  // visits such points routinely, and `status` tells a batch which rows to
  // distrust. Everything but `psi` stays NA.
  if (prescribed != nullptr) {
    if (!route_evaluate_feasible(l, s.model, prescribed->psi)) {
      out.status = Status::Refused;
      out.psi = prescribed->psi;
      return;
    }
  } else {
    route_solve(l, s.model);
  }
  out.psi = route_psi(l, s.model);
  outputs(l, out.value);
  // A prescribed point is prescribed whatever kind its evaluation wrote: a
  // collar past the end of the stem curve seats the shut-down state.
  out.status = prescribed != nullptr ? Status::Prescribed
                                     : status_of(l.operating_point_kind());
  // ProfitMax with no light, or no stream to normalise over, has no objective.
  if (s.model.curve == Leaf::CostCurve::ProfitMax && !l.profitmax_has_objective()) {
    out.status = Status::NoGradient;
  }
  // A collar the caller chose may put the leaf potential in the stem curve's
  // tail, where inverting the flux for it divides by a conductance near
  // e^-40 and every row is rounding. The solving path places such a leaf on
  // its bound instead; a prescribed one has no rows to give.
  if (prescribed != nullptr && s.model.collar() &&
      l.stem_curve_in_tail(l.opt_psi_leaf_)) {
    out.status = Status::NoGradient;
  }
  if (out.status == Status::NoGradient) {
    return;
  }
  Leaf::with_curve(s.model.curve, [&](auto tag) {
    surface_rows<tag.value>(l, pars, npars, single, prescribed, out);
  });
}

// --- the batch ---------------------------------------------------------------
//
// `theta` is a column-major `theta_nrow` x n_theta matrix, as R hands one over:
// one row per observation, or one row shared by all.
//
// Per-row status, not an exception: a proposal during a fit will reach points
// the solve cannot handle, and that must cost those rows, not the likelihood
// evaluation. A failed row's gradient is all NA rather than partly filled, so a
// caller checking status per row cannot read a gradient with a hole in it; its
// value and H stay where they were determined.
//
// `psi` is null to solve, or one decision variable per observation to impose;
// `dpsi_dtheta` is then a column-major n x npars matrix, null meaning zero.
inline std::vector<Result> batch(Leaf& l, const double* theta,
                                std::size_t theta_nrow,
                                const std::vector<Drivers>& drivers,
                                bool single, const int* pars,
                                std::size_t npars, const Settings& s,
                                const double* psi = nullptr,
                                const double* dpsi_dtheta = nullptr) {
  const std::size_t n = drivers.size();
  std::vector<Result> out(n);
  double th[n_theta];
  std::vector<double> dpsi(npars);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t row = theta_nrow == 1 ? 0 : i;
    for (int j = 0; j < n_theta; ++j) {
      th[j] = theta[row + std::size_t(j) * theta_nrow];
    }
    Prescribed p;
    if (psi != nullptr) {
      p.psi = psi[i];
      if (dpsi_dtheta != nullptr) {
        for (std::size_t k = 0; k < npars; ++k) {
          dpsi[k] = dpsi_dtheta[i + k * n];
        }
        p.dpsi_dtheta = dpsi.data();
      }
    }
    try {
      at(l, th, drivers[i], single, pars, npars, s, out[i],
         psi == nullptr ? nullptr : &p);
    } catch (const std::exception& e) {
      out[i].status = Status::Error;
      out[i].message = e.what();
      out[i].grad.assign(npars * n_outputs, util::na_value);
    }
  }
  return out;
}

}  // namespace gradient
}  // namespace phylloptim

#endif
