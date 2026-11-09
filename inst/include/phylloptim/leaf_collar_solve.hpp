// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_COLLAR_SOLVE_HPP_
#define PHYLLOPTIM_LEAF_COLLAR_SOLVE_HPP_

// MASTER SOLVER: optimal root-collar (and stem) water potential
//
// find_root_collar_psi solves soil -> root -> stem -> leaf once per call and
// stores the operating point in opt_psi_leaf_, opt_root_psi_ and profit_. The
// inner level is CONTINUITY (find_root_psi on E_column: soil uptake equals the
// stem's flux); the outer is OPTIMISATION, a root-find of dprofit/dcollar = 0
// over [wet bound, collar_at_dry_end] (maximise_profit_over_collar). Where no
// collar in that interval moves water the leaf is shut (set_shut_state).
//
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>
#include <phylloptim/leaf_cost_curves.hpp>

namespace phylloptim {

// The continuity residual at collar potential `x`: soil uptake less the stem's
// flux to `psi_leaf`, zero where the continuum balances.
inline double Leaf::E_column(double x, const std::vector<double>& psi_soil, double psi_leaf) {

  E_from_Soil_to_Root_Collar(x, psi_soil);
  double E_root_to_leaf = transpiration(psi_leaf, x);
  return E_up_ - E_root_to_leaf;
}

// Soil uptake at collar potential `x`; zero at the zero-uptake collar.
inline double Leaf::E_column_zero(double x, const std::vector<double>& psi_soil) {

  E_from_Soil_to_Root_Collar(x, psi_soil);

  return E_up_;
}

// The collar root on [wettest_soil_layer, end] of E_column with the stem held
// at `end` (find_root_crit == 1) or of E_column_zero (0), to bracket_tol. Both
// are smooth and monotone there, and the caller's early exits guarantee the
// sign change, so a superlinear bracketing solver (TOMS748) is safe.
inline double Leaf::find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit,
                                  double f_wettest, double f_crit, double end) {
  if (find_root_crit == 1) {
    auto target = [&](double x) -> double {
      return E_column(x, psi_soil, end);
    };
    try {
      // The wet end of the final bracket, not its midpoint: the target rises
      // with the collar, and past the root the soil can deliver more than
      // G(inf), which no leaf potential carries and the inverse refuses.
      return util::uniroot_smooth_bracket(target, wettest_soil_layer, end,
                                          f_wettest, f_crit, bracket_tol,
                                          ci_niter)
          .first;
    } catch (const std::exception& e) {
      util::stop_infeasible("collar_bracket",
                 "find_root_psi(find_root_crit=1) failed: " + std::string(e.what()) +
                 "; min=" + util::to_string(wettest_soil_layer) +
                 "; max=" + util::to_string(end));
    }
  }

  auto target = [&](double x) -> double {
    return E_column_zero(x, psi_soil);
  };
  try {
    return util::uniroot_smooth(target, wettest_soil_layer, end, f_wettest,
                                f_crit, bracket_tol, ci_niter);
  } catch (const std::exception& e) {
    util::stop_infeasible("collar_bracket",
               "find_root_psi(find_root_crit=0) failed: " + std::string(e.what()) +
               "; min=" + util::to_string(wettest_soil_layer) +
               "; max=" + util::to_string(end));
  }

}

// Re-solves find_root_psi's 1e-4 MPa bound `x0` to 1e-15, for a bound that is
// the answer. Where a loose root-find stops moves with the parameters through
// its iterates, while the derivative rows (wet_bound_at) differentiate the root
// itself, so an unpolished pinned optimum's root-trait rows were 100x off a
// difference of the solve. Only the answering bound is polished: polishing
// every bound slows every solve.
inline double Leaf::polish_bracket_bound(double x0,
                                         const std::vector<double>& psi_soil,
                                         int find_root_crit) {
  // The wet bound's target is uptake less the residual flux (the uptake itself
  // at g_min = 0).
  auto target = [&](double x) -> double {
    if (find_root_crit == 1) {
      return E_column(x, psi_soil, dry_end_);
    }
    const double up = E_column_zero(x, psi_soil);
    return up - residual_flux(x, up);
  };
  // bound_root_tol, tighter than collar_root_tol (~6e-12 at a 5 MPa collar):
  // near the dry end the steep stem curve turns that into ~2e-11 of psi_leaf, a
  // 2e-6 row error to a finite difference at a 1e-6 step.
  const double f0 = target(x0);
  if (f0 == 0.0 || !std::isfinite(f0)) {
    return x0;
  }
  // Both targets rise with the collar, so f0 > 0 puts the root below x0: that
  // side is probed first, starting at the bracket's own tolerance at x0. The
  // continuity root's probes stop at the dry end, past which the collar is
  // drier than the stem and transpiration() refuses the reversed gradient.
  const double hi = find_root_crit == 1 ? dry_end_
                                        : std::numeric_limits<double>::infinity();
  const double sign = f0 > 0.0 ? -1.0 : 1.0;
  for (double d = bracket_tol_at(x0); d <= 0.1; d *= 4.0) {
    for (double side : {sign, -sign}) {
      const double x = std::min(x0 + side * d, hi);
      if (x == x0) {
        continue;
      }
      const double fx = target(x);
      if (std::isfinite(fx) && fx * f0 <= 0.0) {
        const std::pair<double, double> r =
            side < 0.0
                ? util::uniroot_smooth_bracket(target, x, x0, fx, f0, bound_root_tol,
                                               static_cast<size_t>(ci_niter))
                : util::uniroot_smooth_bracket(target, x0, x, f0, fx, bound_root_tol,
                                               static_cast<size_t>(ci_niter));
        // The continuity root from its feasible side, for the reason
        // find_root_psi takes it; the wet bound from the middle of its bracket.
        return find_root_crit == 1 ? r.first : 0.5 * (r.first + r.second);
      }
    }
  }
  // No sign change near x0: keep the bracket's answer rather than guess.
  return x0;
}

// The leaf potential that carries the soil's uptake at collar `psi_root`.
inline double Leaf::find_psi_leaf_from_psi_root(double psi_root, const std::vector<double>& psi_soil){
  E_from_Soil_to_Root_Collar(psi_root, psi_soil);

  double psi_leaf = transpiration_to_psi_leaf(E_up_, psi_root);
  return psi_leaf;
}

// find_wet_bound: the collar at which soil uptake equals the residual flux the
// leaf realises there, the wet bound and the shut leaf's collar. The target
// rises with the collar, and can change sign past the dry end, so the bracket
// is widened until it does, which it must: the residual flux is capped below
// the roots' ceiling. At g_min = 0 it is the zero-uptake collar, bit-identical
// (hazard 18). `E_up_at_crit` is the uptake at the dry end, NaN if not in hand.
inline double Leaf::find_wet_bound(double wettest_soil_layer,
                                   double E_up_at_wettest,
                                   double E_up_at_crit) {
  const std::vector<double>& soil = supply_psi_soil();
  auto target = [&](double x) -> double {
    E_from_Soil_to_Root_Collar(x, soil);
    return E_up_ - residual_flux(x, E_up_);
  };
  const double f_lo0 =
      E_up_at_wettest - residual_flux(wettest_soil_layer, E_up_at_wettest);
  double lo = wettest_soil_layer, f_lo = f_lo0;
  // The first bracket is [wettest, psi_crit] whatever this solve's dry end, so
  // the wet bound is the same root from the same bracket on every route to it
  // and a replay re-takes it bit for bit.
  double hi = std::max(psi_crit, wettest_soil_layer +
                                     bracket_tol_at(wettest_soil_layer));
  double f_hi = (hi == psi_crit && dry_end_ == psi_crit && std::isfinite(E_up_at_crit))
                    ? E_up_at_crit - residual_flux(psi_crit, E_up_at_crit)
                    : target(hi);
  if (!widen_to_sign_change(target, wettest_soil_layer, lo, f_lo, hi, f_hi)) {
    util::stop_infeasible("collar_bracket",
                          "find_wet_bound: the soil delivers less than the "
                          "residual flux at every collar potential tried; "
                          "wettest=" + util::to_string(wettest_soil_layer) +
                          "; last=" + util::to_string(hi));
  }
  try {
    return util::uniroot_smooth(target, lo, hi, f_lo, f_hi, bracket_tol,
                                ci_niter);
  } catch (const std::exception& e) {
    util::stop_infeasible("collar_bracket",
                          "find_wet_bound failed: " + std::string(e.what()) +
                          "; min=" + util::to_string(lo) +
                          "; max=" + util::to_string(hi));
  }
}

// The shut leaf's wet bound, bracketed from below by supply_wet_floor(), the
// lowest gravity balance among conducting layers: the target is <= 0 there and
// rises above, so one root-find to collar precision places it. Where no layer
// conducts the floor is returned, the bound's limit from the wet side; the soil
// potential instead would leave a step of the gravity head in psi_leaf. The
// open leaf keeps find_wet_bound, whose bracket keeps g_min = 0 bit-identical.
inline double Leaf::find_shut_wet_bound() {
  const std::vector<double>& soil = supply_psi_soil();
  auto target = [&](double x) -> double {
    E_from_Soil_to_Root_Collar(x, soil);
    return E_up_ - residual_flux(x, E_up_);
  };
  const double floor = supply_wet_floor();
  double lo = floor, f_lo = target(lo);
  if (!(f_lo < 0.0)) {
    return lo;  // nothing conducts, or the floor is the bound itself
  }
  double hi = std::max(supply_wet_ceiling(), lo + bracket_tol_at(lo));
  double f_hi = target(hi);
  if (!widen_to_sign_change(target, floor, lo, f_lo, hi, f_hi)) {
    util::stop_infeasible("collar_bracket",
                          "find_shut_wet_bound: the soil delivers less than "
                          "the residual flux at every collar tried; floor=" +
                          util::to_string(floor) + "; last=" + util::to_string(hi));
  }
  if (f_hi == 0.0) {
    return hi;
  }
  try {
    return util::uniroot_smooth(target, lo, hi, f_lo, f_hi, bound_root_tol,
                                static_cast<size_t>(ci_niter));
  } catch (const std::exception& e) {
    util::stop_infeasible("collar_bracket",
                          "find_shut_wet_bound failed: " + std::string(e.what()) +
                          "; min=" + util::to_string(lo) +
                          "; max=" + util::to_string(hi));
  }
}

// The leaf potential carrying exactly the residual flux from `psi_upstream`:
// the wet end of a leaf route. It exists because the residual flux is strictly
// below the stem's ceiling, and may lie past psi_crit. Solved to 1e-15 because
// a leaf that declines to open sits exactly here.
inline double Leaf::leaf_wet_bound(double psi_upstream) {
  // A leaf route has no roots upstream of psi_soil, so it clears the roots'
  // ceiling a previous collar solve seated.
  residual_root_ceiling_ = std::numeric_limits<double>::infinity();
  if (g_min == 0.0) {
    return psi_upstream;
  }
  auto target = [&](double psi) -> double {
    const double E = transpiration(psi, psi_upstream);
    return E - residual_flux(psi_upstream, E);
  };
  double lo = psi_upstream, f_lo = -residual_flux(psi_upstream, 0.0);
  double hi = std::max(psi_crit, psi_upstream + bracket_tol_at(psi_upstream));
  double f_hi = target(hi);
  if (!widen_to_sign_change(target, psi_upstream, lo, f_lo, hi, f_hi)) {
    util::stop_infeasible("leaf_wet_bound",
                          "the stem cannot carry the residual flux from "
                          "psi_upstream=" + util::to_string(psi_upstream));
  }
  return util::uniroot_smooth(target, lo, hi, f_lo, f_hi, bound_root_tol,
                              static_cast<size_t>(ci_niter));
}

// The shut operating point every early exit of the collar solve places
// (hazard 18): the collar where uptake meets the residual flux, the stem where
// it carries it, and every output from the ordinary evaluator there, so profit
// is continuous across the exit. Tagged here, not at the call sites, so no
// exit can shut the leaf unclassified.
template <Leaf::CostCurve K>
inline void Leaf::set_shut_state() {
  const std::vector<double>& soil = supply_psi_soil();
  const double collar = find_shut_wet_bound();
  wet_bound_root_ = collar;
  opt_root_psi_ = collar;
  E_from_Soil_to_Root_Collar(collar, soil);
  // The stem carries the uptake. Where that is zero (g_min = 0) or the bound's
  // root-find leaves it at or below zero (a ~1e-50 residual flux at the end of
  // the stem curve), the stem sits at the collar, the flux's limit from above:
  // inverting a negative flux would put the stem wetter than the collar.
  opt_psi_leaf_ = (residual_flux(collar, E_up_) == 0.0 || !(E_up_ > 0.0))
                      ? collar
                      : find_psi_leaf_from_psi_root(collar, soil);
  if constexpr (Curve<K>::prepares_normalisers) {
    if (!profitmax_has_objective()) {
      set_profitmax_degenerate_state(opt_psi_leaf_, collar);
      operating_point_kind_ = OperatingPointKind::HydraulicShutdown;
      return;
    }
  }
  profit_ = profit_psi_leaf_for<K>(opt_psi_leaf_, collar);
  // NA where the supply's conductance is unavailable (a one-layer wet bound is
  // its gravity-balance kink).
  lambda_emergent_ = lambda_collar_for<K>();
  operating_point_kind_ = OperatingPointKind::HydraulicShutdown;
  if (!std::isfinite(profit_)) {
    util::stop_infeasible("collar_solve",
                          "profit is not finite at the shut state; collar=" +
                              util::to_string(collar) +
                              "; psi_leaf=" + util::to_string(opt_psi_leaf_));
  }
}

// Setup and feasibility for every collar solve. Returns false when the leaf is
// shut, with its operating point already placed; otherwise true with
// [bound_a, bound_b] the feasible collar interval. `first_bracket` lets a curve
// not bounded at psi_crit try psi_crit as a provisional dry end.
template <Leaf::CostCurve K>
inline bool Leaf::prepare_collar_solve(double& bound_a, double& bound_b,
                                       bool first_bracket){

  // ⚠️ Cleared first, so a path that declines to classify reports Unsolved
  // rather than the previous plant's kind (hazard 8).
  operating_point_kind_ = OperatingPointKind::Unsolved;
  wet_bound_root_ = std::numeric_limits<double>::quiet_NaN();
  dry_bound_root_ = std::numeric_limits<double>::quiet_NaN();

  // Builds the supply's per-solve caches; returns the wettest rooted layer, the
  // smallest suction and the lower bracket end.
  const double wettest_soil_layer = supply_begin_solve();
  // Seated before anything reads the residual flux. Only the layered supply
  // bounds its uptake; a series resistance delivers without limit.
  residual_root_ceiling_ =
      (g_min > 0.0 && supply_kind_ == SupplyKind::MultiLayer)
          ? roots_.uptake_limit<double>(roots_.held_supply())
          : std::numeric_limits<double>::infinity();
  const double E_up_at_wettest =
      E_column_zero(wettest_soil_layer, supply_psi_soil());

  bracket_wettest_ = wettest_soil_layer;
  bracket_E_up_at_wettest_ = E_up_at_wettest;

  // One bracket with the stem held at `end` (seated in dry_end_ first, which the
  // helpers read): false if no collar short of it both supplies the cuticle and
  // leaves the stem flow, else [wet bound, continuity root]. The two root-finds
  // share endpoint values, handed over rather than re-taken: at the end the
  // column term is exactly zero, and at the wettest layer both read one uptake.
  auto bracket = [&](double end) -> bool {
    dry_end_ = end;
    if (wettest_soil_layer >= end) {
      return false;
    }
    // The most the soil delivers with the stem short of `end`.
    const double E_at_end = E_column(end, supply_psi_soil(), end);
    if (E_at_end < residual_flux(end, E_at_end)) {
      return false;
    }
    const double collar_at_dry_end = find_root_psi(
        wettest_soil_layer, supply_psi_soil(), 1,
        E_up_at_wettest - transpiration(end, wettest_soil_layer), E_at_end, end);
    if (collar_at_dry_end >= end) {
      return false;
    }
    bound_a = find_wet_bound(wettest_soil_layer, E_up_at_wettest, E_at_end);
    // The root's own critical potential bounds nothing: a clamp to it puts a
    // kink in the stem and root traits wherever the two tie.
    bound_b = collar_at_dry_end;
    // Ends of an interval narrower than bracket_tol cannot be trusted, so they
    // are solved to collar precision before anything is decided on them.
    if (bound_b - bound_a < bracket_tol_at(bound_b)) {
      bound_a = polish_bracket_bound(bound_a, supply_psi_soil(), 0);
      bound_b = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
    }
    return bound_b - bound_a > collar_root_tol;
  };

  // For a curve not bounded at psi_crit, psi_crit is only a first bracket: it
  // holds the optimum wherever profit already falls there, and
  // maximise_profit_over_collar extends it otherwise. Where it holds no interval
  // the end of the stem curve is tried before the leaf is called shut.
  if constexpr (!cost_domain_ends_at_psi_crit<K>()) {
    if (first_bracket && bracket(psi_crit)) {
      dry_end_provisional_ = true;
      return true;
    }
  }
  dry_end_provisional_ = false;
  if (!bracket(dry_end<K>())) {
    set_shut_state<K>();
    return false;
  }
  return true;
}

template <Leaf::CostCurve K>
inline bool Leaf::seat_leaf_at_dry_pin() {
  if constexpr (cost_domain_ends_at_psi_crit<K>()) {
    return false;
  } else {
    if (operating_point_kind_ != OperatingPointKind::BoundaryDry) {
      return false;
    }
    opt_psi_leaf_ = vulnerability_psi_end(stem_b, stem_c);
    return true;
  }
}

inline double Leaf::extend_collar_bracket(double bound_a) {
  const double end = vulnerability_psi_end(stem_b, stem_c);
  dry_end_ = end;
  dry_end_provisional_ = false;
  const double E_at_end = E_column(end, supply_psi_soil(), end);
  const double collar_at_dry_end = find_root_psi(
      bracket_wettest_, supply_psi_soil(), 1,
      bracket_E_up_at_wettest_ - transpiration(end, bracket_wettest_), E_at_end,
      end);
  // A narrow interval's dry end is polished, as in prepare_collar_solve, or a
  // root within bracket_tol_at of the wet bound can land on it and report shut
  // a leaf whose profit still rises. The wet end is not: moving it would strand
  // the caller's `lo`, already stepped in from it.
  const double bound_b = std::max(collar_at_dry_end, bound_a);
  if (bound_b - bound_a < bracket_tol_at(bound_b)) {
    return std::max(polish_bracket_bound(bound_b, supply_psi_soil(), 1),
                    bound_a);
  }
  return bound_b;
}

// The profit-maximising collar on [bound_a, bound_b], by a bracketing
// root-find of dprofit/dpsi == 0; prepare_collar_solve must have seated the
// soil caches. The first-order condition, not a search of profit: a search's
// argmax wanders within its tolerance as parameters move, an offset that
// swamps, and can reverse, the hydraulic traits' true response (hazard 3).
// It is util::maximise_over_closed_interval_foc at n = 0 plus the
// OperatingPointKind every exit sets, which pinned golden rows read: the one
// solver plus a classification (.claude/CLAUDE.md "One solver").
template <Leaf::CostCurve K>
inline double Leaf::maximise_profit_over_collar(double bound_a, double bound_b) {
  double width = bound_b - bound_a;
  const auto dprofit = [&](double psi, bool* feasible = nullptr) {
    return dprofit_at_collar_psi<K>(psi, feasible);
  };

  // Places lo at the wet bound: the answer where profit falls from the wet end,
  // and for a leaf with no objective to maximise (below).
  double lo = bound_a;
  double hi = bound_b;
  // Every exit returns through this, so none returns a collar unclassified.
  const auto classified = [this](OperatingPointKind kind, double collar) {
    operating_point_kind_ = kind;
    return collar;
  };
  const auto place_at_wet_bound = [&]() {
    const double a = polish_bracket_bound(bound_a, supply_psi_soil(), 0);
    wet_bound_root_ = a;
    if (g_min > 0.0) {
      // With a residual conductance the bound is a regular point (no stomatal
      // flux, the cuticle conducts), and a leaf that declines to open sits on it.
      lo = a;
    } else if (a < hi) {
      // With nothing conducting at the bound a dark leaf has no finite ci there
      // (ci = ca + R_d/g_c), so the state is taken a step inside, always from
      // the polished root and on the bound's own scale: a fraction of the width
      // moves with the dry bound, a dependence the rows do not have. Not below
      // 1e-6: the rows' error grows with the step, but at 1e-7 the stem
      // inversion reaches its rounding floor. Only the no-flow test is repeated:
      // f_lo's sign cannot flip across the ~1e-4 the bound moved.
      const double scale = std::max(a, 1e-3);
      const double limit = 0.5 * (bound_b - a);
      for (double frac = 1e-6; frac * scale < limit; frac *= 10.0) {
        const double x = a + frac * scale;
        const double s = find_psi_leaf_from_psi_root(x, supply_psi_soil());
        if (std::isfinite(s) && x < s) {
          lo = x;
          break;
        }
      }
    }
  };

  // A leaf with no objective is shut: a log-link one with A <= 0 at the dry end
  // (objective_has_log_link), or ProfitMax without its normalisers. Decided
  // before either end's gradient is read: through the log link at a negative A
  // they say nothing.
  bool shut = false;
  if constexpr (objective_has_log_link<K>()) {
    shut = dark_under_log_link<K>(
        find_psi_leaf_from_psi_root(bound_b, supply_psi_soil()), bound_b);
  } else if constexpr (Curve<K>::prepares_normalisers) {
    shut = !profitmax_has_objective();
  }
  if (shut) {
    // Where the open path would start: the bound at g_min > 0, a step inside
    // the no-flow sentinel at g_min = 0, should placement find no step.
    lo = g_min > 0.0 ? bound_a : bound_a + 1e-6 * width;
    place_at_wet_bound();
    return classified(OperatingPointKind::BoundarySoil, lo);
  }

  // The WET endpoint. With g_min > 0 it is an ordinary point, tried as-is. At
  // g_min = 0 it is the zero-uptake collar, where dprofit's 0.0 is a sentinel a
  // bracketing solver would report as the optimum, so the read steps inside
  // until the gradient is real (the infeasible region is under 1e-6 of the
  // bracket, so the first step lands). Starting at the bound there instead
  // moves every interior row by the solver floor.
  bool ok_lo = false;
  double f_lo = 0.0;
  for (double frac = (g_min > 0.0 ? 0.0 : 1e-6); frac < 0.5;
       frac = (frac == 0.0 ? 1e-6 : frac * 10.0)) {
    lo = bound_a + frac * width;
    f_lo = dprofit(lo, &ok_lo);
    if (ok_lo && std::isfinite(f_lo)) {
      break;
    }
  }

  // The DRY endpoint. For a curve bounded at psi_crit it is tried as-is. For any
  // other it is the continuity root with the stem at the end of its curve,
  // placed only to bracket_tol_at, and a read there can fall past the root,
  // where transpiration() refuses; so the first read is a step inside.
  // Profit can still rise at that end, which is the model: where the stem has
  // little left to carry and A climbs steeply in E, the leaf transpires all the
  // stem can carry, a pin at the dry end. The probes only say where to look;
  // the pin is decided at the end itself, solved to collar precision (hazard 19).
  // On a first bracket (stem at psi_crit) the end is an ordinary point: profit
  // falling there keeps the bracket, anything else extends it.
  constexpr bool bounded = cost_domain_ends_at_psi_crit<K>();
  bool ok_hi = false;
  bool rising = false;
  double f_hi = 0.0;
  bool settled = false;
  if (dry_end_provisional_) {
    f_hi = dprofit(hi, &ok_hi);
    settled = ok_hi && std::isfinite(f_hi) && f_hi < 0.0;
    if (!settled) {
      bound_b = extend_collar_bracket(bound_a);
      width = bound_b - bound_a;
    }
  }
  // The bracket the interior root-find runs on: the step inside the wet bound,
  // and `hi` wherever the probes leave it.
  double root_lo = lo;
  double f_root_lo = f_lo;
  // Whether the dry end has decided the point: an interior maximum between
  // the first probe and the end, or the pin. Neither needs the probes below.
  bool decided = false;
  if (!settled) {
    hi = bounded ? bound_b : bound_b - 1e-6 * width;
    f_hi = dprofit(hi, &ok_hi);
  }
  // ⚠️ Where the first probe reads a rising profit the end decides (hazard 19):
  // the marginal falls steeply over the last micro-MPa, so a probe's sign is
  // not the end's. Solved to collar precision, falling at the end is an
  // interior maximum between probe and end, rising is the pin; deciding on the
  // probes made the solve jump. Read before the deeper probes, each a ci solve.
  // Only where profit rises from the wet end too: falling from it is two
  // maxima, which the refusal below compares.
  if (!bounded && !settled && f_lo > 0.0 && ok_hi && std::isfinite(f_hi) &&
      f_hi >= 0.0) {
    dry_bound_root_ = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
    bool ok_end = false;
    const double f_end = dprofit(dry_bound_root_, &ok_end);
    if (ok_end && std::isfinite(f_end) && f_end < 0.0 && dry_bound_root_ > hi) {
      root_lo = hi;
      f_root_lo = f_hi;
      hi = dry_bound_root_;
      f_hi = f_end;
      decided = true;
    } else if (ok_end && std::isfinite(f_end) && f_end >= 0.0) {
      rising = true;
      decided = true;
    }
    bound_b = dry_bound_root_;
  }
  if (!settled && !decided &&
      !(ok_hi && std::isfinite(f_hi) && (bounded || f_hi < 0.0))) {
    rising = ok_hi && std::isfinite(f_hi);
    for (double frac = bounded ? 1e-6 : 1e-5; frac < 0.5; frac *= 10.0) {
      hi = bound_b - frac * width;
      f_hi = dprofit(hi, &ok_hi);
      if (ok_hi && std::isfinite(f_hi) && (bounded || f_hi < 0.0)) {
        rising = false;
        break;
      }
      rising = rising || (ok_hi && std::isfinite(f_hi));
    }
  }
  if (!bounded && rising) {
    // Rising at the end, or at every probe where the end was not read: the pin,
    // on the end solved to collar precision.
    if (!std::isfinite(dry_bound_root_)) {
      dry_bound_root_ = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
      bound_b = dry_bound_root_;
    }
    hi = dry_bound_root_;
    f_hi = 1.0;
    ok_hi = true;
  }

  // No usable gradient at one end, which no golden-grid row reaches: a signal,
  // answered by golden section rather than a failure. `!ok` is the supply
  // reporting no informative gradient (SolverRefused); `ok && !isfinite` is a
  // partial feasibility said would be finite (NonFiniteGradient). A NaN would
  // otherwise fail both pin tests below and reach the root-find unbracketed.
  if (!ok_lo || !ok_hi || !std::isfinite(f_lo) || !std::isfinite(f_hi)) {
    const OperatingPointKind kind =
        ((ok_lo && !std::isfinite(f_lo)) || (ok_hi && !std::isfinite(f_hi)))
            ? OperatingPointKind::NonFiniteGradient
            : OperatingPointKind::SolverRefused;
    return classified(kind, util::golden_section_max(
                                [&](double bound) {
                                  const double psi_leaf = find_psi_leaf_from_psi_root(
                                      bound, supply_psi_soil());
                                  return profit_psi_leaf_for<K>(psi_leaf, bound);
                                },
                                bound_a, bound_b, GSS_tol_abs));
  }

  // A CONSTRAINED optimum: profit still climbs at one end, so the answer is that
  // end, which is common on dry soil. Return lo/hi, never the raw bound_a/bound_b:
  // at g_min = 0 the wet bound is the sentinel stepped over above. Where a row
  // crosses between interior and pinned the argmax has a real kink in traits.
  // A pin's bound is the answer, so only a side that would pin is polished
  // (polish_bracket_bound says why); an interior solve pays nothing.
  // dprofit <= 0 at the wet end and >= 0 at the dry end is not a third pin: the
  // interior stationary point is a minimum, so it is SolverRefused and the
  // larger of the two ends' profits is returned.
  if (f_lo <= 0.0) {
    place_at_wet_bound();
  }
  if (f_hi >= 0.0 && hi == bound_b && !std::isfinite(dry_bound_root_)) {
    const double b = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
    // f_hi's sign cannot flip across the ~1e-4 the bound moved.
    if (b != bound_b && b > lo) {
      hi = b;
    }
    dry_bound_root_ = hi;
  }

  if (f_lo <= 0.0 && f_hi >= 0.0) {
    const double p_lo = profit_psi_leaf_for<K>(
        find_psi_leaf_from_psi_root(lo, supply_psi_soil()), lo);
    const double p_hi = profit_psi_leaf_for<K>(
        find_psi_leaf_from_psi_root(hi, supply_psi_soil()), hi);
    const double larger = !std::isfinite(p_hi)   ? lo
                          : !std::isfinite(p_lo) ? hi
                          : p_hi > p_lo          ? hi
                                                 : lo;
    return classified(OperatingPointKind::SolverRefused, larger);
  }
  if (f_lo <= 0.0) {
    return classified(OperatingPointKind::BoundarySoil, lo);
  }
  if (f_hi >= 0.0) {
    return classified(OperatingPointKind::BoundaryDry, hi);
  }

  // Interior stationary point: f_lo > 0 > f_hi, a valid bracket.
  double root;
  if (!bounded && decided && std::isfinite(dry_bound_root_)) {
    // Near the end the marginal is close to linear in the log of the distance
    // to it (dry_end_neighbourhood), so the root is found in that coordinate: a
    // few superlinear steps where the potential takes a dozen ci solves. A root
    // within dry_pin_resolution of the end is the pin.
    const double end = dry_bound_root_;
    bool ok_in = false;
    const double f_in = end - root_lo > dry_pin_resolution
                            ? dprofit(end - dry_pin_resolution, &ok_in)
                            : 0.0;
    if (!(ok_in && std::isfinite(f_in) && f_in < 0.0)) {
      return classified(OperatingPointKind::BoundaryDry, end);
    }
    const auto in_log_distance = [&](double s) {
      return dprofit(end - std::exp(s));
    };
    // 1e-10 in the logarithm is 1e-12 MPa at the neighbourhood's outer edge.
    const double s = util::uniroot_smooth(
        in_log_distance, std::log(dry_pin_resolution), std::log(end - root_lo),
        f_in, f_root_lo, 1e-10, static_cast<size_t>(ci_niter));
    root = end - std::exp(s);
  } else {
    root = util::uniroot_smooth(dprofit, root_lo, hi, f_root_lo, f_hi,
                                collar_root_tol, static_cast<size_t>(ci_niter));
  }
  if constexpr (!bounded) {
    // A root within dry_end_neighbourhood of the dry end reads that end for its
    // curvature, so the end is polished; within dry_pin_resolution it is the
    // pin. Only where the bracket's end is the domain's: a first bracket's end
    // (stem at psi_crit) is an ordinary point, and a replay, bracketing the
    // whole domain, would take a different stencil.
    if (!dry_end_provisional_ &&
        bound_b - root < dry_end_neighbourhood + bracket_tol_at(bound_b)) {
      if (!std::isfinite(dry_bound_root_)) {
        dry_bound_root_ = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
      }
      if (dry_bound_root_ - root < dry_pin_resolution) {
        return classified(OperatingPointKind::BoundaryDry, dry_bound_root_);
      }
    }
  }
  return classified(OperatingPointKind::Interior, root);
}

template <Leaf::CostCurve K>
inline void Leaf::begin_collar_point() {
  // Without the check CF77 with an unset lambda maximises a NaN objective and
  // returns a bracket end as an operating point. The normalisers' reference is
  // psi_soil where there is one and the zero-flow collar where there is not
  // (profitmax_reference_psi); unseeded, a normalised profit is not finite.
  check_cost_parameters<K>();
  if constexpr (Curve<K>::prepares_normalisers) {
    prepare_profitmax();
  }
}

template <Leaf::CostCurve K>
inline void Leaf::find_root_collar_psi_for(){
  begin_collar_point<K>();
    double bound_a, bound_b;
    if (!prepare_collar_solve<K>(bound_a, bound_b, true)) {
      return;
    }

    const double opt_root_psi = maximise_profit_over_collar<K>(bound_a, bound_b);

    opt_psi_leaf_ = find_psi_leaf_from_psi_root(opt_root_psi, supply_psi_soil());
    seat_leaf_at_dry_pin<K>();

    opt_root_psi_ = opt_root_psi;
    // The leaf route's degenerate exit: its outputs, at the wet bound placed.
    if constexpr (Curve<K>::prepares_normalisers) {
      if (!profitmax_has_objective()) {
        set_profitmax_degenerate_state(opt_psi_leaf_, opt_root_psi);
        return;
      }
    }
    profit_ = profit_psi_leaf_for<K>(opt_psi_leaf_, opt_root_psi);
    // On the single-layer solvers' axis: the collar is the free variable, so the
    // marginal cost carries the soil-to-collar series-resistance correction.
    lambda_emergent_ = lambda_collar_for<K>();

    if(!std::isfinite(profit_)){
        util::stop_infeasible("collar_solve", "non-finite profit; opt_psi_leaf_=" + util::to_string(opt_psi_leaf_) +
             "; opt_root_psi_=" + util::to_string(opt_root_psi_) +
             "; bound_a=" + util::to_string(bound_a) +
             "; bound_b=" + util::to_string(bound_b) +
             "; E_up_=" + util::to_string(E_up_) +
             "; assim_colimited_=" + util::to_string(assim_colimited_) +
             "; hydraulic_cost_=" + util::to_string(hydraulic_cost_));
    }
}

// The collar solve on the seated curve, whatever route and method are seated;
// optimise() solves the seated model as seated. It writes the fallback flag as
// optimise() does, since a caller may read it after either (hazard 8).
inline void Leaf::find_root_collar_psi() {
  last_solve_fell_back_ = false;
  with_curve(cost_curve_, [&](auto tag) { find_root_collar_psi_for<tag.value>(); });
}

// The operating point at a prescribed collar, on the seated curve. A collar
// outside the feasible interval is refused as infeasible (profit_at_collar_psi),
// so a caller tracking the collar as a state can reject the step.
inline double Leaf::evaluate_root_collar_psi(double target_opt_root_psi) {
  return with_curve(cost_curve_, [&](auto tag) {
    return evaluate_root_collar_psi_for<tag.value>(target_opt_root_psi);
  });
}

inline double Leaf::dprofit_droot_collar_psi(double opt_root_psi, bool* feasible) {
  return with_curve(cost_curve_, [&](auto tag) {
    return dprofit_droot_collar_psi_for<tag.value>(opt_root_psi, feasible);
  });
}

template <Leaf::CostCurve K>
inline double Leaf::evaluate_root_collar_psi_for(double target_opt_root_psi){
    begin_collar_point<K>();
    double bound_a, bound_b;
    if (!prepare_collar_solve<K>(bound_a, bound_b)) {
      // Shut: profit_ is already set.
      return profit_;
    }
    // The bracket's dry end is placed only to its tolerance, and an interior
    // maximum can sit between it and the continuity root: a target past the
    // bracket is read against the root, solved to collar precision.
    if (target_opt_root_psi > bound_b) {
      const double b = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
      if (b > bound_b) {
        dry_bound_root_ = b;
        bound_b = b;
      }
    }

    return profit_at_collar_psi<K>(target_opt_root_psi, bound_a, bound_b);
}

template <Leaf::CostCurve K>
inline double Leaf::replay_operating_point_for(double collar,
                                               OperatingPointKind kind) {
    // Seeded as the solve seeds them (hazard 8): on a leaf that solved another
    // point since, a normalised curve's profit would read that point's.
    begin_collar_point<K>();
    double bound_a, bound_b;
    double out;
    if (!prepare_collar_solve<K>(bound_a, bound_b)) {
      out = profit_;
    } else {
      if (kind == OperatingPointKind::BoundarySoil) {
        wet_bound_root_ = polish_bracket_bound(bound_a, supply_psi_soil(), 0);
      }
      // A dry pin sits on the polished end, and an interior point near it reads
      // that end for its curvature: seat it as the solve did, under the solve's
      // own test, bracket tolerance included, so the two take the same stencil.
      if (kind == OperatingPointKind::BoundaryDry ||
          (kind == OperatingPointKind::Interior &&
           bound_b - collar < dry_end_neighbourhood + bracket_tol_at(bound_b))) {
        dry_bound_root_ = polish_bracket_bound(bound_b, supply_psi_soil(), 1);
      }
      // A pinned solve's collar is a polished bound, up to ~1e-4 outside the
      // bracket returned here; clamped back it would replay a different point,
      // so the interval is widened to hold it.
      bool evaluate = true;
      if constexpr (Curve<K>::prepares_normalisers) {
        // The solve's degenerate exit, placed without evaluating an objective
        // whose normalisers it cannot divide by.
        if (!profitmax_has_objective()) {
          opt_psi_leaf_ = find_psi_leaf_from_psi_root(collar, supply_psi_soil());
          opt_root_psi_ = collar;
          set_profitmax_degenerate_state(opt_psi_leaf_, collar);
          out = profit_;
          evaluate = false;
        }
      }
      if (evaluate) {
        out = profit_at_collar_psi<K>(collar, std::min(bound_a, collar),
                                      std::max(bound_b, collar));
      }
    }
    operating_point_kind_ = kind;
    if (seat_leaf_at_dry_pin<K>()) {
      profit_ = profit_psi_leaf_for<K>(opt_psi_leaf_, opt_root_psi_);
      lambda_emergent_ = lambda_collar_for<K>();
      out = profit_;
    }
    return out;
}

// The operating point at a prescribed collar, after prepare_collar_solve: a
// caller evaluating several collars in one step reuses the soil-side caches.
// A collar outside [bound_a, bound_b] is refused, not projected: projected, the
// outputs follow the bound rather than the collar asked for. Infeasible-class,
// so a consumer integrating the collar can reject the step, not end the run.
template <Leaf::CostCurve K>
inline double Leaf::profit_at_collar_psi(double target_opt_root_psi,
                                  double bound_a, double bound_b){
    if (!(target_opt_root_psi >= bound_a && target_opt_root_psi <= bound_b)) {
      util::stop_infeasible(
          "prescribed_potential",
          "the prescribed collar potential " +
              util::format_double(target_opt_root_psi) +
              " MPa is outside the feasible collar interval [" +
              util::format_double(bound_a) + ", " +
              util::format_double(bound_b) + "] MPa");
    }
    const double opt_root_psi = target_opt_root_psi;

    opt_psi_leaf_ = find_psi_leaf_from_psi_root(opt_root_psi, supply_psi_soil());
    opt_root_psi_ = opt_root_psi;
    profit_ = profit_psi_leaf_for<K>(opt_psi_leaf_, opt_root_psi);
    // dC/dE is a property of the point, defined off the optimum, where it is not dA/dE.
    lambda_emergent_ = lambda_collar_for<K>();
    // Here a non-zero dprofit is the answer (TF24f's acclimation rate), not a
    // residual to check against zero.
    operating_point_kind_ = OperatingPointKind::Prescribed;

    if(!std::isfinite(profit_)){
        util::stop_infeasible("collar_solve",
             "non-finite profit in evaluate_root_collar_psi; "
             "target=" + util::to_string(target_opt_root_psi) +
             "; opt_root_psi=" + util::to_string(opt_root_psi) +
             "; bound_a=" + util::to_string(bound_a) +
             "; bound_b=" + util::to_string(bound_b));
    }
    return profit_;
}

// Exact d(profit)/d(collar): A'(ci) dci/dpsi - C'(psi_leaf) dpsi_leaf/dpsi,
// with dci/dpsi from the implicit-function theorem on the ci residual and A', C'
// by forward AD of the model's own kernels; vignette("the-models"), "The
// benefit, and the root-find inside it". This wrapper seats the supply caches on
// psi_soil_; the collar solve calls dprofit_at_collar_psi, seated once.
template <Leaf::CostCurve K>
inline double Leaf::dprofit_droot_collar_psi_for(double opt_root_psi, bool* feasible) {
  supply_begin_solve();
  return dprofit_at_collar_psi<K>(opt_root_psi, feasible);
}

inline std::vector<double> Leaf::dprofit_droot_collar_psi_checked(
    double opt_root_psi) {
  bool feasible = false;
  const double d = dprofit_droot_collar_psi(opt_root_psi, &feasible);
  return {d, feasible ? 1.0 : 0.0};
}

template <Leaf::CostCurve K>
inline double Leaf::dprofit_at_collar_psi(double opt_root_psi, bool* feasible) {
  const double psi = opt_root_psi;
  // Infeasible until the exits below are passed: their 0.0 must be
  // distinguishable from a stationary point.
  if (feasible != nullptr) {
    *feasible = false;
  }

  // Past the continuity root with the stem at the end of its curve the soil
  // delivers more than the stem can carry: no operating point, reported like the
  // other no-gradient exits rather than thrown (a difference across a dry pin reaches it).
  // The uptake's collar slope comes from the same walk (dpsi_leaf/dpsi below).
  const double dEup_dpsi = E_and_slope_from_soil(psi, supply_psi_soil());
  if (!stem_carries(E_up_, psi)) {
    return 0.0;
  }
  const double psi_leaf = transpiration_to_psi_leaf(E_up_, psi);
  // psi >= psi_leaf is the no-flow / reversed case, caught before the ci solve:
  // there psi_leaf_to_ci throws (transpiration() refuses the pair) rather than
  // returning non-finite, and the throw would end a consumer's run.
  if (!std::isfinite(psi_leaf) || psi >= psi_leaf) {
    return 0.0;  // shut-down / infeasible: no informative gradient
  }

  // With the energy balance on, the temperature parameters are seated at this
  // candidate's Tleaf(E), as the objective is: otherwise this is the first-order
  // condition of the air-temperature model. Same order as
  // set_leaf_states_rates_from_psi_leaf: transpiration, temperature, ci solve.
  const double E = transpiration(psi_leaf, psi);
  double dT_dE = 0.0;
  double Tleaf_here = leaf_temp_;
  if (use_energy_balance_) {
    Tleaf_here = leaf_temp_from_E(E, &dT_dE);
    update_temperature_dependent_params(Tleaf_here);
    set_leaf_vpd(Tleaf_here);
  }

  // Below the wet bound the stomatal flux is zero, the corner the solve keeps
  // (E_s >= 0); with nothing conducting (the g_min = 0 wet bound) there is no
  // derivative to take.
  const double r_star = residual_flux(psi, E);
  const double E_s = E - r_star < 0.0 ? 0.0 : E - r_star;
  const double cuticle = cuticular_cond_CO2();
  if (E_s == 0.0 && cuticle == 0.0) {
    return 0.0;
  }

  const double ci = psi_leaf_to_ci(psi_leaf, psi);
  if (!std::isfinite(ci)) {
    return 0.0;
  }

  // The stem table's slope at either end, each read once: a read is an
  // out-of-line table lookup, and the compiler does not share two of them.
  const double f_psi = stem_curve_integral_deriv(psi);
  const double f_leaf = stem_curve_integral_deriv(psi_leaf);
  // dpsi_leaf/dpsi, shared by both branches below.
  double dpsileaf_dpsi;
  if (std::isfinite(dEup_dpsi)) {
    const double dEpsileaf_dpsi =
        dEup_dpsi / leaf_specific_conductance_max_ + f_psi;
    // dpsi_leaf/dE_psi_leaf = 1/G'(psi_leaf), at the psi_leaf inverted above
    // from the same E_up_ and psi.
    dpsileaf_dpsi = (1.0 / f_leaf) * dEpsileaf_dpsi;
  } else {
    // NaN from the supply means "use finite differences" (hazard 6).
    const double h = 1e-6;
    dpsileaf_dpsi =
        (find_psi_leaf_from_psi_root(psi + h, supply_psi_soil()) -
         find_psi_leaf_from_psi_root(psi - h, supply_psi_soil())) / (2.0 * h);
  }


  // Past both exits: whatever is returned below is a real derivative, so a zero
  // from here IS a stationary point.
  if (feasible != nullptr) {
    *feasible = true;
  }

  // A'(ci) by forward AD of the kernel assim_colimited() instantiates, so it is
  // the derivative of the function evaluated, not of a hand-kept mirror.
  const double A_prime = assim_slope(ci);
  const double C_prime = cost_deriv<K>(psi_leaf, psi);

  // gc and its partials, on the stem table's own slope. The residual flux reads
  // psi_leaf only through Tleaf (the energy-balance term) and the collar through
  // its ceiling (below); both are exact zeros at g_min = 0.
  const double gc_const =
      atm_kpa_ * kg_to_mol_h2o / vpd_leaf_ / H2O_CO2_stom_diff_ratio_;
  double gc = gc_const * E_s + cuticle;
  double dgc_dpsileaf =
      gc_const * leaf_specific_conductance_max_ * f_leaf;
  double dgc_dpsi =
      gc_const * leaf_specific_conductance_max_ * (-f_psi);
  const double drstar_dpsi = residual_flux_dpsi_upstream(psi, Tleaf_here);
  dgc_dpsi = dgc_dpsi - gc_const * drstar_dpsi;
  // The boundary layer's share of the CO2 path: the value and every partial
  // through the split's one slope. [[unlikely]] is measured: without it this
  // branch slows the gate-off collar solve, and moving it out of line does not help.
  double bl_slope = 1.0;
  if (use_energy_balance_) [[unlikely]] {
    const std::pair<double, double> split = boundary_layer_split_CO2(gc);
    gc = split.first;
    bl_slope = split.second;
    dgc_dpsileaf *= bl_slope;
    dgc_dpsi *= bl_slope;
  }

  // IFT on g(ci; psi_leaf, psi): dci/dp = -(dg/dp)/(dg/dci).
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const double g_ci = A_prime * umol_to_mol + gc * inv_atm;      // dg/dci
  const double dci_dpsileaf = -(-dgc_dpsileaf * (ca_ - ci) * inv_atm) / g_ci;
  const double dci_dpsi_expl = -(-dgc_dpsi * (ca_ - ci) * inv_atm) / g_ci;


  const double dci_dpsi = dci_dpsileaf * dpsileaf_dpsi + dci_dpsi_expl;
  // The identity arm keeps every difference objective's operation order, which
  // the golden files pin; only a non-identity link takes the general path.
  double base;
  if constexpr (benefit_link<K>() == BenefitLink::Identity) {
    base = A_prime * dci_dpsi - C_prime * dpsileaf_dpsi;
  } else {
    const double A = assim_colimited_kernel(ci);
    base = benefit_link_deriv<K>(A) * (A_prime * dci_dpsi) - C_prime * dpsileaf_dpsi;
    // A log link has no derivative where the leaf is a net source (A <= 0), an
    // ordinary point once the cuticle conducts CO2: reported as absent, never NaN.
    if (!std::isfinite(base)) {
      if (feasible != nullptr) {
        *feasible = false;
      }
      return 0.0;
    }
  }
  // The cost's direct dependence on the collar, appended so a curve that does not
  // read upstream keeps the arm above operation for operation.
  if constexpr (cost_reads_upstream<K>()) {
    base = base - cost_deriv_upstream<K>(psi_leaf, psi);
  }
  // Gated at the CALL SITE, not just inside the callee: the block is out of line
  // (deliberately, so adding it cannot change FMA contraction in this inlined
  // body), and an out-of-line call costs even when it returns 0.0 immediately.
  if (!use_energy_balance_) {
    return base;
  }
  return base + dprofit_energy_balance_term<K>(
                    ci, gc, g_ci, inv_atm, gc_const, A_prime, dgc_dpsileaf,
                    dgc_dpsi, dpsileaf_dpsi, dT_dE, Tleaf_here, E_s, r_star,
                    drstar_dpsi,
                    residual_cond_temp_deriv(psi, Tleaf_here, gc_const), bl_slope);
}

} // namespace phylloptim

#endif
