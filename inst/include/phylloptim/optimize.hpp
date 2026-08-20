// -*-c++-*-
#ifndef PHYLLOPTIM_OPTIMIZE_HPP_
#define PHYLLOPTIM_OPTIMIZE_HPP_

// Robust 1-D function minimiser: Brent's method (golden-section search with
// parabolic interpolation). This is the algorithm behind R's optimize() /
// Forsythe-Malcolm-Moler `fmin`. It converges super-linearly on a smooth
// objective near its optimum but falls back to a golden-section step whenever
// the parabolic fit is unreliable, so it keeps the bracketing robustness that
// the leaf hydraulic solvers depend on (cf. the TOMS748 caveat in uniroot.h:
// here we only ever evaluate STRICTLY INTERIOR points, never the clamped
// bracket endpoints).

#include <cfloat>  // DBL_EPSILON
#include <cmath>
#include <limits>
#include <vector>

namespace phylloptim {
namespace util {

// Minimise f over [ax, bx] (requires ax <= bx). Returns the argmin; on return
// fmin (if non-null) holds f at the argmin. `tol` is the absolute tolerance on
// the location of the minimum; values below ~sqrt(DBL_EPSILON)*|x| are not
// useful. A direct transcription of R's Brent_fmin (src/appl/fmin.c).
template <typename Function>
double brent_fmin(Function f, double ax, double bx, double tol,
                  double* fmin = nullptr) {
  // c is the squared inverse of the golden ratio (~0.3819660).
  const double c = (3.0 - std::sqrt(5.0)) * 0.5;
  const double eps = std::sqrt(DBL_EPSILON);

  double a = ax, b = bx;
  double v = a + c * (b - a);
  double w = v, x = v;
  double d = 0.0, e = 0.0;
  double fx = f(x);
  double fv = fx, fw = fx;
  const double tol3 = tol / 3.0;

  for (;;) {
    const double xm = (a + b) * 0.5;
    const double tol1 = eps * std::abs(x) + tol3;
    const double tol2 = tol1 * 2.0;

    // Convergence check.
    if (std::abs(x - xm) <= tol2 - (b - a) * 0.5)
      break;

    double p = 0.0, q = 0.0, r = 0.0;
    bool use_golden = true;

    if (std::abs(e) > tol1) {
      // Fit a parabola through (x,fx), (v,fv), (w,fw).
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.0;
      if (q > 0.0)
        p = -p;
      else
        q = -q;
      r = e;
      e = d;
      // Accept the parabolic step only if it is well inside (a,b) and shrinks.
      if (std::abs(p) < std::abs(0.5 * q * r) &&
          p > q * (a - x) && p < q * (b - x)) {
        d = p / q;
        const double u = x + d;
        // Keep the new point away from the bracket endpoints.
        if (u - a < tol2 || b - u < tol2)
          d = (x < xm) ? tol1 : -tol1;
        use_golden = false;
      }
    }

    if (use_golden) {
      e = (x < xm) ? (b - x) : (a - x);
      d = c * e;
    }

    // Evaluate f at a point at least tol1 away from x.
    double u;
    if (std::abs(d) >= tol1)
      u = x + d;
    else
      u = (d > 0.0) ? (x + tol1) : (x - tol1);
    const double fu = f(u);

    // Update the bracket and the three best points.
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }

  if (fmin != nullptr)
    *fmin = fx;
  return x;
}

// Golden-section search for the MAXIMUM of a unimodal f over [ax, bx]. Returns
// the argmax (midpoint of the final bracket); terminates when the bracket width
// falls to `tol`. Reuses one interior golden point per iteration, so it costs a
// single new f() evaluation per step after the initial two.
//
// Why this exists alongside brent_fmin: brent_fmin converges faster but its
// parabolic step makes the argmax a *non-smooth* function of the inputs. Where
// the argmax does not feed a gradient (the single-layer leaf optimisers), prefer
// brent_fmin, which was measured ~2.3-2.6x faster there.
//
// ⚠️ **The production collar solver no longer uses this** -- PLAN 11a replaced it
// with a safeguarded root-find on the first-order condition
// (Leaf::maximise_profit_over_collar), and this is now only that solver's fallback
// for the case where neither bracket endpoint has a usable gradient. The reasoning
// this comment used to give -- that the collar argmax feeds the demographic
// growth-rate gradient and so must vary smoothly with plant state -- was right
// about the requirement and wrong about which solver meets it best:
//
//   * A fixed iteration COUNT is not smoothness. Golden section terminates on
//     bracket WIDTH, so it resolves the argmax only to `tol` and the residual
//     offset wanders discontinuously as the comparison sequence flips. Measured:
//     6 distinct answers across 11 trait steps, tread width ~GSS_tol_abs.
//   * That made the argmax piecewise constant at fine scales, so trait
//     derivatives came back exactly zero -- or, for traits in the hydraulic path,
//     smooth, plausible and SIGN-INVERTED.
//   * Solving dprofit == 0 instead resolves the argmax to solver precision and
//     measured ~1000x smoother second differences in a trait, and 24.5% faster
//     (2.65 vs 3.51 us/solve, interleaved at reps=2000).
//
// So the constraint stands and the conclusion inverted. Keep the constraint in
// mind before changing the collar solver again; do not read this function's
// existence as evidence that a comparison-based search is the safe choice.
template <typename Function>
double golden_section_max(Function f, double ax, double bx, double tol) {
  const double gr = (std::sqrt(5.0) + 1.0) / 2.0;  // ~1.6180339...
  double a = ax, b = bx;
  double c = b - (b - a) / gr;
  double d = a + (b - a) / gr;
  double fc = f(c);
  double fd = f(d);
  while (std::abs(b - a) > tol) {
    if (fc > fd) {
      b  = d;
      d  = c;
      fd = fc;                 // reuse
      c  = b - (b - a) / gr;
      fc = f(c);               // 1 new eval
    } else {
      a  = c;
      c  = d;
      fc = fd;                 // reuse
      d  = a + (b - a) / gr;
      fd = f(d);               // 1 new eval
    }
  }
  return (a + b) / 2.0;
}

// Maximise f over the CLOSED interval [lo, hi], where the maximum may sit at an
// endpoint and f need not be unimodal. Returns the argmax; `fmax`, if non-null,
// receives f there.
//
// ⚠️ WHY A BARE brent_fmin IS NOT THIS FUNCTION, and why every caller that can be
// maximised at a constraint must use this one instead. Brent steps in from the
// bounds, so it can return neither endpoint, and it follows one basin, so it
// cannot see past a local maximum. A bracketing optimiser answers "where is the
// interior maximum", which is a different question from "where is the maximum".
//
// Three parts, and each is load-bearing:
//
//   1. both endpoints are evaluated, so a constrained optimum is reachable;
//   2. an `n`-cell scan locates the basin, so a second interior hump cannot hide
//      the global one;
//   3. Brent refines inside the winning cell with a tolerance scaled to that
//      CELL, not to the whole interval.
//
// ⚠️ Part 3 is not a detail. Brent terminates on bracket width, so a fixed
// absolute tolerance comparable to the cell width leaves the answer at
// essentially the grid point -- and adding grid points then makes the result
// WORSE, because the cells get narrower while the tolerance does not. Measured
// over a 1728-row driver sweep against a 2001-point reference: with a fixed
// tolerance the shortfall count rose from 64 (n=8) to 104 (n=32); with the
// cell-scaled tolerance n=64 matches the reference on every row of both
// single-layer objectives.
//
// The grid argmax is always kept as a candidate, so the refinement can only
// improve on it.
template <typename Function>
double maximise_over_closed_interval(Function f, double lo, double hi, int n,
                                     double* fmax = nullptr) {
  double best_x = lo;
  double best_f = -std::numeric_limits<double>::infinity();
  auto consider = [&](double x, double fx) {
    if (std::isfinite(fx) && fx > best_f) {
      best_f = fx;
      best_x = x;
    }
  };

  if (!(hi > lo) || n < 2) {
    // Degenerate interval: the two endpoints are all there is to compare.
    consider(lo, f(lo));
    if (hi > lo) consider(hi, f(hi));
    if (fmax != nullptr) *fmax = best_f;
    return best_x;
  }

  // The scan, endpoints included -- so parts 1 and 2 share one loop rather than
  // evaluating the ends twice.
  int arg = 0;
  std::vector<double> xs(static_cast<std::size_t>(n) + 1);
  for (int i = 0; i <= n; ++i) {
    xs[static_cast<std::size_t>(i)] = lo + (hi - lo) * double(i) / double(n);
    const double fx = f(xs[static_cast<std::size_t>(i)]);
    if (std::isfinite(fx) && fx > best_f) arg = i;
    consider(xs[static_cast<std::size_t>(i)], fx);
  }

  // An endpoint argmax IS the answer, not a failure to search -- returning it
  // unrefined is the whole point of parts 1 and 2.
  if (arg == 0 || arg == n) {
    if (fmax != nullptr) *fmax = best_f;
    return best_x;
  }

  const double a = xs[static_cast<std::size_t>(arg - 1)];
  const double b = xs[static_cast<std::size_t>(arg + 1)];
  double neg = 0.0;
  const double x = brent_fmin([&](double v) { return -f(v); }, a, b,
                              (b - a) * 1e-4, &neg);
  consider(x, -neg);
  if (fmax != nullptr) *fmax = best_f;
  return best_x;
}

}
}

#endif
