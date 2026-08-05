# Trait gradients (issue #4, PLAN 11d stage 1).
#
# dA/dtheta, dgs/dtheta and dpsi_leaf/dtheta at a solved operating point, for the
# traits. The customer is `leaf-calibration`, whose hierarchical fit spent 97.6%
# of 10,045 model evaluations on finite-difference gradients.
#
# WHAT THIS FILE IS FOR, in one paragraph, because the arithmetic that motivated
# it did not survive being measured. The plan was that the implicit function
# theorem would replace 2N re-solves per gradient with 2N evaluations of the exact
# `dprofit`, at 7% of the cost of a solve, for a 10x win in the R layer. The
# theorem is right and it is implemented below. The 10x is not there, and the
# reason is worth stating once: in R the per-call overhead (~1.8 us) is seven
# times the C++ work it wraps (0.26 us), so replacing a 3.6 us solve with a
# 0.26 us gradient evaluation saves almost nothing measurable from R. What DOES
# dominate is that perturbing a trait used to mean constructing a new Leaf, which
# costs ~204 us -- 55 solves. `set_traits()` removed that, and it is where the
# speedup actually came from. See PLAN 11d for the numbers and for what stage 2
# has to do instead.
#
# So the value of the code below is EXACTNESS and the ACTIVE-SET CLASSIFICATION,
# not speed, and the classification is the part that took the work.

##' Replace the traits on an existing leaf
##'
##' Changes the traits of a `Leaf` in place, leaving [leaf_control()]'s numerical
##' settings alone. This is the fast path for anything that varies traits --
##' calibration, sensitivity analysis, a trait response curve -- because
##' constructing a `Leaf` from R costs about fifty times as much as solving one.
##'
##' @section You must set the drivers again afterwards:
##' This returns the leaf to its just-constructed state, so [set_drivers()] has to
##' be called before the next solve. That is not caution: `vcmax_`, `jmax_` and
##' `R_d_` are derived from the traits inside `set_physiology()`, so they are
##' genuinely unknown until the drivers are re-supplied.
##'
##' The reason it is a function rather than fifteen assignable fields is that
##' `leaf$vcmax_25 <- x` could not be made correct. Three separate pieces of
##' derived state go stale on a bare trait write: the two pre-integrated
##' vulnerability splines, the solved operating point, and -- least visibly --
##' `vcmax_`/`jmax_`/`R_d_`, which `set_physiology()` computes behind a cache keyed
##' on leaf temperature and O2 alone. That last one means the obvious repair,
##' "change the trait and then set the drivers again", silently does not work.
##'
##' @param x a `Leaf`, from [leaf_model()]
##' @param traits a [leaf_traits()] object
##'
##' @return `x`, invisibly.
##' @seealso [leaf_model()], [set_drivers()], [leaf_gradient()]
##' @examples
##' l <- leaf_model()
##' set_traits(l, leaf_traits(vcmax_25 = 120))
##' set_drivers(l, psi_soil = 2.0, PPFD = 900)
##' l$find_root_collar_psi()
##' operating_point(l)
##' @export
set_traits <- function(x, traits) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }
  # Positional, in the C++ argument order, for the same reason leaf_model() is:
  # the ordering is written down once rather than at every call site.
  x$set_traits(traits$vcmax_25, traits$stem_c, traits$stem_b, traits$psi_crit,
               traits$root_c, traits$root_b, traits$root_psi_crit, traits$beta2,
               traits$jmax_25, traits$a, traits$curv_fact_elec_trans,
               traits$curv_fact_colim, traits$cost_scale_TF24, traits$beta_R_H,
               traits$beta_R_V)
  invisible(x)
}

##' Trait gradients of a solved operating point
##'
##' The derivatives of the solved outputs with respect to the traits: `dA/dtheta`,
##' `dgc/dtheta`, `dpsi_stem/dtheta` and `dcollar/dtheta`, at one operating point.
##'
##' @section The maths, and why it is not just a finite difference:
##' The reported outputs are evaluated at the profit-maximising collar potential
##' `psi*`, so a trait moves them two ways: directly, and by moving `psi*`. At an
##' interior optimum `dprofit/dpsi = 0`, and differentiating that condition gives
##' `psi*`'s response without re-solving the model:
##'
##' \deqn{d\psi^*/d\theta = -M / H, \quad
##'       M = \partial^2 profit/\partial\psi\partial\theta, \quad
##'       H = \partial^2 profit/\partial\psi^2}
##'
##' \deqn{dY/d\theta = \partial Y/\partial\theta|_\psi +
##'                    (\partial Y/\partial\psi)(d\psi^*/d\theta)}
##'
##' Both second derivatives are differences of `dprofit_droot_collar_psi()`, which
##' is exact in `psi` and smooth in the traits -- `psi` is one of its arguments, so
##' no argmax is involved and `M` is a first difference of an exact quantity. It is
##' stable to seven significant figures across five decades of step size.
##'
##' The second term is not a correction. For `cost_scale_TF24`, `beta2`, `stem_b`
##' and `stem_c` it is 100% of the answer, and for `vcmax_25` 52%.
##'
##' @section The active set, which is the reason this function is careful:
##' Stationarity is the premise of the whole derivation, and it fails when the
##' optimum is pinned to an end of the feasible collar interval. There `psi*` is a
##' bound, `dprofit` is not zero at the answer, and `-M/H` is not the bound's
##' derivative. The formula does not fail loudly: measured at a wet-pinned point
##' the true gradient is about `1e-08` and the bare composite returns `O(1)` --
##' plausible-looking and wrong by seven orders of magnitude. About one operating
##' point in six is pinned across this package's test grid, all at soil suctions of
##' 3-4 MPa, which is the dry end a drought calibration will visit.
##'
##' So the premise is **tested rather than assumed**. The test is the implied
##' Newton step `|dprofit(psi*) / H|`, a distance in MPa: at a stationary point it
##' is zero to solver precision, and at a pinned one it is not. Over the package's
##' 288-point grid the two populations are five orders of magnitude apart -- worst
##' interior `5e-11`, smallest pinned `6e-06` -- so `stationarity_tol` sits in an
##' empty band rather than on a judgement call.
##'
##' When the premise does not hold, this falls back to a central finite difference
##' of the whole solve, which respects the constraint by construction because it
##' differences the constrained answer. That is the one thing the crude method does
##' better, and `method` in the result says which was used.
##'
##' It is better in a second way, which is sharper than accuracy. `psi_crit` does
##' not appear in the profit function at all -- it only sets the dry end of the
##' feasible collar interval -- so the composite necessarily returns **exactly
##' zero** for it. At an interior optimum that is correct. At a dry-pinned optimum
##' `psi_crit` *is* the binding constraint, the true `dA/dθ` is about `1.26`, and a
##' zero is arguably a worse answer than a wildly wrong one: it tells an optimiser
##' the parameter does nothing, and nothing about it looks suspicious. The
##' difference of the solve gets it for free.
##'
##' There is a second guard behind that one, found by measurement rather than
##' designed: a pinned optimum sits so close to its bound that the difference in
##' `psi` cannot be centred there without being clamped, and that is detected too.
##' Over this package's grid it catches every pinned point on its own, so forcing
##' `method = "ift"` at one fails loudly instead of returning the wrong number.
##' Do not read that as making the stationarity test redundant -- it fires only
##' when the feasible interval is narrow, which is true of these points but not of
##' pinned optima in general.
##'
##' @section Two of these are not traits:
##' `pars` also accepts `leaf_specific_conductance_max` and, on the
##' single-potential path, `resistance`. They are here because a calibration fits
##' them: of `leaf-calibration`'s four free parameters two are traits and two are
##' these, so restricting `pars` to [leaf_traits()] left half of that fit without
##' an exact gradient.
##'
##' Nothing in the derivation cares that the parameter is a trait — the implicit
##' function theorem is applied to `dprofit/dpsi = 0`, and any parameter the profit
##' function depends on goes through it identically. What differs is only which
##' setter applies the perturbation. `leaf_specific_conductance_max` is also how
##' plant's height reaches the leaf, so the same route gives a gradient with
##' respect to plant state and not only with respect to traits.
##'
##' ⚠️ **Both move the feasible collar interval, not just the profit function** —
##' `resistance` directly and `leaf_specific_conductance_max` through the supply
##' curve — so a perturbation can push a bound past `psi*` more readily than a
##' trait perturbation can. That is the case the clamp detector below exists for,
##' and it fails loudly rather than quietly.
##'
##' @section Precision:
##' Do not ask for more than about `1e-09` from any of this. That is the
##' achievable precision of the solved outputs themselves, set by the tolerance of
##' the intercellular-CO2 root-find.
##'
##' @inheritParams set_drivers
##' @param traits a [leaf_traits()] object: the point in trait space to
##'   differentiate at
##' @param control a [leaf_control()] object
##' @param supply how water reaches the root collar: [leaf_supply_multilayer()]
##'   (the default) or [leaf_supply_single()]
##' @param pars what to differentiate with respect to. Any of the fifteen
##'   [leaf_traits()] names, plus `"leaf_specific_conductance_max"` and — on the
##'   single-potential path only — `"resistance"`. Defaults to all of them.
##' @param step relative step for the trait difference. The default `1e-06` is
##'   near the middle of the five decades over which the mixed partial was
##'   measured stable; it is also used, relative to the collar potential, for the
##'   two differences in `psi`.
##' @param stationarity_tol how far from stationary, as an implied Newton step in
##'   MPa, the solved point may be before the optimum is called pinned and the
##'   finite-difference fallback is used. The default `1e-08` sits in the empty
##'   band between the two measured populations; there is no operating point in
##'   this package's grid within four orders of magnitude of it.
##' @param fast_stem_curve perturb `stem_b` without rebuilding the stem
##'   vulnerability spline (`TRUE`, the default). The curve is homogeneous of
##'   degree 1 in `stem_b`, so the spline for a perturbed `stem_b` is the existing
##'   one with its argument rescaled — which makes the rebuild, otherwise most of
##'   the cost of that gradient, unnecessary. `FALSE` rebuilds, and exists so the
##'   equivalence is checkable rather than assumed; the package's tests require
##'   the two routes to agree. `stem_c` has no such identity and always rebuilds.
##' @param method which route to take: `"auto"` (the default) uses the
##'   implicit-function composite where its premise holds and the finite
##'   difference where it does not. `"ift"` and `"fd"` force one route, for
##'   checking one against the other. Forcing `"ift"` at a pinned optimum returns
##'   a confidently wrong number -- that is the whole reason `"auto"` exists --
##'   and `status` still reports the truth about the point.
##'
##' @return A list with
##'   \describe{
##'     \item{`gradient`}{a matrix, one row per trait in `pars`, with columns
##'       `A` (umol C m^-2 s^-1 per trait unit), `gc`, `psi_stem` and `collar`}
##'     \item{`value`}{the solved outputs the gradient is taken at}
##'     \item{`method`}{`"ift"` if the implicit-function composite was used,
##'       `"fd"` if the fallback was}
##'     \item{`status`}{`"interior"`, `"pinned"` or `"no-gradient"`}
##'     \item{`H`}{the curvature of profit in the collar potential at `psi*`}
##'     \item{`stationarity`}{the implied Newton step, in MPa, that `method` was
##'       decided on}
##'   }
##'
##' @seealso [leaf_solve()] for the operating point itself, [set_traits()].
##' @examples
##' g <- leaf_gradient(psi_soil = 2.0, PPFD = 900)
##' g$method
##' g$gradient[c("vcmax_25", "stem_b"), ]
##'
##' # In dry soil at high VPD the optimum pins against the edge of the feasible
##' # range, the composite's premise fails, and the fallback is used instead.
##' dry <- leaf_gradient(psi_soil = 4.5, PPFD = 900, atm_vpd = 3.0)
##' dry$method       # "fd"
##' dry$status       # "pinned"
##' dry$stationarity # how far from stationary the optimum is, in MPa
##' @export
leaf_gradient <- function(psi_soil,
                          PPFD = 900,
                          soil_depth = NULL,
                          root_carbon_per_leaf_area = NULL,
                          leaf_specific_conductance_max = 3.14e-5,
                          atm_vpd = 2.0,
                          ca = 40.0,
                          leaf_temp = 25.0,
                          atm_o2_kpa = 21.0,
                          atm_kpa = 101.3,
                          traits = leaf_traits(),
                          control = leaf_control(),
                          supply = leaf_supply_multilayer(),
                          pars = NULL,
                          step = 1e-6,
                          stationarity_tol = 1e-8,
                          method = c("auto", "ift", "fd"),
                          fast_stem_curve = TRUE) {
  method <- match.arg(method)
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }
  if (!(is.numeric(step) && length(step) == 1L && step > 0)) {
    stop("`step` must be a single positive number", call. = FALSE)
  }

  drivers <- list(psi_soil = psi_soil, PPFD = PPFD, soil_depth = soil_depth,
                  root_carbon_per_leaf_area = root_carbon_per_leaf_area,
                  leaf_specific_conductance_max = leaf_specific_conductance_max,
                  atm_vpd = atm_vpd, ca = ca, leaf_temp = leaf_temp,
                  atm_o2_kpa = atm_o2_kpa, atm_kpa = atm_kpa)

  # The differentiable parameters: the fifteen traits, plus the two that are not
  # traits and that a calibration nonetheless fits (#44). See .gradient_theta.
  theta <- .gradient_theta(traits, leaf_specific_conductance_max, supply)
  if (is.null(pars)) {
    pars <- names(theta)
  }
  unknown <- setdiff(pars, names(theta))
  if (length(unknown)) {
    # Name the multi-layer case specially: `resistance` is a real parameter on
    # the other supply path, so "not differentiable" would be misleading rather
    # than merely unhelpful.
    why <- if (identical(supply$kind, "multilayer") &&
               "resistance" %in% unknown) {
      paste(" `resistance` is a parameter of the single-potential path only;",
            "on the multi-layer path the soil-to-collar resistances are derived",
            "from root carbon.")
    } else {
      ""
    }
    stop("`pars` names things this cannot differentiate: ",
         paste(unknown, collapse = ", "),
         ". Available: the traits from leaf_traits(), plus ",
         "`leaf_specific_conductance_max`",
         if (identical(supply$kind, "single")) " and `resistance`" else "",
         ".", why, call. = FALSE)
  }

  # ONE leaf for the whole gradient, re-traited rather than reconstructed. This is
  # the measurement that reordered PLAN 11d: a fresh Leaf costs ~204 us against
  # ~4 us to re-trait and re-drive this one, so reconstructing per perturbation
  # would swamp any difference between the two gradient routes below.
  l <- leaf_model(traits, control, supply)
  reset <- .gradient_setter(l, traits, drivers, supply, fast_stem_curve)
  do.call(set_drivers, c(list(l), drivers))
  l$find_root_collar_psi()

  psi_star <- l$opt_root_psi_
  value <- .gradient_outputs(l)

  # Is the premise true here? See the active-set section: the test is the implied
  # Newton step, which is a distance in MPa and so needs no scale of its own.
  h_psi <- max(abs(psi_star), 1) * step
  resid <- l$dprofit_droot_collar_psi(psi_star)
  H <- (l$dprofit_droot_collar_psi(psi_star + h_psi) -
        l$dprofit_droot_collar_psi(psi_star - h_psi)) / (2 * h_psi)
  # H == 0 with resid == 0 is the shut-down signature: dprofit returns a sentinel
  # zero there rather than a derivative, so the ratio below would be 0/0. H > 0
  # would not be a maximum. Both mean the composite has nothing to stand on.
  usable <- is.finite(H) && H < 0 && is.finite(resid)
  stationarity <- if (usable) abs(resid / H) else Inf
  status <- if (!usable) "no-gradient" else
    if (stationarity > stationarity_tol) "pinned" else "interior"

  # `status` describes the POINT and is reported whichever route runs; `use_ift`
  # is the route. They differ only when the caller has forced one.
  use_ift <- switch(method,
                    auto = identical(status, "interior"),
                    ift = TRUE,
                    fd = FALSE)
  if (use_ift && !usable) {
    stop("leaf_gradient(): method = \"ift\" was asked for at a point with no ",
         "usable curvature (H = ", format(H), "), so -M/H has nothing to stand ",
         "on. This is a shut-down or otherwise determined operating point; use ",
         "method = \"auto\".", call. = FALSE)
  }

  if (use_ift) {
    # dY/dpsi at fixed traits. evaluate_root_collar_psi CLAMPS its target into the
    # feasible interval, so check it landed where it was asked to: a clamped
    # evaluation would silently make this a one-sided difference over a shorter
    # interval, which is the same class of error as the pinned case.
    #
    # This turns out to be a SECOND, INDEPENDENT detector of a pinned optimum
    # rather than the unreachable guard it was written as, and the measurement is
    # worth recording. At a pinned point psi* sits one step-in fraction (1e-06 of
    # the bracket width) from the bound, so a step of `step` * psi in psi crosses
    # it whenever the bracket is narrower than psi -- which every pinned row in
    # the package's grid is, being at the dry end where the feasible interval has
    # nearly closed. Measured: forcing method = "ift" fails here on all 42 pinned
    # rows and on all 48 shut-down ones, so the composite's silently-wrong answer
    # is not reachable through this function at all.
    #
    # It is NOT a substitute for the stationarity test, and reading it as one
    # would be the mistake: it fires only when the bracket is narrow, so a pinned
    # optimum on a wide bracket would pass it. The stationarity test is the one
    # that is scale-free and the one the default relies on.
    hi <- .gradient_outputs_at(l, psi_star + h_psi)
    lo <- .gradient_outputs_at(l, psi_star - h_psi)
    if (is.null(hi) || is.null(lo)) {
      if (identical(method, "ift")) {
        stop("leaf_gradient(): method = \"ift\" was asked for at a point whose ",
             "feasible collar interval is narrower than one step, so dY/dpsi ",
             "cannot be centred on psi*. Use method = \"auto\".", call. = FALSE)
      }
      use_ift <- FALSE
      status <- "pinned"
    } else {
      dY_dpsi <- (hi - lo) / (2 * h_psi)
    }
  }

  grad <- if (use_ift) {
    .gradient_ift(l, reset, theta, pars, psi_star, H, dY_dpsi, step)
  } else {
    .gradient_fd(l, reset, theta, pars, step)
  }

  list(gradient = grad,
       value = value,
       method = if (use_ift) "ift" else "fd",
       status = status,
       H = H,
       stationarity = stationarity)
}

# --- internals ---------------------------------------------------------------

# The four differentiated outputs, in one place so the composite and the fallback
# cannot disagree about what they are. `collar` is psi* itself, which makes
# dcollar/dtheta equal to dpsi*/dtheta -- so the two routes compute the same
# quantity by different means, and a test can compare them.
.gradient_output_names <- c("A", "gc", "psi_stem", "collar")

.gradient_outputs <- function(l) {
  c(A = l$assim_colimited_, gc = l$stom_cond_CO2_, psi_stem = l$opt_psi_stem_,
    collar = l$opt_root_psi_)
}

# Outputs with the collar held at `psi` rather than optimised. NULL when the
# clamp moved the target, because then this is not the evaluation that was asked
# for. Exact equality is the right test: the clamp is a min/max, so an unclamped
# target comes back bit-identical.
.gradient_outputs_at <- function(l, psi) {
  l$evaluate_root_collar_psi(psi)
  if (l$opt_root_psi_ != psi) {
    return(NULL)
  }
  .gradient_outputs(l)
}

# Everything this can differentiate, and its current value, as one named vector.
#
# The fifteen traits, plus the two quantities a calibration fits that are NOT
# traits (#44) and that `pars` therefore used to reject:
#
#   * `leaf_specific_conductance_max` -- a DRIVER, set through set_drivers(). It
#     is also how plant's height reaches the leaf, so the same code path answers
#     plant #537's question about differentiating w.r.t. state.
#   * `resistance` -- the single-potential path's whole soil-to-collar
#     resistance, set through $set_supply_single(). Only present on that path: on
#     the multi-layer path the resistances are derived from root carbon and there
#     is no such parameter, which is why it appears here conditionally rather
#     than being rejected later with a worse message.
#
# Nothing in the derivation cares that theta is a trait -- `dpsi*/dtheta = -M/H`
# and `dY/dtheta = dY/dtheta|_psi + (dY/dpsi)(dpsi*/dtheta)` hold for any
# parameter profit depends on. What differs per parameter is only which setter
# applies it, which is .gradient_setter's job.
.gradient_theta <- function(traits, kmax, supply) {
  theta <- c(unlist(traits), leaf_specific_conductance_max = kmax)
  if (identical(supply$kind, "single")) {
    theta <- c(theta, resistance = supply$resistance)
  }
  theta
}

# Push a parameter vector back onto the leaf, in the one order that is correct.
#
# ⚠️ ORDER IS LOAD-BEARING and the reason this is a function rather than three
# lines at each call site. `set_traits()` and `$set_supply_single()` both return
# the leaf to its just-constructed state, so the drivers have to be re-supplied
# AFTER them -- and `set_drivers()` is what re-derives vcmax_/jmax_/R_d_ behind
# the temperature cache that `set_traits()` has just cleared.
.gradient_setter <- function(l, traits, drivers, supply,
                             fast_stem_curve = TRUE) {
  trait_names <- names(traits)
  trait_class <- class(traits)
  single <- identical(supply$kind, "single")
  function(theta, only = NULL) {
    # THE FAST PATH FOR stem_b, and it is the whole of PLAN 11f.
    #
    # stem_b owns the pre-integrated stem vulnerability spline, so moving it
    # through set_traits() rebuilds it -- 11.9 us of incomplete gammas plus two
    # interpolator inits, which is essentially the entire cost of its gradient.
    # But the curve is homogeneous of degree 1 in stem_b, so the spline for
    # another stem_b IS this one with its argument rescaled, and no rebuild is
    # needed at all.
    #
    # ⚠️ stem_c is NOT here, and that is a measured decision rather than an
    # omission: it has no such identity, and reading the curve from its closed
    # form instead -- the obvious alternative -- differentiates a slightly
    # different model and disagrees by 3e-4. PLAN 11f has the numbers.
    #
    # ⚠️ Sound only because `only` names a SINGLE parameter, so everything else in
    # `theta` is still what the object was last set to. That is true here and
    # nowhere else, which is why the argument exists rather than the function
    # guessing: .gradient_ift and .gradient_fd perturb exactly one parameter at a
    # time, and the next parameter's first call goes through the full path below,
    # rebuilding the spline on its way.
    if (fast_stem_curve && identical(only, "stem_b")) {
      l$perturb_stem_b(theta[["stem_b"]])
      return(invisible(l))
    }
    set_traits(l, structure(as.list(theta[trait_names]), class = trait_class))
    if (single) {
      l$set_supply_single(theta[["resistance"]], supply$gravity_head)
    }
    drivers$leaf_specific_conductance_max <-
      theta[["leaf_specific_conductance_max"]]
    do.call(set_drivers, c(list(l), drivers))
    invisible(l)
  }
}

# Parameters whose step is RELATIVE rather than floored at 1. See .gradient_step.
.gradient_relative_pars <- c("leaf_specific_conductance_max", "resistance")

# The step: relative to the parameter for values above 1, and plain `step` below
# it. Not a relative step with an epsilon floor -- the floor is at 1, which is
# deliberate. Traits here span `a` = 0.3 to `beta_R_V` = 9.4e3, and a strictly
# relative step would perturb the small ones so little that the difference is
# dominated by the solve's ~1e-09 noise. It also means a trait sitting at zero is
# still perturbed rather than not differentiated at all.
#
# ⚠️ That floor is WRONG for a parameter whose natural magnitude is far below 1,
# and `leaf_specific_conductance_max` is: it defaults to 3.14e-05, so flooring at
# 1 would perturb it by **3%** and measure a secant across a range over which the
# model is visibly nonlinear, not a derivative. Those parameters get a plain
# relative step. For `resistance` (~1e3 and up) the two rules coincide; it is
# listed for the reason rather than for the arithmetic.
.gradient_step <- function(par, value, step) {
  floor <- if (par %in% .gradient_relative_pars) 0 else 1
  max(abs(value), floor) * step
}

# The implicit-function composite. Two perturbed evaluations per parameter,
# neither of which re-solves the model: `dprofit` at the UNPERTURBED psi* gives
# the mixed partial, and the outputs at that same psi* give the direct term.
.gradient_ift <- function(l, reset, theta, pars, psi_star, H, dY_dpsi, step) {
  out <- t(vapply(pars, function(p) {
    h <- .gradient_step(p, theta[[p]], step)
    side <- function(sign) {
      th <- theta
      th[[p]] <- th[[p]] + sign * h
      reset(th, only = p)
      y <- .gradient_outputs_at(l, psi_star)
      if (is.null(y)) {
        stop("leaf_gradient(): perturbing `", p, "` moved the feasible collar ",
             "interval past psi*, so the operating point could not be ",
             "evaluated there. This point is on an active-set boundary; ",
             "lower `stationarity_tol` or difference the solve directly.",
             call. = FALSE)
      }
      c(dprofit = l$dprofit_droot_collar_psi(psi_star), y)
    }
    up <- side(1)
    dn <- side(-1)
    # M = d2profit/dpsi dtheta, with psi held FIXED at psi*.
    dpsi_dtheta <- -((up[["dprofit"]] - dn[["dprofit"]]) / (2 * h)) / H
    direct <- (up[-1L] - dn[-1L]) / (2 * h)
    # `collar` is not an output of the evaluation -- it IS psi*, held fixed, so
    # its direct term is zero by construction and the composite reduces to
    # dpsi*/dtheta. Setting it explicitly says so, rather than relying on the
    # difference of two identical numbers.
    g <- direct + dY_dpsi * dpsi_dtheta
    g[["collar"]] <- dpsi_dtheta
    g[.gradient_output_names]
  }, numeric(length(.gradient_output_names))))
  reset(theta)
  out
}

# The fallback: a central difference of the whole solve. Correct at a pinned
# optimum because it differences the constrained answer, which is exactly what the
# composite cannot do.
.gradient_fd <- function(l, reset, theta, pars, step) {
  out <- t(vapply(pars, function(p) {
    h <- .gradient_step(p, theta[[p]], step)
    side <- function(sign) {
      th <- theta
      th[[p]] <- th[[p]] + sign * h
      reset(th, only = p)
      l$find_root_collar_psi()
      .gradient_outputs(l)
    }
    ((side(1) - side(-1)) / (2 * h))[.gradient_output_names]
  }, numeric(length(.gradient_output_names))))
  reset(theta)
  out
}
