# Trait gradients at a solved operating point: set_traits() and leaf_gradient().
#
# leaf_gradient() is leaf_gradient_batch() with one observation: this file
# handles its arguments and the shape of its result, and the gradient itself is
# inst/include/phylloptim/gradient.hpp's, through `.gradient_run()`.
#
# When the exact gradient beats differencing is ?leaf_gradient's "What it
# costs"; the measured per-call costs are README's "Performance from R".

##' Replace the traits on an existing leaf
##'
##' Changes the traits of a `Leaf` in place, leaving [leaf_control()]'s numerical
##' settings alone. This is the fast path for anything that varies traits --
##' calibration, sensitivity analysis, a trait response curve -- because
##' constructing a `Leaf` from R costs about 45 times as much as solving one.
##'
##' @section You must set the drivers again afterwards:
##' This returns the leaf to its just-constructed state, so [set_drivers()] has to
##' be called before the next solve. That is not caution: `vcmax_`, `jmax_` and
##' `R_d_` are derived from the traits inside `set_physiology()`, so they are
##' genuinely unknown until the drivers are re-supplied.
##'
##' The reason it is a function rather than sixteen assignable fields is that
##' `leaf$vcmax_25 <- x` could not be made correct. Three separate pieces of
##' derived state go stale on a bare trait write: the two pre-integrated
##' vulnerability tables, the solved operating point, and -- least visibly --
##' `vcmax_`/`jmax_`/`R_d_`, which `set_physiology()` computes behind a cache. The
##' cache is keyed on the traits it reads, so setting the drivers again refreshes
##' it, but nothing short of `set_traits()` refreshes the tables or the solved point.
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
##' l$optimise()
##' operating_point(l)
##' @export
set_traits <- function(x, traits) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }
  # Positional, in the C++ setter's argument order; a trait added there is an
  # argument added here.
  x$set_traits(traits$vcmax_25, traits$stem_c, traits$stem_P50,
               traits$root_c, traits$root_P50, traits$TF24_beta2,
               traits$jmax_25, traits$phi_J, traits$curv_fact_elec_trans,
               traits$curv_fact_colim, traits$TF24_cost_scale, traits$R_d_25,
               traits$JS22_gamma, traits$CMax_a, traits$CMax_b, traits$g_min)
  invisible(x)
}

##' Trait gradients of a solved operating point
##'
##' The derivatives of the solved outputs with respect to the traits:
##' `dA/dtheta`, `dgc/dtheta`, `dpsi_leaf/dtheta`, `dcollar/dtheta` and
##' `dprofit/dtheta`, at one operating point.
##'
##' @section Choosing which model to differentiate:
##' A model is a cost curve and a route, and every curve has a gradient on both
##' routes. `model` names the curve and `route` the potential it is maximised
##' over: `"collar"`, the root-collar potential, on either supply path; or
##' `"leaf"`, `psi_leaf` with the upstream potential pinned at `psi_soil`, which
##' needs the single-potential supply path. The defaults are one set,
##' everywhere: `"TF24_floor"` (plant's model; at its default price,
##' `TF24_floor_lambda_o = 0`, it is TF24) on the leaf route over the
##' single-potential supply. The collar route and the multi-layer supply are
##' asked for.
##'
##' The product objectives, `SOX` and `JW26`, need no special case. They maximise
##' `A * g(psi)`, which is the `Log` benefit link applied to the same
##' difference form everything else uses: `d/dpsi[h(A) - C] =
##' h'(A) * dA/dpsi - dC/dpsi` with `h' = 1/A`. One derivative serves all nine.
##'
##' ⚠️ **The two routes are different models, not two ways of asking one
##' question**, even for one curve. They optimise different variables over
##' different supply topologies, for TF24 disagree on 20 of 30 driver rows, and
##' the gap does not close as the root resistance goes to zero.
##'
##' ⚠️ **`model` SEATS the leaf, including one you passed as `x`.** The curve
##' and route are applied with `$set_model()` and left there when this returns,
##' alongside the base-point traits and solve. A loop that alternates models on one `x` is fine;
##' a loop that takes a gradient and then reads `$optimise()` off the same object
##' expecting the default is not.
##'
##' ⚠️ **On a leaf route the decision variable is `psi_leaf`, so `collar` comes
##' back non-finite.** The collar is not solved for at all there, and an absent
##' number is the truthful answer rather than a gap -- reading one would be reading
##' state from whichever solve ran last.
##'
##' **Every route reaches the composite under `method = "auto"`, and every route
##' uses the same finite-difference step (1e-06).** Both were once untrue.
##' `ProfitMax` normalises by `|A|max`, and while that came from a 500-point scan
##' its argmax was piecewise constant in the parameters -- so the composite had to
##' pin it and returned a PARTIAL (0.0766 against a true 0.0579), and differencing
##' needed a step large enough to cross a scan cell. `|A|max` is found by a
##' root-find now, so it is re-solved per perturbation, the composite returns the
##' TOTAL, and no model carries that split.
##'
##' ⚠️ **A zero column for another curve's parameter means two different things,
##' and they are not distinguishable from the number.** On a `collar` or `TF24`
##' route, `JS22_gamma`, `CMax_a`, `CMax_b` and the rest come back exactly zero,
##' and read one way that is correct: those parameters do not enter the TF24
##' cost, so their derivative is genuinely zero. Read the other way it is a trap:
##' if you are *fitting* one of those curves, the zero is telling you this
##' differentiated a different model, not that the parameter does nothing.
##' Measured on the single-potential path at `psi_soil` 1.5, the TF24 collar
##' optimum sits at 2.09 MPa and the JS22 optimum at 3.15.
##'
##' An optimiser handed such a column concludes the parameter is inert and leaves
##' it at its starting value, which looks like convergence. **Name the curve you
##' are fitting in `model`**; `.gradient_available_pars()` reports which
##' parameters a given route can move.
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
##' The second term is not a correction. For `TF24_cost_scale`, `TF24_beta2` and
##' `stem_c` it is 100% of the answer, and for `vcmax_25` 52%.
##'
##' @section profit, which is the one output the envelope theorem reaches:
##' The first four outputs are what a gas-exchange calibration observes. `profit`
##' is here because it is what a **demographic** consumer bills: plant's carbon
##' budget reads `leaf.profit_`, not `assim_colimited_`, so before it was added
##' the two sets were disjoint and no gradient from this package reached a
##' demographic model at all.
##'
##' It is also the cheapest of the five, because it is the objective. At an
##' interior optimum `dprofit/dpsi = 0`, so the indirect term above vanishes
##' identically and
##'
##' \deqn{dprofit/d\theta = \partial profit/\partial\theta|_{\psi^*}}
##'
##' -- the direct partial alone, with no `dY/dpsi` and no `-M/H`. That is the
##' envelope theorem, and it is the only place this package uses it.
##'
##' ⚠️ **The term this drops is noise, not an `h^2` truncation, and that is why
##' dropping it matters.** `profit` is the maximum, so it is flat, and a central
##' difference of it at `psi*` divides the solve's ~1e-09 floor by a ~1e-06 step.
##' Over the golden grid's 136 interior rows:
##'
##' | | median | max |
##' |---|---|---|
##' | `\|dprofit/dpsi\|`, exact (forward AD) | 4.8e-15 | 5.4e-10 |
##' | `\|dprofit/dpsi\|`, central difference | 7.8e-10 | 2.1e-04 |
##' | relative move in `dprofit/dtheta` if kept | 2.7e-10 | 8.0e-05 |
##'
##' Eleven orders between the two instruments at the median, and the worst row
##' sits in the band this package calls a real difference. So `dprofit/dtheta` is
##' set from the direct term rather than computed through the composite, exactly
##' as `collar` is set from `dpsi*/dtheta` rather than differenced.
##'
##' ⚠️ **The envelope does not survive the active set.** At a pinned optimum
##' `psi*` is a trait-dependent bound and `dprofit/dpsi` is not zero there, so the
##' indirect term is real. This applies the identity only where `status` is
##' `"interior"`; everywhere else `profit` goes through the same finite-difference
##' fallback as the other four.
##'
##' @section A collar potential you supply, instead of the one this solves for:
##' Pass `psi` and the outputs are evaluated at that collar rather than at the
##' argmax. This is for a caller whose model **tracks** the optimum instead of
##' finding it — plant's TF24f carries the collar potential as an ODE state,
##' `dpsi/dt = k * dprofit/dpsi`, so at finite gain its `dprofit/dpsi` is
##' deliberately non-zero and *is* the acclimation rate.
##'
##' There the derivation above simplifies rather than breaks. `psi` is exogenous,
##' so the indirect term is whatever the caller says it is:
##'
##' \deqn{dY/d\theta = \partial Y/\partial\theta|_\psi +
##'                    (\partial Y/\partial\psi)(d\psi/d\theta)}
##'
##' with `dpsi_dtheta` supplied — defaulting to zero, the partial at fixed
##' collar. Nothing is derived from `-M/H`, so nothing needs stationarity, and
##' `method` is refused because the two routes it chooses between are both about
##' a solved optimum.
##'
##' `M`, `H` and `dY_dpsi` come back in the result, because a caller with a
##' *dynamic* `psi` cannot supply `dpsi_dtheta` as a constant: for the
##' gradient-ascent law above it obeys `ds/dt = k(M + H s)`, and those are its
##' coefficients. Its fixed point is `-M/H`, which is what the solving path
##' returns — the two agree in the large-gain limit.
##'
##' ⚠️ **`stationarity` is still reported, and is used for exactly one thing.**
##' It no longer routes anything; it now measures how far the collar you supplied
##' sits from the optimum. The one decision it makes is `profit`'s: at a
##' stationary point the envelope theorem applies and the analytic zero is used,
##' and away from one the exact `dprofit/dpsi` is used instead of differencing
##' it. So `psi = <the solved psi*>` with `dpsi_dtheta = -M/H` reproduces the
##' solving path **bit-for-bit**, which the tests assert.
##'
##' ⚠️ **An infeasible `psi` gets no gradient either**, reported as
##' `"no-gradient"`. `dprofit/dpsi` is a sentinel `0.0` on the shut-down and
##' reversed-gradient exits rather than a derivative, and this path never divides
##' by `H`, so it would otherwise adopt that zero as the real thing.
##'
##' ⚠️ **A `psi` outside the feasible interval is refused, and gets no
##' gradient.** The leaf does not project it onto the interval's end, which
##' would make the outputs move with the *bound* rather than with `dpsi_dtheta`.
##' `status` reports `"refused"`, `psi` in the result is the one you gave, and
##' everything else is `NA`. This is reported rather than thrown because a fit
##' reaches these points routinely: a trait perturbation moves the interval past
##' a fixed `psi`. A `psi` inside the interval but within one step of an end is
##' `"refused"` too, with its value, because `dY/dpsi` cannot be centred there.
##'
##' @section The active set, which is the reason this function is careful:
##' Stationarity is the premise of the whole derivation, and it fails when the
##' optimum is pinned to an end of the feasible collar interval. There `psi*` is a
##' bound, `dprofit` is not zero at the answer, and `-M/H` is not the bound's
##' derivative. The formula does not fail loudly: measured at a boundary-soil point
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
##' It is better in a second way, which is sharper than accuracy. The dry end of
##' the feasible collar interval -- `psi_crit` for the curves it bounds, the end of
##' the stem curve for the rest -- does not appear in the profit function at all,
##' so the composite necessarily returns **exactly zero** for it. At an interior
##' optimum that is correct. At a dry pin it *is* the binding constraint, and a
##' zero is arguably a worse answer than a wildly wrong one: it tells an optimiser
##' the parameter does nothing, and nothing about it looks suspicious. The
##' difference of the solve gets it for free.
##'
##' There is a second guard behind that one, found by measurement rather than
##' designed: a pinned optimum sits so close to its bound that the difference in
##' `psi` cannot be centred there without stepping out of the interval, and that is
##' detected too.
##' Over this package's grid it catches every pinned point on its own, so forcing
##' `method = "ift"` at one fails loudly instead of returning the wrong number.
##' Do not read that as making the stationarity test redundant -- it fires only
##' when the feasible interval is narrow, which is true of these points but not of
##' pinned optima in general.
##'
##' @section Two of these are not traits:
##' `pars` also accepts `leaf_specific_conductance_max` and, on the
##' single-potential path's collar route, `resistance`. They are here because a calibration fits
##' them: of a companion calibration study's four free parameters two are traits and two are
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
##' trait perturbation can. That is the case the refusal detector below exists for,
##' and it fails loudly rather than quietly.
##'
##' @section What it costs, and when that beats differencing:
##' Both routes are linear in a parameter count, and it is a **different** count
##' for each. With `N` observations:
##'
##' \deqn{T_{FD} = 2\,P_{fit}\,T_{objective}}
##' \deqn{T_{exact} = N\,(t_{construct} + t_{solve}) + 2\,P_{model}\,N\,t_{eval}}
##'
##' An optimiser differences the `P_fit` parameters it is moving. This function is
##' asked for the `P_model` parameters the leaf has — the length of `pars`. Those
##' are equal only when the fit varies traits directly; a parameterisation in
##' between (pooling, a hierarchy, any shared or derived parameter) makes
##' `P_fit > P_model`, and then this route's cost does not follow the optimiser's
##' dimension.
##'
##' ⚠️ **So "is the exact gradient faster" has no answer independent of your
##' parameterisation, and both regimes are measured.** `vignette("fitting")` fits
##' 72 observations with `P_fit == P_model` and this loses at every count it
##' measures.
##' A companion calibration study has `N = 1327`, `P_fit = 40`,
##' `P_model = 4`, wins 3.4x, and reaches the same optimum — with its
##' 57-parameter variant costing the same as its 40-parameter one, because
##' `P_model` is 4 in both.
##'
##' ⚠️ **Always pass `pars`.** It is `P_model`, so the default — every trait plus
##' `leaf_specific_conductance_max` on the multi-layer path — is the most expensive
##' thing you can ask for, and a fit that reads four of them pays for the rest.
##'
##' The intercept is a fresh `Leaf` per call. Pass `x` to reuse one across calls
##' (#52); without it, construction is most of a one-observation gradient.
##'
##' @section Precision:
##' Do not ask for more than about `1e-09` from any of this. That is the
##' achievable precision of the solved outputs themselves, set by the tolerance of
##' the intercellular-CO2 root-find.
##'
##' @section Reusing a leaf across gradients:
##' Constructing a `Leaf` is the largest single cost on this surface, and this
##' function does it once per call: it is **most of a one-parameter gradient** on a
##' fresh leaf, paid once per observation by a fit that differentiates per
##' observation. The measured ratios are in the README's "Performance from R", the
##' one place the package quotes them. Pass `x` to reuse one:
##'
##' ```r
##' l <- leaf_model(traits, control, supply)
##' for (i in seq_len(n)) leaf_gradient(psi_soil = ps[i], x = l, traits = traits, pars = FIT)
##' ```
##'
##' `x` is a **vessel, not the point**: `traits` still says where the gradient is
##' taken and is applied to `x`, rather than read from it (a `Leaf` does not expose
##' its traits). That is why `traits` is required with `x` — defaulting it would
##' differentiate at the package defaults on somebody else's leaf and return a
##' plausible answer for a point nobody asked about.
##'
##' `control` and `supply` are fixed when `x` is built, so passing them alongside it
##' is an error rather than silently resolved.
##'
##' **`x` is left solved at the base point**, not at the last perturbation — the
##' state this function seats it in, restored on the way out including on an error.
##' The restore costs one solve, a small fraction of the construction it saves.
##'
##' @inheritParams set_drivers
##' @param x an existing `Leaf` from [leaf_model()] to reuse instead of building
##'   one, for a loop that takes many gradients. `traits` is required with it, and
##'   `control`/`supply` must be omitted. See the section above.
##' @param traits a [leaf_traits()] object: the point in trait space to
##'   differentiate at
##' @param control a [leaf_control()] object
##' @param supply how water reaches the root collar: [leaf_supply_singlelayer()]
##'   (the default) or [leaf_supply_multilayer()]
##' @param model which optimality model to differentiate: any name from
##'   [cost_curve_names()], on `route`; `"TF24_floor"` by default. See *Choosing
##'   which model to differentiate*.
##' @param route `"collar"` or `"leaf"`, the potential the curve is maximised
##'   over; `NULL` (the default) is the leaf route. See *Choosing which model
##'   to differentiate*.
##' @param pars what to differentiate with respect to. Any of the sixteen
##'   [leaf_traits()] names, plus `"leaf_specific_conductance_max"`, on the
##'   single-potential path's collar route only `"resistance"` (the leaf route
##'   never reads it), on `model = "CF77"` only
##'   `"CF77_lambda_"`, on `model = "TF24_floor"` only `"TF24_floor_lambda_o"`,
##'   and on `model = "LeastCost"` only `"LeastCost_beta"`.
##'   On a reused leaf `x` with `use_thermal_damage_` on, also the layer's inputs
##'   that reach the solve: its states `"A_"`, `"I_"`, `"D_"`, the optimum shift
##'   `"topt_offset_"`, `"dTopt_max_"`, `"K_A_"`, and every quantity a maintenance
##'   cost multiplies (`"k_i_"`, `"k_i_baseline_"`, `"k_rec_"`, `"k_rec_Ea_"`,
##'   `"c_acclim_maint_"`, `"c_repair_maint_"`,
##'   `"c_protect_maint_"`, `"c_repair_flux_"`). The layer's dynamics parameters
##'   move no solved output and are refused; their gradient is
##'   [leaf_damage_gradient()]'s. Defaults to all the route allows.
##'
##'   ⚠️ At 25 C the peaked Arrhenius is normalised to its reference, so
##'   `topt_offset_` and `dTopt_max_` have an exactly zero effect there -- a true
##'   zero, not a missing row.
##'   The protection cost is signed, `k_i_baseline_ - k_i_`, so its columns are
##'   the same linear slope on both sides of the default, where the two are
##'   equal.
##'
##'   `dY/dvcmax_25` is a PARTIAL at fixed respiration: `R_d_25` is its own trait,
##'   so a fit that wants respiration to follow Vcmax moves both and adds the two
##'   columns.
##'
##'   `beta_R_H` and `beta_R_V` were here until #33 and are not any more: they
##'   parameterise the root-architecture model, which the leaf no longer runs, so
##'   this function has no way to reach them. Getting a gradient in either means
##'   solving at two networks yourself. The perturbation itself is cheap —
##'   [root_network_from_carbon()] is homogeneous of degree 1 in each constant, so
##'   scaling `beta_R_H` scales `r_R_H_min` by the same factor and needs no
##'   rebuild — but the two solves are still two solves.
##' @param psi a collar water potential to evaluate at, in MPa as a positive
##'   magnitude, instead of solving for the profit-maximising one. `NULL` (the
##'   default) solves. See the section above; `method` cannot be given with it.
##' @param dpsi_dtheta how the prescribed `psi` itself responds to each
##'   parameter, as one value per entry of `pars` (or one value recycled, or a
##'   vector named by `pars`). Defaults to zero, which is the partial derivative
##'   at fixed collar. Must be finite — an infinite value would return an
##'   all-infinite gradient row with a finite, plausible-looking `profit` in it,
##'   because the envelope theorem assigns that column rather than multiplying
##'   through. Only meaningful with `psi`.
##' @param step relative step for the trait difference. The default `1e-06` is
##'   near the middle of the five decades over which the mixed partial was
##'   measured stable; it is also used, relative to the collar potential, for the
##'   two differences in `psi`.
##' @param stationarity_tol how far from stationary, as an implied Newton step in
##'   MPa, the solved point may be before the optimum is called pinned and the
##'   finite-difference fallback is used. The default `1e-08` sits in the empty
##'   band between the two measured populations; there is no operating point in
##'   this package's grid within four orders of magnitude of it.
##' @param fast_stem_curve perturb `stem_P50` through `perturb_stem_P50()`, which
##'   moves `stem_b` and leaves the seated physiology alone (`TRUE`, the
##'   default), rather than through `set_traits()`. Neither route rebuilds the
##'   stem's vulnerability table, which depends on `stem_c` alone; `FALSE` exists
##'   so the equivalence is checkable rather than assumed, and the package's tests
##'   require the two routes to agree.
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
##'       `A` (umol C m^-2 s^-1 per trait unit), `gc`, `psi_leaf`, `collar` and
##'       `profit` (also umol C m^-2 s^-1 per trait unit)}
##'     \item{`value`}{the solved outputs the gradient is taken at}
##'     \item{`method`}{`"ift"` if the implicit-function composite was used,
##'       `"fd"` if the fallback was, `"prescribed"` if `psi` was given}
##'     \item{`status`}{`"interior"`, `"pinned"` or `"no-gradient"` on the
##'       solving path; `"prescribed"`, `"refused"` or `"no-gradient"` with
##'       `psi`}
##'     \item{`H`}{the curvature of profit in the collar potential}
##'     \item{`stationarity`}{the implied Newton step, in MPa: what `method` was
##'       decided on, or with `psi` how far the collar you gave is from the
##'       optimum}
##'     \item{`M`}{the mixed partial `d2profit/dpsi dtheta`, one per entry of
##'       `pars`, or `NA` where the fallback ran}
##'     \item{`dY_dpsi`}{`dY/dpsi` at fixed traits, one per output. `profit`'s
##'       entry is the analytic zero at a stationary point and the exact
##'       `dprofit/dpsi` elsewhere -- see the envelope section}
##'     \item{`psi`}{the collar the outputs were evaluated at: `psi*` on the
##'       solving path, and the `psi` argument on the prescribed path, refused
##'       or not}
##'     \item{`model`, `route`}{the model differentiated and its route}
##'   }
##'
##'   An entry that is not determined, such as a refused point's gradient or a
##'   leaf route's `collar`, is `NaN`.
##'
##' @seealso [leaf_solve()] for the operating point itself, [set_traits()].
##' @examples
##' g <- leaf_gradient(psi_soil = 2.0, PPFD = 900)
##' g$method
##' g$gradient[c("vcmax_25", "stem_P50"), ]
##'
##' # On the collar route, in dry soil at high VPD, the optimum pins against the
##' # edge of the feasible range, the composite's premise fails, and the fallback
##' # is used instead.
##' dry <- leaf_gradient(psi_soil = 4.5, PPFD = 900, atm_vpd = 3.0,
##'                      supply = leaf_supply_multilayer(), route = "collar")
##' dry$method       # "fd"
##' dry$status       # "pinned"
##' dry$stationarity # how far from stationary the optimum is, in MPa
##' @export
leaf_gradient <- function(psi_soil,
                          PPFD = 900,
                          soil_depth = NULL,
                          root_network = NULL,
                          leaf_specific_conductance_max = 3.14e-5,
                          atm_vpd = 2.0,
                          ca = 40.0,
                          leaf_temp = 25.0,
                          atm_o2_kpa = 21.0,
                          atm_kpa = 101.3,
                          x = NULL,
                          traits = leaf_traits(),
                          control = leaf_control(),
                          supply = leaf_supply_singlelayer(),
                          pars = NULL,
                          psi = NULL,
                          dpsi_dtheta = NULL,
                          step = 1e-6,
                          stationarity_tol = 1e-8,
                          method = c("auto", "ift", "fd"),
                          fast_stem_curve = TRUE,
                          model = "TF24_floor",
                          route = NULL) {
  method <- match.arg(method)
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }
  if (!(is.numeric(step) && length(step) == 1L && step > 0)) {
    stop("`step` must be a single positive number", call. = FALSE)
  }
  .gradient_check_settings(stationarity_tol, fast_stem_curve)
  # A caller who names `step` means it for the finite-difference arm too;
  # otherwise that arm takes the step matched to the solver it differences.
  fd_step <- if (missing(step)) .gradient_fd_step else step
  psi <- .gradient_check_psi(psi, dpsi_dtheta, method)

  # --- reusing a Leaf --------------------------------------------------------
  # `x` is a vessel: `traits` says where the gradient is taken and is applied to
  # it, so `traits` is required. `control` and `supply` are fixed in `x` and are
  # refused; `supply` is rebuilt from it because only its `kind` is read.
  if (!is.null(x)) {
    if (!inherits(x, "Leaf")) {
      stop("`x` must be a Leaf, from leaf_model(); got ", class(x)[[1]],
           call. = FALSE)
    }
    if (missing(traits)) {
      stop("`traits` must be given with `x`: it is the point the gradient is ",
           "taken at, and it is applied to `x` rather than read from it. Pass ",
           "the same leaf_traits() the leaf was built with.", call. = FALSE)
    }
    if (!missing(control)) {
      stop("`control` comes from `x` and cannot be given alongside it -- the ",
           "tolerances and the integrator were fixed when `x` was built.",
           call. = FALSE)
    }
    if (!missing(supply)) {
      stop("`supply` comes from `x` and cannot be given alongside it. `x` is on ",
           "the ", x$supply_kind, " path.", call. = FALSE)
    }
    supply <- if (identical(x$supply_kind, "single")) {
      leaf_supply_singlelayer()
    } else {
      leaf_supply_multilayer()
    }
  }
  single <- identical(supply$kind, "single")
  m <- .resolve_model(model, route, supply$kind)
  # One observation: every driver but the soil profile is one number.
  .check_scalars(list(PPFD = PPFD, leaf_specific_conductance_max =
                        leaf_specific_conductance_max, atm_vpd = atm_vpd,
                      ca = ca, leaf_temp = leaf_temp, atm_o2_kpa = atm_o2_kpa,
                      atm_kpa = atm_kpa), "leaf_gradient")
  # The thermal-damage layer's inputs are parameters ONLY on a leaf with its gate
  # on, which only a reused leaf can have: one built here has it off.
  damage <- !is.null(x) && isTRUE(x$use_thermal_damage_)

  # Resolved here rather than in set_drivers() because `resistance` can be
  # differentiated, so theta must hold the value the solve uses. The same
  # function set_drivers() calls, so there is one default.
  if (single && is.null(root_network)) {
    root_network <- .default_series_resistance()
  }

  # Model-owned fields are read off the leaf when there is one, so theta holds
  # the value the solve uses; NA on a model that does not own the field, so it
  # is left unwritten (see .owned_value()).
  owned <- unname(.owned_fields())
  theta <- .gradient_theta(traits, leaf_specific_conductance_max, supply,
                           root_network,
                           owned = stats::setNames(vapply(owned, function(f) {
                             .owned_value(model, f, x)
                           }, numeric(1)), owned),
                           damage = if (damage) .gradient_damage_values(x))
  if (is.null(pars)) {
    pars <- .gradient_available_pars(single, model, damage, m$route)
  }
  # Shared with leaf_gradient_batch(), so the two entry points cannot disagree
  # about which parameters exist or explain a rejection differently.
  .gradient_check_pars(pars, single, model, damage, m$route)
  # Checked here, unconditionally: the composite reads it only at a usable point,
  # so checking it there would accept a bad value at some psi and not others.
  dpsi_dtheta <- .gradient_dpsi_dtheta(dpsi_dtheta, pars, !is.null(psi))

  # One observation of the batch, on the caller's leaf when there is one:
  # construction is the largest cost on this surface (#52).
  l <- if (is.null(x)) leaf_model(traits, control, supply) else x
  resolved <- .resolve_drivers(l, psi_soil, PPFD, soil_depth, root_network,
                               leaf_specific_conductance_max, atm_vpd, ca,
                               leaf_temp, atm_o2_kpa, atm_kpa)
  par_names <- .gradient_par_names()
  row <- stats::setNames(rep(NA_real_, length(par_names)), par_names)
  row[names(theta)] <- theta
  if (!is.null(x)) {
    # The batch leaves the leaf at the base point's parameters; solved there,
    # on the model it seated, it is not left at the last perturbation (hazard 8).
    # `on.exit`, because an error row is rethrown below.
    on.exit(try(l$optimise(), silent = TRUE), add = TRUE)
  }
  res <- .gradient_run(l, .gradient_prepare(list(resolved)),
                       matrix(unname(row), nrow = 1L), pars, step, fd_step,
                       stationarity_tol, method, fast_stem_curve, psi,
                       if (is.null(dpsi_dtheta)) NULL else
                         matrix(dpsi_dtheta, nrow = 1L),
                       m)
  # One observation is one call: a row's failure is the call's, classified as
  # it would have been thrown.
  if (identical(res$status[[1L]], "error")) {
    with_phylloptim_conditions(stop(res$message[[1L]], call. = FALSE))
  }

  outs <- .gradient_output_names()
  list(gradient = matrix(res$gradient[1L, , , drop = FALSE],
                         nrow = length(pars), ncol = length(outs),
                         dimnames = list(pars, outs)),
       value = res$value[1L, ],
       method = res$method[[1L]],
       status = res$status[[1L]],
       H = res$H[[1L]],
       stationarity = res$stationarity[[1L]],
       model = model,
       route = m$route,
       M = stats::setNames(res$M[1L, ], pars),
       dY_dpsi = res$dY_dpsi[1L, ],
       psi = res$psi[[1L]])
}

# --- the model: which optimum is being differentiated --------------------------
#
# A (curve, route) pair, resolved by `.resolve_model()`. On the collar route the
# decision variable is the collar; on a leaf route it is `psi_leaf`, `collar` is
# not solved for, and its NA gradient is the truthful answer (hazard 11).
# Every model name: the cost curves, each of which takes either route. The one
# list `.gradient_model_pars()` and `.resolve_model()` validate against.
.gradient_models <- function() cost_curve_names()

# The finite-difference arm's step when the caller names none. Both routes
# root-find dprofit == 0, so psi* is stationary to ~1e-15 and a 1e-06 parameter
# step is resolvable; the floor moves with the solver, so re-sweep it before
# changing this or the solver (hazard 15). test-gradient-referee.R sweeps it.
.gradient_fd_step <- 1e-6

# Observations' resolved drivers, held C++-side: what `leaf_batch()` builds and
# `leaf_gradient()` builds for its one.
.gradient_prepare <- function(resolved) {
  pick <- function(field) vapply(resolved, function(a) a[[field]], numeric(1))
  gradient_batch_prepare(
    root_network = lapply(resolved, `[[`, "root_network"),
    psi_soil = lapply(resolved, `[[`, "psi_soil"),
    soil_depth = lapply(resolved, `[[`, "soil_depth"),
    PPFD = pick("PPFD"), atm_vpd = pick("atm_vpd"), ca = pick("ca"),
    leaf_temp = pick("leaf_temp"), atm_o2_kpa = pick("atm_o2_kpa"),
    atm_kpa = pick("atm_kpa"))
}

# The one crossing: every gradient, of one observation or many, is
# `gradient::batch` in C++. `m` is `.resolve_model()`'s pair.
.gradient_run <- function(l, drivers, theta, pars, step, fd_step,
                          stationarity_tol, method, fast_stem_curve, psi,
                          dpsi_dtheta, m) {
  res <- with_phylloptim_conditions(gradient_batch_run(
    l, drivers, theta, match(pars, .gradient_par_names()) - 1L, step,
    stationarity_tol, method, fast_stem_curve, psi, dpsi_dtheta,
    match(m$curve, cost_curve_names()) - 1L,
    identical(m$route, "collar"), fd_step))
  dimnames(res$gradient) <- list(NULL, pars, .gradient_output_names())
  dimnames(res$value) <- list(NULL, .gradient_output_names())
  dimnames(res$M) <- list(NULL, pars)
  dimnames(res$dY_dpsi) <- list(NULL, .gradient_output_names())
  res
}

# --- internals ---------------------------------------------------------------

# The differentiated outputs, read from C++ rather than restated here, so adding
# one is one edit in gradient.hpp. `collar` is psi* itself, so dcollar/dtheta is
# dpsi*/dtheta. Cached at first call: calling into the shared library while the
# package is being sourced is not safe to rely on.
.gradient_output_names <- local({
  nms <- NULL
  function() {
    if (is.null(nms)) {
      nms <<- gradient_output_names()
    }
    nms
  }
})

# The differentiable parameters, in the order C++ indexes them: leaf_traits()'
# arguments, then the non-traits. C++ derives its own (`gradient::par_names`)
# and test-gradient-batch.R compares the two: R passes positions, so appending
# is safe and reordering differentiates the wrong parameter (hazard 13).
# Computed at first call, for the collation reason above.
.gradient_par_names <- local({
  nms <- NULL
  function() {
    if (is.null(nms)) {
      nms <<- c(names(formals(leaf_traits)), .gradient_non_traits())
    }
    nms
  }
})

# The differentiable parameters that are not traits, in enumeration order after
# the traits: the drivers that are parameters, the fields the curves own (their
# definitions', in curve order), and the thermal-damage inputs.
.gradient_non_traits <- function() {
  c("leaf_specific_conductance_max", "resistance", unname(.owned_fields()),
    .gradient_damage_pars)
}

# The thermal-damage layer's inputs that reach the solve, in enumeration order:
# Leaf fields, read and written by name, available only with
# `use_thermal_damage_` on. Its dynamics parameters move no solved output and are
# absent (their gradient is leaf_damage_gradient()'s), except repair's two, which
# the repair-flux cost reads.
.gradient_damage_pars <- c("A_", "I_", "D_", "topt_offset_", "dTopt_max_", "K_A_",
                           "k_i_", "k_i_baseline_", "k_rec_", "k_rec_Ea_",
                           "c_acclim_maint_", "c_repair_maint_",
                           "c_protect_maint_", "c_repair_flux_")

# Their current values, off a leaf, in that order.
.gradient_damage_values <- function(l) {
  vapply(.gradient_damage_pars, function(nm) l[[nm]], numeric(1))
}

# What `pars` may name, shared by leaf_gradient() and leaf_gradient_batch(). A
# parameter is unavailable on one of two axes, the supply path and route
# (`resistance` is the single-potential supply's soil-to-collar resistance, which
# the leaf route never reads) or the model (an owned parameter belongs to one
# curve), and the refusal names which. Refusing is not a zero: a zero means "in
# the objective but inactive here", which is information.
#
# A model-owned parameter belongs to one curve, and `.owned_fields()`, read from
# the curves' definitions, says which (hazard 13). Not all are prices:
# `LeastCost_beta` is a dimensionless cost ratio.
.gradient_model_pars <- function(model) {
  if (!(model %in% .gradient_models())) {
    return(NULL)
  }
  owned <- .owned_fields()
  unname(owned[names(owned) != model])
}

.gradient_available_pars <- function(single, model = "TF24_floor", damage = FALSE,
                                     route = "leaf") {
  nms <- .gradient_par_names()
  if (!single || !identical(route, "collar")) nms <- setdiff(nms, "resistance")
  if (!damage) nms <- setdiff(nms, .gradient_damage_pars)
  excluded <- .gradient_model_pars(model)
  if (is.null(excluded)) {
    stop("no gradient route for model \"", model, "\". Available: ",
         paste(.gradient_models(), collapse = ", "), ".", call. = FALSE)
  }
  setdiff(nms, excluded)
}

.gradient_check_pars <- function(pars, single, model = "TF24_floor", damage = FALSE,
                                 route = "leaf") {
  unknown <- setdiff(pars, .gradient_available_pars(single, model, damage, route))
  if (length(unknown)) {
    off <- intersect(unknown, .gradient_damage_pars)
    if (length(off) && !damage) {
      stop("`pars` names `", off[[1L]], "`, an input of the thermal-damage layer, ",
           "which is off on this leaf. Set use_thermal_damage_ to TRUE on the leaf ",
           "you pass (`x`, or the batch's leaf) to differentiate it.",
           call. = FALSE)
    }
    why <- if (!single && "resistance" %in% unknown) {
      paste(" `resistance` is a parameter of the single-potential path only;",
            "on the multi-layer path the soil-to-collar resistances are derived",
            "from root carbon.")
    } else if ("resistance" %in% unknown) {
      paste(" `resistance` reaches only the collar route: the leaf route pins",
            "the potential upstream of the stem at `psi_soil`, so it never reads",
            "the soil-to-collar resistance. Pass route = \"collar\" to",
            "differentiate it.")
    } else {
      ""
    }
    owned <- intersect(unknown, .owned_fields())
    if (length(owned)) {
      par <- owned[[1L]]
      owner <- names(.owned_fields())[match(par, .owned_fields())]
      stop("`pars` names `", par, "`, which is not a parameter of the ", model,
           " model: it is the ", owner, " curve's PRESCRIBED price of water, ",
           "and on every other cost curve the price is emergent -- derived from ",
           "that curve's own parameters rather than set. Differentiate those ",
           "instead, or pass model = \"", owner, "\" (with route = ",
           "\"collar\" on the multi-layer supply).", call. = FALSE)
    }
    stop("`pars` names things this cannot differentiate: ",
         paste(unknown, collapse = ", "),
         ". Available: the traits from leaf_traits(), plus ",
         "`leaf_specific_conductance_max`",
         if (single && identical(route, "collar")) " and `resistance`" else "",
         ".", why, call. = FALSE)
  }
  invisible(pars)
}

# Everything this can differentiate, with its current value, as one named vector:
# the traits, then the non-traits. `leaf_specific_conductance_max` and
# `resistance` are drivers (set_drivers()); `resistance` is read from
# `r_R_V_sum[1]` and exists only on the single-potential path, where it is one
# scalar. The derivation does not care whether a parameter is a trait; only the
# setter differs (`gradient::apply()`).
.gradient_theta <- function(traits, kmax, supply, root_network,
                            owned = NULL, damage = NULL) {
  theta <- c(unlist(traits), leaf_specific_conductance_max = kmax)
  if (identical(supply$kind, "single")) {
    theta <- c(theta, resistance = root_network$r_R_V_sum[[1]])
  }
  # Each owned field, named by field: NA unless the caller is on the route that
  # owns it. C++ reads every slot, so each holds the value the solve uses; off
  # its route the value is never read.
  fields <- unname(.owned_fields())
  theta <- c(theta, stats::setNames(
    if (is.null(owned)) rep(NA_real_, length(fields)) else owned[fields], fields))
  # The thermal-damage inputs, read off the leaf when its gate is on and NA
  # otherwise -- in which case neither side reads them.
  if (is.null(damage)) {
    damage <- stats::setNames(rep(NA_real_, length(.gradient_damage_pars)),
                              .gradient_damage_pars)
  }
  c(theta, damage)
}


# The two settings C++ reads as a double and a bool, checked here: an NA
# tolerance would class every point interior, and an NA bool reads as TRUE.
.gradient_check_settings <- function(stationarity_tol, fast_stem_curve) {
  if (!(is.numeric(stationarity_tol) && length(stationarity_tol) == 1L &&
        !is.na(stationarity_tol) && stationarity_tol >= 0)) {
    stop("`stationarity_tol` must be a single non-negative number",
         call. = FALSE)
  }
  if (!(is.logical(fast_stem_curve) && length(fast_stem_curve) == 1L &&
        !is.na(fast_stem_curve))) {
    stop("`fast_stem_curve` must be TRUE or FALSE", call. = FALSE)
  }
  invisible(TRUE)
}

# `psi` and `dpsi_dtheta` only mean anything together, and neither means anything
# alongside `method`.
.gradient_check_psi <- function(psi, dpsi_dtheta, method) {
  if (is.null(psi)) {
    if (!is.null(dpsi_dtheta)) {
      stop("`dpsi_dtheta` is the trait response of a collar potential YOU ",
           "imposed, so it needs `psi`. Without one the collar is solved for ",
           "and its response is derived by the implicit function theorem.",
           call. = FALSE)
    }
    return(NULL)
  }
  if (!(is.numeric(psi) && length(psi) == 1L && is.finite(psi) && psi > 0)) {
    stop("`psi` must be a single finite positive number: it is a collar water ",
         "potential in MPa, as a positive magnitude (see the sign convention ",
         "in ?leaf_model).", call. = FALSE)
  }
  if (!identical(method, "auto")) {
    stop("`method` cannot be given with `psi`. The two routes it chooses ",
         "between are about a SOLVED operating point -- \"fd\" differences the ",
         "solve, which would answer about the optimum rather than about the ",
         "collar you gave -- so at a prescribed psi neither applies and the ",
         "reported method is \"prescribed\".", call. = FALSE)
  }
  # Coerced, and the caller must use the return: C++ takes a double, as
  # `.gradient_check_psi_batch()` hands it.
  as.numeric(psi)
}

# dpsi/dtheta as one value per parameter, in `pars` order. NULL on the solving
# path means "derive it"; NULL on the prescribed path means "the collar does not
# move with theta", which is zero and not the same statement.
.gradient_dpsi_dtheta <- function(dpsi_dtheta, pars, prescribed) {
  if (!prescribed) {
    return(NULL)
  }
  if (is.null(dpsi_dtheta)) {
    return(stats::setNames(numeric(length(pars)), pars))
  }
  # Finite, not just non-NA: an infinite value makes four columns +-Inf while
  # `profit`, assigned by the envelope theorem, stays finite and plausible.
  if (!is.numeric(dpsi_dtheta) || !all(is.finite(dpsi_dtheta))) {
    stop("`dpsi_dtheta` must be numeric and finite -- no NA, NaN or Inf. It ",
         "is dpsi/dtheta for a collar YOU imposed; a non-finite value returns ",
         "an all-infinite gradient row with a plausible `profit` in it.",
         call. = FALSE)
  }
  # Names are matched before recycling, or a named value would be applied to
  # whichever parameter recycling put it on.
  if (!is.null(names(dpsi_dtheta))) {
    if (!setequal(names(dpsi_dtheta), pars)) {
      stop("`dpsi_dtheta` is named, so its names must be exactly `pars`. ",
           "Missing: ", paste(setdiff(pars, names(dpsi_dtheta)),
                              collapse = ", "),
           "; unexpected: ", paste(setdiff(names(dpsi_dtheta), pars),
                                   collapse = ", "), call. = FALSE)
    }
    return(dpsi_dtheta[pars])
  }
  if (length(dpsi_dtheta) == 1L) {
    return(stats::setNames(rep(dpsi_dtheta, length(pars)), pars))
  }
  if (length(dpsi_dtheta) != length(pars)) {
    stop("`dpsi_dtheta` must be length 1 or one value per parameter (",
         length(pars), " for this `pars`); got ", length(dpsi_dtheta),
         call. = FALSE)
  }
  stats::setNames(dpsi_dtheta, pars)
}
