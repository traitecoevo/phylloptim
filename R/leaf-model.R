# The friendly R surface (issue #5, stage 2).
#
# Everything here is hand-written and sits ABOVE the generated glue in
# R/RcppR6.R, which is what makes it rewritable without touching the build, the
# YAML or CI. Nothing in this file reaches into C++ except through the generated
# `Leaf` R6 class.
#
# Two things it is for:
#
#   1. Separating traits from tolerances. The C++ constructor takes 19 positional
#      arguments and four of them -- GSS_tol_abs, ci_abs_tol, ci_niter,
#      vulnerability_curve_ncontrol -- are numerical settings sitting among the
#      physiology. A trait-calibration loop should not have to know which of 19
#      arguments are not traits. `leaf_traits()` and `leaf_control()` split them.
#   2. A one-call entry point. `leaf_solve()` takes drivers and returns a
#      data.frame, which is the shape someone reaching for plantecophys expects.

# The C++ default constructor's values, in one place. These are Eucalyptus
# saligna, and they are the defaults for every function below.
#
# ⚠️ These must equal leaf::Leaf's default constructor, and it is not enough to
# believe so: tests/testthat/test-golden.R solves through THESE defaults and
# compares bit-exactly against a golden file generated from a default-constructed
# C++ Leaf. If the two ever drift apart, that file fails.
#
# root_c / root_b / root_psi_crit / beta_R_H / beta_R_V default in MultiLayerRoots
# rather than in Leaf, because a second copy of the root Weibull pair inside Leaf
# is exactly hazard 1 in the developer guide -- it has already cost this project a
# wrong exponent in a manuscript draft. They are restated here because the R
# constructor has to pass them, not because Leaf owns them.
.leaf_trait_defaults <- list(
  vcmax_25 = 96,
  stem_c = 2.680147,
  stem_b = 3.898245,
  psi_crit = 5.870283,
  root_c = 2.680147,
  root_b = 3.898245,
  root_psi_crit = 5.870283,
  beta2 = 1.5,
  jmax_25 = 157.44,
  a = 0.30,
  curv_fact_elec_trans = 0.7,
  curv_fact_colim = 0.99,
  cost_scale_TF24 = 7.5,
  beta_R_H = 3.4e2,
  beta_R_V = 9.4e3
)

.leaf_control_defaults <- list(
  GSS_tol_abs = 1e-3,
  vulnerability_curve_ncontrol = 100,
  ci_abs_tol = 1e-3,
  ci_niter = 1000,
  integration_rule = 21,
  integration_tol = 1e-8
)

##' Physiological traits for a leaf
##'
##' The traits, and only the traits. Numerical settings live in
##' [leaf_control()]; the C++ constructor mixes the two and this is the split.
##'
##' Defaults are *Eucalyptus saligna*, and are the same values
##' `leaf::Leaf`'s default constructor uses.
##'
##' @section Two vulnerability curves, not one:
##' `stem_b`/`stem_c` describe the STEM Weibull curve, which drives the
##' hydraulic cost; `root_b`/`root_c` describe the ROOT curve, which drives
##' uptake. They are separate parameters with separate meanings and they are
##' named accordingly, because they were once an unmarked `b`/`c` pair alongside
##' `root_b`/`root_c` and an analysis used the root parameters for the stem cost.
##'
##' @param vcmax_25 maximum carboxylation rate at 25 C (umol m^-2 s^-1)
##' @param stem_c shape parameter of the stem vulnerability curve (unitless)
##' @param stem_b sensitivity parameter of the stem vulnerability curve (MPa)
##' @param psi_crit critical stem water potential (MPa, positive magnitude)
##' @param root_c shape parameter of the root vulnerability curve (unitless)
##' @param root_b sensitivity parameter of the root vulnerability curve (MPa)
##' @param root_psi_crit critical root water potential (MPa, positive magnitude)
##' @param beta2 exponent for the effect of hydraulic risk (unitless)
##' @param jmax_25 maximum electron transport rate at 25 C (umol m^-2 s^-1)
##' @param a quantum yield of photosynthetic electron transport (mol mol^-1)
##' @param curv_fact_elec_trans curvature of the light response curve (unitless)
##' @param curv_fact_colim curvature of the colimited photosynthesis equation
##' @param cost_scale_TF24 cost parameter for the TF24 profit model
##'   (umol m^-2 s^-1)
##' @param beta_R_H proportionality constant between minimum horizontal
##'   (intralayer) root hydraulic resistance and C_r^-1
##'   (MPa s mol C / mol H2O)
##' @param beta_R_V proportionality constant between minimum vertical
##'   (interlayer) root hydraulic resistance and dz^2/C_r
##'   (MPa mol C s / mol H2O / m^2)
##'
##' @return A `leaf_traits` object; a named list.
##' @seealso [leaf_control()], [leaf_model()], [leaf_solve()]
##' @examples
##' leaf_traits()
##' leaf_traits(vcmax_25 = 120, stem_b = 2.5)
##' @export
leaf_traits <- function(vcmax_25 = 96,
                        stem_c = 2.680147,
                        stem_b = 3.898245,
                        psi_crit = 5.870283,
                        root_c = 2.680147,
                        root_b = 3.898245,
                        root_psi_crit = 5.870283,
                        beta2 = 1.5,
                        jmax_25 = 157.44,
                        a = 0.30,
                        curv_fact_elec_trans = 0.7,
                        curv_fact_colim = 0.99,
                        cost_scale_TF24 = 7.5,
                        beta_R_H = 3.4e2,
                        beta_R_V = 9.4e3) {
  out <- list(vcmax_25 = vcmax_25, stem_c = stem_c, stem_b = stem_b,
              psi_crit = psi_crit, root_c = root_c, root_b = root_b,
              root_psi_crit = root_psi_crit, beta2 = beta2,
              jmax_25 = jmax_25, a = a,
              curv_fact_elec_trans = curv_fact_elec_trans,
              curv_fact_colim = curv_fact_colim,
              cost_scale_TF24 = cost_scale_TF24,
              beta_R_H = beta_R_H, beta_R_V = beta_R_V)
  .check_scalars(out, "leaf_traits")
  structure(out, class = c("leaf_traits", "list"))
}

##' Numerical settings for a leaf solve
##'
##' Tolerances and iteration counts. These are not traits: they control how
##' precisely the solve is carried out, not what is being solved, and separating
##' them from [leaf_traits()] is the point of having two functions.
##'
##' @section What is deliberately not here:
##' The internal leaf-temperature clamp (`-40 C` to `70 C`) is a compile-time
##' constant and stays one. It is a guard that keeps an extreme non-equilibrium
##' transpiration on the energy-balance path from driving the Arrhenius block
##' non-finite -- not a tolerance anyone tunes. Making it settable would let a
##' caller produce non-finite photosynthetic parameters and get NaNs back with no
##' indication of why.
##'
##' @param GSS_tol_abs absolute tolerance for the golden-section search over
##'   stem water potential. The nested solvers amplify perturbations up to about
##'   this size, so it also sets how well determined the reported operating
##'   point is.
##' @param vulnerability_curve_ncontrol number of control points used to
##'   pre-integrate the two Weibull vulnerability curves into splines. Higher is
##'   more accurate and slower to construct; it does not affect solve speed.
##' @param ci_abs_tol absolute tolerance for root-solving intercellular CO2
##' @param ci_niter maximum iterations for root-solving intercellular CO2
##' @param integration_rule accepted and ignored. It selected a Gauss-Kronrod
##'   rule order in plant's compiled quadrature; the header-only adaptive
##'   Simpson that replaced it has no rule order. Retained because plant's R
##'   bindings expose it.
##' @param integration_tol absolute tolerance for the direct quadrature in
##'   `transpiration_full_integration()`. Not on the production path, which
##'   reads the pre-integrated spline instead.
##'
##' @return A `leaf_control` object; a named list.
##' @seealso [leaf_traits()], [leaf_model()]
##' @examples
##' leaf_control()
##' leaf_control(GSS_tol_abs = 1e-5)
##' @export
leaf_control <- function(GSS_tol_abs = 1e-3,
                         vulnerability_curve_ncontrol = 100,
                         ci_abs_tol = 1e-3,
                         ci_niter = 1000,
                         integration_rule = 21,
                         integration_tol = 1e-8) {
  out <- list(GSS_tol_abs = GSS_tol_abs,
              vulnerability_curve_ncontrol = vulnerability_curve_ncontrol,
              ci_abs_tol = ci_abs_tol, ci_niter = ci_niter,
              integration_rule = integration_rule,
              integration_tol = integration_tol)
  .check_scalars(out, "leaf_control")
  structure(out, class = c("leaf_control", "list"))
}

##' Choose how water reaches the root collar
##'
##' Two supply paths, and the choice belongs with the model rather than with the
##' drivers, so it is made once in [leaf_model()].
##'
##' `leaf_supply_multilayer()` is the default and is what plant uses: a soil
##' profile of one or more layers, each with its own water potential, feeding a
##' root network whose resistances are derived from a root carbon profile.
##'
##' `leaf_supply_single()` collapses the whole soil-to-collar path to **one
##' resistance**. This is what a leaf physiologist arriving from `plantecophys`
##' or `tealeaves` actually has: a soil water potential and no root-mass profile.
##' It is also what makes comparison against other optimality models meaningful,
##' since Medlyn, Prentice least-cost and Cowan-Farquhar are all formulated
##' against a single soil potential — comparing them against a multi-layer root
##' network compares two things at once.
##'
##' @section Why this is not a settable field:
##' The obvious R interface would be `leaf$supply_kind <- "single"`. It is not
##' offered, because flipping the tag on its own leaves the other path's state
##' configured and silently ignored, and flipping back makes it stale rather than
##' absent. Both entry points here reconfigure the object completely, so the tag
##' and the supply state can never disagree. The cost is that switching clears
##' the solved state and the drivers must be set again — which is honest, since
##' the two paths read different inputs.
##'
##' @param resistance the whole soil-to-collar path as one series resistance,
##'   **per unit leaf area**, MPa s (mol H2O)^-1 m^2 leaf. Must be positive: a
##'   zero resistance is an infinite flux. Per unit leaf area because the leaf is
##'   purely intensive — nothing here may scale with plant size.
##' @param gravity_head head required to lift water to the collar, MPa. Zero by
##'   default, which is right for a bare leaf that is not thinking about rooting
##'   depth.
##'
##' @return A `leaf_supply` object, for the `supply` argument of [leaf_model()]
##'   and [leaf_solve()].
##' @examples
##' # a bare leaf: one soil potential, one resistance, no root profile
##' leaf_solve(psi_soil = 1.5, PPFD = 900,
##'            supply = leaf_supply_single(resistance = 1e3))
##' @export
leaf_supply_single <- function(resistance, gravity_head = 0) {
  out <- list(kind = "single", resistance = resistance,
              gravity_head = gravity_head)
  .check_scalars(out[c("resistance", "gravity_head")], "leaf_supply_single")
  if (resistance <= 0) {
    stop("leaf_supply_single(): `resistance` must be positive; a zero ",
         "resistance is an infinite flux", call. = FALSE)
  }
  if (gravity_head < 0) {
    stop("leaf_supply_single(): `gravity_head` must be non-negative (MPa)",
         call. = FALSE)
  }
  structure(out, class = c("leaf_supply", "list"))
}

##' @rdname leaf_supply_single
##' @export
leaf_supply_multilayer <- function() {
  structure(list(kind = "multilayer"), class = c("leaf_supply", "list"))
}

##' Build a leaf
##'
##' The recommended way to construct a leaf. [Leaf()] is the raw C++ constructor,
##' with all nineteen arguments positional and no defaults; this splits them into
##' traits and numerical settings, defaults both, and initialises the integrator.
##'
##' The result is a stateful R6 object: set drivers with [set_drivers()], solve
##' with `$find_root_collar_psi()`, then read the operating point off the object
##' or with [operating_point()]. For a single solve you probably want
##' [leaf_solve()] instead, which does all of that in one call.
##'
##' @param traits a [leaf_traits()] object
##' @param control a [leaf_control()] object
##' @param supply how water reaches the root collar: [leaf_supply_multilayer()]
##'   (the default) or [leaf_supply_single()]
##'
##' @return A `Leaf` R6 object.
##' @seealso [leaf_traits()], [leaf_control()], [leaf_supply_single()],
##'   [set_drivers()], [leaf_solve()]
##' @examples
##' l <- leaf_model()
##' set_drivers(l, psi_soil = 2.0, PPFD = 900)
##' l$find_root_collar_psi()
##' operating_point(l)
##' @export
leaf_model <- function(traits = leaf_traits(), control = leaf_control(),
                       supply = leaf_supply_multilayer()) {
  if (!inherits(traits, "leaf_traits")) {
    stop("`traits` must come from leaf_traits()", call. = FALSE)
  }
  if (!inherits(control, "leaf_control")) {
    stop("`control` must come from leaf_control()", call. = FALSE)
  }
  if (!inherits(supply, "leaf_supply")) {
    stop("`supply` must come from leaf_supply_multilayer() or ",
         "leaf_supply_single()", call. = FALSE)
  }

  # Positional, because that is what the generated constructor takes. The whole
  # value of this function is that the ordering is written down once, here,
  # instead of at every call site.
  l <- Leaf(
    vcmax_25 = traits$vcmax_25,
    stem_c = traits$stem_c,
    stem_b = traits$stem_b,
    psi_crit = traits$psi_crit,
    root_c = traits$root_c,
    root_b = traits$root_b,
    root_psi_crit = traits$root_psi_crit,
    beta2 = traits$beta2,
    jmax_25 = traits$jmax_25,
    a = traits$a,
    curv_fact_elec_trans = traits$curv_fact_elec_trans,
    curv_fact_colim = traits$curv_fact_colim,
    GSS_tol_abs = control$GSS_tol_abs,
    vulnerability_curve_ncontrol = control$vulnerability_curve_ncontrol,
    ci_abs_tol = control$ci_abs_tol,
    ci_niter = control$ci_niter,
    cost_scale_TF24 = traits$cost_scale_TF24,
    beta_R_H = traits$beta_R_H,
    beta_R_V = traits$beta_R_V
  )
  l$initialize_integrator(control$integration_rule, control$integration_tol)
  # After the integrator, because set_supply_single clears the solved state --
  # not the integrator tolerance, but relying on that ordering would be a
  # dependency on an implementation detail of setup_clean_leaf.
  if (identical(supply$kind, "single")) {
    l$set_supply_single(supply$resistance, supply$gravity_head)
  }
  l
}

##' Set the drivers for a leaf
##'
##' Named and defaulted, over the C++ `set_physiology()`'s ten positional
##' arguments. Modifies `x` in place and returns it invisibly, so it can be
##' piped.
##'
##' @section Water potentials are positive magnitudes:
##' `psi_soil` is a positive number of MPa. There is one representation for water
##' potential throughout this package and it is asserted at the input boundary
##' rather than documented: a negative entry is an error. That check is the only
##' thing standing between a script written against the old signed convention and
##' a plausible wrong number, so this wrapper surfaces it rather than smoothing
##' it over -- it does not take `abs()`.
##'
##' @section Multiple soil layers:
##' Pass vectors for `psi_soil`, and optionally for `soil_depth` and
##' `root_carbon_per_leaf_area`. All three must end up the same length. By
##' default `soil_depth` becomes 1 m layers and the root carbon is split evenly.
##'
##' @param x a `Leaf`, from [leaf_model()]
##' @param psi_soil soil water potential, MPa, **positive magnitude**. Length
##'   sets the number of soil layers.
##' @param PPFD photosynthetic photon flux density, umol m^-2 s^-1
##' @param soil_depth cumulative depth to the bottom of each layer, m. Defaults
##'   to 1 m layers.
##' @param root_carbon_per_leaf_area root carbon **per unit leaf area**, kg C
##'   m^-2 leaf, per layer. Defaults to 20 split evenly across layers. This is
##'   the one argument with no good default -- a bare-leaf user does not have a
##'   root carbon profile, which is why 20 is a stand-in rather than a
##'   recommendation. It is also why the single-potential supply path exists.
##' @param leaf_specific_conductance_max maximum leaf-specific hydraulic
##'   conductance, kg m^-2 s^-1 MPa^-1
##' @param atm_vpd atmospheric vapour pressure deficit, kPa
##' @param ca atmospheric CO2 partial pressure, Pa
##' @param leaf_temp leaf temperature, deg C
##' @param atm_o2_kpa atmospheric O2 partial pressure, kPa
##' @param atm_kpa atmospheric pressure, kPa. Note this is not decorative: the
##'   ppm-to-Pa conversion is derived from it, and getting it wrong moves
##'   assimilation by percent, not by rounding.
##'
##' @return `x`, invisibly.
##' @seealso [leaf_model()], [leaf_solve()]
##' @examples
##' l <- leaf_model()
##' set_drivers(l, psi_soil = 2.0, PPFD = 900, atm_vpd = 1.5)
##' # three drying layers
##' set_drivers(l, psi_soil = c(1.0, 1.5, 2.0), PPFD = 900)
##' @export
set_drivers <- function(x,
                        psi_soil,
                        PPFD = 900,
                        soil_depth = NULL,
                        root_carbon_per_leaf_area = NULL,
                        leaf_specific_conductance_max = 3.14e-5,
                        atm_vpd = 2.0,
                        ca = 40.0,
                        leaf_temp = 25.0,
                        atm_o2_kpa = 21.0,
                        atm_kpa = 101.3) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  psi_soil <- as.numeric(psi_soil)
  n <- length(psi_soil)
  if (n < 1L) {
    stop("`psi_soil` must have at least one layer", call. = FALSE)
  }

  # The C++ boundary rejects this too, and it is the enforcement. This check is
  # here only to name the fix, because "psi_soil must be positive magnitudes"
  # arriving from inside a solve is a less useful thing to read than this is.
  if (any(psi_soil < 0)) {
    stop("`psi_soil` must be positive magnitudes in MPa, not signed ",
         "potentials. If this came from an older script, drop the minus sign ",
         "rather than negating the result.", call. = FALSE)
  }

  single <- identical(x$supply_kind, "single")

  if (single) {
    # The single-potential path reads psi_soil[1] and nothing else: no depth
    # profile, no root carbon. Rather than quietly ignoring extra layers -- which
    # would let someone pass a profile and believe it was used -- say so.
    if (n != 1L) {
      stop("`psi_soil` must be a single value on the single-potential supply ",
           "path; got ", n, " layers. Build the leaf with ",
           "leaf_supply_multilayer() for a layered profile.", call. = FALSE)
    }
    if (!is.null(soil_depth) || !is.null(root_carbon_per_leaf_area)) {
      stop("`soil_depth` and `root_carbon_per_leaf_area` are not used on the ",
           "single-potential supply path -- the whole soil-to-collar path is ",
           "the `resistance` given to leaf_supply_single(). Drop them, or ",
           "build the leaf with leaf_supply_multilayer().", call. = FALSE)
    }
    # set_physiology keeps one signature for both paths and requires the three
    # vectors to match in length; on this path it reads only psi_soil[1], so
    # these two are placeholders rather than inputs.
    soil_depth <- 1.0
    root_carbon_per_leaf_area <- 0.0
  } else {
    soil_depth <- if (is.null(soil_depth)) {
      seq_len(n) * 1.0
    } else {
      as.numeric(soil_depth)
    }
    root_carbon_per_leaf_area <- if (is.null(root_carbon_per_leaf_area)) {
      rep(20 / n, n)
    } else {
      as.numeric(root_carbon_per_leaf_area)
    }

    if (length(soil_depth) != n || length(root_carbon_per_leaf_area) != n) {
      stop("`psi_soil` (", n, "), `soil_depth` (", length(soil_depth),
           ") and `root_carbon_per_leaf_area` (",
           length(root_carbon_per_leaf_area),
           ") must all have one entry per soil layer", call. = FALSE)
    }
  }

  x$set_physiology(
    root_carbon_per_leaf_area = root_carbon_per_leaf_area,
    PPFD = PPFD,
    psi_soil = psi_soil,
    soil_depth = soil_depth,
    leaf_specific_conductance_max = leaf_specific_conductance_max,
    atm_vpd = atm_vpd,
    ca = ca,
    leaf_temp = leaf_temp,
    atm_o2_kpa = atm_o2_kpa,
    atm_kpa = atm_kpa
  )
  invisible(x)
}

##' The solved operating point, as one row
##'
##' Reads the outputs off a leaf that has been solved. Call after
##' `$find_root_collar_psi()`; before that the values are the missing-value
##' sentinels the object was constructed with.
##'
##' @param x a solved `Leaf`
##' @return A one-row data.frame.
##' @seealso [leaf_solve()], which does the whole thing in one call.
##' @examples
##' l <- leaf_model()
##' set_drivers(l, psi_soil = 2.0, PPFD = 900)
##' l$find_root_collar_psi()
##' operating_point(l)
##' @export
operating_point <- function(x) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  consumption <- x$soil_consumption_
  data.frame(
    psi_stem = x$opt_psi_stem_,        # MPa, positive magnitude
    collar = x$opt_root_psi_,          # MPa, positive magnitude
    ci = x$ci_,                        # Pa
    A = x$assim_colimited_,            # umol C m^-2 s^-1
    E = x$transpiration_,              # kg H2O m^-2 s^-1
    gc = x$stom_cond_CO2_,             # mol CO2 m^-2 s^-1
    profit = x$profit_,                # umol C m^-2 s^-1
    hydraulic_cost = x$hydraulic_cost_,
    E_up = x$E_up_,
    uptake = sum(consumption[is.finite(consumption)]),
    lambda = x$lambda,                 # dA/dE
    g1_eff = x$g1_eff                  # the Medlyn g1 this leaf implies
  )
}

##' Solve a leaf, in one call
##'
##' Drivers in, operating point out. Vectorised over the scalar drivers, so a
##' response curve is one call rather than a loop.
##'
##' This is the entry point for someone who would otherwise reach for
##' `plantecophys::Photosyn()`. The difference in what comes back is the point of
##' the package: `gc` here is not from a fitted conductance model, it is what
##' falls out of maximising profit over an explicit hydraulic path, and `g1_eff`
##' reports the Medlyn `g1` that would have been needed to reproduce it.
##'
##' @section Recycling:
##' `PPFD`, `atm_vpd`, `ca`, `leaf_temp`, `atm_o2_kpa`, `atm_kpa` and
##' `leaf_specific_conductance_max` are recycled to a common length, the usual R
##' way. `psi_soil` is recycled too **when it is a plain numeric vector** --
##' each element is then a separate single-layer solve. For multi-layer solves
##' pass a list of numeric vectors, one per row; `soil_depth` and
##' `root_carbon_per_leaf_area` follow the same rule.
##'
##' @inheritParams set_drivers
##' @param traits a [leaf_traits()] object
##' @param control a [leaf_control()] object
##' @param supply how water reaches the root collar: [leaf_supply_multilayer()]
##'   (the default) or [leaf_supply_single()]. On the single-potential path
##'   `soil_depth` and `root_carbon_per_leaf_area` must be omitted, and each
##'   `psi_soil` is one value rather than a profile.
##' @param reuse solve every row with one `Leaf` object (`TRUE`) or construct a
##'   fresh one per row (`FALSE`, the default). Reuse is faster because the two
##'   vulnerability splines are built once, and it is safe -- every exit from the
##'   solve writes its own outputs. `FALSE` exists so that can be checked rather
##'   than assumed, and the package's tests do check it.
##'
##' @return A data.frame with one row per driver combination: the drivers, then
##'   the columns of [operating_point()].
##' @seealso [leaf_model()] for the stateful interface.
##' @examples
##' # one point
##' leaf_solve(psi_soil = 2.0, PPFD = 900)
##'
##' # a light response curve
##' leaf_solve(psi_soil = 1.0, PPFD = seq(100, 1800, length.out = 6))
##'
##' # a drought response, and the effective g1 it implies
##' d <- leaf_solve(psi_soil = seq(0.5, 5, length.out = 6), PPFD = 900)
##' d[, c("psi_soil", "A", "gc", "g1_eff")]
##'
##' # two three-layer profiles
##' leaf_solve(psi_soil = list(c(1, 1.5, 2), c(3, 3.5, 4)), PPFD = 900)
##' @export
leaf_solve <- function(psi_soil,
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
                       reuse = TRUE) {
  layered <- .as_layer_list(psi_soil, "psi_soil")
  scalars <- list(PPFD = PPFD, atm_vpd = atm_vpd, ca = ca,
                  leaf_temp = leaf_temp, atm_o2_kpa = atm_o2_kpa,
                  atm_kpa = atm_kpa,
                  leaf_specific_conductance_max = leaf_specific_conductance_max)

  n <- max(length(layered), vapply(scalars, length, integer(1)))
  if (n == 0L) {
    stop("nothing to solve: every driver has length zero", call. = FALSE)
  }
  layered <- .recycle_to(layered, n, "psi_soil")
  scalars <- lapply(seq_along(scalars), function(i) {
    .recycle_to(scalars[[i]], n, names(scalars)[[i]])
  })
  names(scalars) <- c("PPFD", "atm_vpd", "ca", "leaf_temp", "atm_o2_kpa",
                      "atm_kpa", "leaf_specific_conductance_max")

  depths <- if (is.null(soil_depth)) NULL else
    .recycle_to(.as_layer_list(soil_depth, "soil_depth"), n, "soil_depth")
  carbon <- if (is.null(root_carbon_per_leaf_area)) NULL else
    .recycle_to(.as_layer_list(root_carbon_per_leaf_area,
                               "root_carbon_per_leaf_area"),
                n, "root_carbon_per_leaf_area")

  shared <- if (reuse) leaf_model(traits, control, supply) else NULL

  rows <- lapply(seq_len(n), function(i) {
    l <- if (reuse) shared else leaf_model(traits, control, supply)
    set_drivers(l,
                psi_soil = layered[[i]],
                PPFD = scalars$PPFD[[i]],
                soil_depth = if (is.null(depths)) NULL else depths[[i]],
                root_carbon_per_leaf_area =
                  if (is.null(carbon)) NULL else carbon[[i]],
                leaf_specific_conductance_max =
                  scalars$leaf_specific_conductance_max[[i]],
                atm_vpd = scalars$atm_vpd[[i]],
                ca = scalars$ca[[i]],
                leaf_temp = scalars$leaf_temp[[i]],
                atm_o2_kpa = scalars$atm_o2_kpa[[i]],
                atm_kpa = scalars$atm_kpa[[i]])
    l$find_root_collar_psi()

    # One row of drivers, then one row of outputs. psi_soil is reported as the
    # wettest layer when there are several, with the layer count alongside, so
    # the column stays numeric and the frame stays rectangular.
    cbind(
      data.frame(psi_soil = min(layered[[i]]),
                 layers = length(layered[[i]]),
                 PPFD = scalars$PPFD[[i]],
                 atm_vpd = scalars$atm_vpd[[i]],
                 ca = scalars$ca[[i]],
                 leaf_temp = scalars$leaf_temp[[i]],
                 atm_kpa = scalars$atm_kpa[[i]]),
      operating_point(l)
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# --- internals ---------------------------------------------------------------

.check_scalars <- function(x, what) {
  bad <- names(x)[!vapply(x, function(v) {
    is.numeric(v) && length(v) == 1L && is.finite(v)
  }, logical(1))]
  if (length(bad)) {
    stop(what, "(): each argument must be a single finite number; ",
         "the following ", if (length(bad) == 1L) "is" else "are", " not: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Normalise a driver that may describe soil layers into a list of numeric
# vectors, one per row. A plain numeric vector means "one single-layer solve per
# element"; a list means "one multi-layer solve per element".
.as_layer_list <- function(x, what) {
  if (is.list(x)) {
    out <- lapply(x, as.numeric)
  } else {
    out <- as.list(as.numeric(x))
  }
  if (any(vapply(out, length, integer(1)) == 0L)) {
    stop("`", what, "` has an entry with no layers", call. = FALSE)
  }
  out
}

.recycle_to <- function(x, n, what) {
  len <- length(x)
  if (len == n) {
    return(x)
  }
  if (len == 1L) {
    return(rep(x, n))
  }
  stop("`", what, "` has length ", len,
       ", which cannot be recycled to ", n,
       ". Driver lengths must be 1 or the longest driver's length.",
       call. = FALSE)
}
