# leaf (development version)

## The AD replicas are gone, and the model's precision floor moved (#4, PLAN 11b)

Two small results-moving changes, landed in this order because the order matters:

**`psi_stem_to_ci`'s tolerance is 1e-10, was 1e-7** — moves results by
**1.82e-07**, lands **335× closer to a converged solve**, costs **+3.4%**. This
tolerance was invisible while the collar solve's own `GSS_tol_abs` (1e-3) dominated;
once PLAN 11a removed that, it became the model's dominant amplifier and therefore
the floor of what every reported output *means*. The working magnitudes for reading
a diff are now **~1e-16 reassociation, ~1e-9 solver floor, ~1e-4 a real
difference**. Not the same knob as the settable `ci_abs_tol` (default 1e-3), which
reaches only the off-path `optimise_psi_stem_*` solvers.

**`namespace detail`'s hand-maintained AD replicas are deleted** — moves results by
**3.51e-10**. `assim_colimited`, its two components, `hydraulic_cost_TF` and
`proportion_of_conductivity` are now `T = double` instantiations of scalar-generic
member templates, and `dprofit_droot_collar_psi` differentiates the *same code*. So
the forward model and its derivative can no longer be derivatives of different
functions — which they were: the deleted replica associated the electron-limited
term differently, making the AD derivative exact for a function the model did not
evaluate.

**The `double` path is bit-identical**, checked directly across all five entry
points. The whole 3.51e-10 is the derivative changing.

Doing the tolerance first is what makes the unification a 3.51e-10 change rather
than a 4.98e-07 one — a factor of 1400 for free, purely from ordering.

To be clear about what this second change is *not*: it does **not** make the solved
operating point more accurate. Both the old and new builds locate their respective
roots to ~1e-14, and the two derivatives differ so little that their roots differ by
~1e-15 MPa. It is a maintainability and correctness-of-construction fix.

## ⚠️ The collar solve changed, and results moved (#4, PLAN 11a)

**`find_root_collar_psi` no longer maximises profit by golden-section search. It
solves the first-order condition, `dprofit/dpsi == 0`, by a safeguarded
root-find.** Every solved operating point moves by roughly the argmax correction —
worst change over the 288-point golden grid **1.5e-03** relative, 240 rows, with
`profit` itself moving at most 1.4e-03 because it is the maximum and therefore
flat. Shut-down rows are unchanged. **If you have recorded numbers from this
package, they will differ.**

This is a correctness fix, and the evidence is the residual rather than the
values: `|dprofit|` at the returned collar improved on **240 of 240** feasible
grid rows with none worse, from a median of 7.8e-04 to **5.6e-15** on the 198 rows
with an interior optimum. The remaining 42 have a *constrained* optimum pinned to
a bracket bound, where the gradient is genuinely non-zero.

Why it was worth moving results:

- **Trait derivatives were unusable, and worse than noisy.** Golden section
  resolved the argmax only to `GSS_tol_abs`, so a finite difference in a
  photosynthetic trait returned exactly **zero** below a relative step of 1e-4 —
  silently dropping a whole term of `dA/dtheta`. For traits in the hydraulic path
  it was worse: the argmax came back smooth, plausible and **sign-inverted**
  (`root_b` gave −2.6e-03 where the truth is +2.6e-04). A gradient-based
  calibration would have walked those traits the wrong way with nothing to show
  it. Finite differences now agree to ~4 digits at any relative step from 1e-8 to
  1e-2, which unblocks #6.
- **It is faster: 24.5%**, 2.65 against 3.51 µs/solve, interleaved at reps=2000.
- **The argmax got smoother, not rougher** — ~1000× smaller second differences in
  a trait. That is the opposite of what the guide's hazard 3 would lead you to
  expect, and it is measured. ⚠️ In a *trait*; smoothness in plant state, which is
  what hazard 3 is actually about, needs re-measuring on the plant side once
  plant #591 clears.

`GSS_tol_abs` is still a control and still has two jobs — the "interval too narrow
to solve over" threshold and the off-path single-layer optimisers — but it no
longer sets how well the reported operating point is determined.

Also in this change: `dprofit_droot_collar_psi` gained an optional `bool*
feasible` out-parameter, so a caller searching for a zero can tell a real
stationary point from the `0.0` it returns on its shut-down exits. The default
leaves every existing caller, including plant's TF24f, untouched.

## The package is callable from R (#5, stage 1)

`leaf::Leaf` now has an R interface, so this is no longer a `LinkingTo`-only
package. `library(leaf)` gives you an R6 `Leaf` with the drivers, the solve and
the whole operating point; `vignette("leaf")` walks through a solve, a drought
response and a light response. **λ and `g1_eff` are exposed for the first time**
— they existed in C++ and were unreachable from R.

**The C++ headers are unchanged, and still need no R.** The model stays a set of
self-contained headers under `inst/include` that use no R and no Rcpp; the R
layer sits on top of them and is never included by them. So a `LinkingTo: leaf`
consumer sees nothing new, and the model remains linkable straight into a C++
program or a Python extension. There is now a **CMake package** (`leaf::leaf`)
for exactly that, and CI builds and installs it on runners with no R, so the
claim is tested rather than asserted. See PLAN item 6a.

Two things worth knowing if you are working on this package:

- **RcppR6 is not a declared dependency.** The generated glue is committed, so
  only a developer regenerating it needs the generator. CI regenerates and diffs
  to catch a stale commit. PLAN item 6c has the reasoning, including what would
  make it worth swapping for odelia's hand-written style.
- **The golden file's bit-exactness turns out to depend on the optimisation
  level, not only on the platform**: identical at `-O1`/`-O2`/`-O3`, off by
  3.47e-15 at `-O0`, which does not contract `a*b + c` into an FMA. A debug build
  failing `test_golden` by ~1e-15 has found nothing.

## A surface you can type at a console (#5, stage 2)

The bindings above are a faithful translation of the C++, which is what let them
be checked against the golden file — and it means nineteen positional arguments.
On top of them:

- **`leaf_solve()`** — drivers in, operating point out as a data.frame,
  vectorised, so a response curve is one call. This is the entry point for
  someone who would otherwise reach for `plantecophys::Photosyn()`.
- **`leaf_traits()` and `leaf_control()`** — the split the issue asked for. Four
  of the constructor's nineteen arguments are tolerances sitting among the
  physiology, and a trait-calibration loop should not have to know which. A test
  asserts the two functions partition the constructor exactly.
- **`leaf_model()`** — named and defaulted, and the recommended constructor.
  `Leaf()` remains exported as the raw one.
- **`set_drivers()`** and **`operating_point()`** for the stateful path.
- **`vignette("leaf")`** — a solve, a drought response, the light/VPD surface,
  and the profit function being maximised.

Two decisions worth knowing:

- **The `psi_soil >= 0` rejection is surfaced, not smoothed over.** A friendly
  wrapper is exactly where someone would be tempted to `abs()` it, and that check
  is the only thing between a script written against the old signed convention
  and a plausible wrong number.
- **The leaf-temperature clamp is NOT in `leaf_control()`**, though the issue
  listed it. It is a guard keeping the Arrhenius block finite on the
  energy-balance path, not a tolerance anyone tunes; making it settable would let
  a caller get NaNs back with no indication why.

## A bare leaf needs no root carbon profile (#5 stage 3, #32)

`leaf_supply_single()` collapses the whole soil-to-collar path to one series
resistance, so a leaf physiologist with a soil water potential and no root-mass
profile can use the model without going through a plant-shaped one to get at a
leaf. It is also what makes the optimality-model comparison meaningful, since
Medlyn, Prentice least-cost and Cowan-Farquhar are all written against a single
soil potential.

```r
leaf_solve(psi_soil = 1.5, PPFD = 900,
           supply = leaf_supply_single(resistance = 1e3))
```

The path is chosen when the leaf is built, and **there is no
`leaf$supply_kind <- "single"`.** Flipping a tag would leave the other path's
state configured and silently ignored — and flipping back would make it stale
rather than absent. `leaf_supply_multilayer()` and `leaf_supply_single()`
reconfigure the object completely instead, so it can never be in a state where
the tag and the supply disagree. `supply_kind`, `single_resistance_` and
`single_gravity_head_` are readable but not settable, for the same reason.

On the multi-layer path `root_carbon_per_leaf_area` still has no good default and
the R layer supplies a stand-in. Taking resistances there too is #33, which is
coupled with plant.

# leaf 0.1.0

**The first version bump since the package was created, and it exists so downstream can
pin.** `0.0.1` has been the version through four merged PRs that changed the API and the
results, so a `leaf (>= 0.0.1)` requirement was satisfied by the original extraction —
which `plant` cannot compile against at all. That is the same trap odelia 0.2.1 was cut
for: a version that names several different header sets is not something a consumer can
depend on.

A minor bump rather than a patch, because the C++ API broke in ways a consumer sees at
compile time, and the numbers moved.

## Breaking API changes

- **`set_physiology` takes 10 arguments, not 14** (#15). `area_leaf`, `rho`, `a_bio` and
  `sapwood_volume_per_leaf_area` were dead stores — assigned, never read — and are gone,
  along with the four members of the same names. `mass_root_prop` became
  `root_carbon_per_leaf_area`, and it is the old value **divided by `area_leaf`**: the
  model is purely intensive, and uptake is exactly homogeneous in that ratio. Passing the
  old absolute carbon compiles and runs, and gives a root system too weak by a factor of
  `1/area_leaf`.
- **`root_collar_psi_` is now `opt_root_psi_`, and it is a positive magnitude** (#25).
  Renamed rather than reused deliberately: keeping the name with a flipped sign is the one
  outcome where an old analysis reads the wrong value in silence, and a rename gives a
  binding error instead.
- **Every ψ in the package is a positive magnitude in MPa** (#25).
  `psi_soil_inverted_`, `psi_soil_inverted_vec_` and `supply_psi_soil_inverted()` are
  deleted. `E_from_Soil_to_Root_Collar`, `find_root_psi`, `find_psi_stem_from_psi_root`,
  `dE_from_soil_dpsi_collar` and `transpiration_to_psi_stem` keep their signatures but
  take magnitudes, and `find_root_psi`'s bracket ends swap (wettest layer first). They
  validate the soil vector and stop on a negative entry, so a pre-#25 caller fails loudly
  rather than returning a wrong number.
- **`dE_from_soil_dpsi_collar` returns a positive conductance** (#25), so callers that
  negated it to recover one must stop.
- **Renames** (#15): `b`/`c` → `stem_b`/`stem_c` (there are two Weibull curves and the
  unmarked pair was the source of a real error), `g1_TF24` → `cost_scale_TF24`.
- **`umol_per_mol_to_Pa` is no longer a namespace-scope constant** (#15). It was
  `0.1013` = 101.3 kPa in disguise; it is now the member `Leaf::umol_per_mol_to_Pa_`,
  derived per call from `atm_kpa_`, so the model is self-consistent away from sea level.

## Behaviour changes

Each landed with its own measured blast radius against the golden file; see the PRs.

- **Four exits no longer leave stale state** (#15, #26). `Leaf` is reused for every
  individual in a patch, so any output a branch declined to write became the previous
  plant's value. Fixed: the three `set_shutdown_state` call sites, the `assim_max_ < 0`
  exit, `soil_consumption_` cleared with `.assign` not `.resize`, and
  `dprofit_droot_collar_psi` (which segfaulted without a prior solve and returned a
  spurious gradient in the reversed-gradient state).
- **The collar bracket is clamped to `root_psi_crit`** (#24, plant #584). The clamp
  compared a magnitude against a signed potential and could never bind. The window is
  empty at this package's defaults and 1.2 MPa wide at plant's.

## Downstream

`plant` requires `leaf (>= 0.1.0)`; anything below it will not compile against plant's
`feature/consume-leaf-package`.
