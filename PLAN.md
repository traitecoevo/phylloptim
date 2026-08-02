# Next steps

## Status, 2026-08-03

Tracked as issues in [traitecoevo/leaf_cpp](https://github.com/traitecoevo/leaf_cpp/issues);
this file keeps the reasoning behind each one.

**Done**

| item | what |
|---|---|
| **1** | Validated against plant. The swap is bit-identical: plant's full suite 2364 pass / 0 fail / 0 error on both builds, `test-leaf.r` 218 expectations, SCM regression bit-identical across 78/78 nodes. Harnesses live in `tests/validate/`. |
| **8** | λ and `g1_eff` reported as outputs, λ verified against finite-difference `dA/dE`, and the multi-layer λ identity implemented and verified. |
| **9** | Closed-form fast path (`leaf/closed_form.hpp`): 6.3× with one Newton step, 27× for the explicit β₂=1/c form. Default off; accuracy characterised. |
| **10a** | Renames done: `stem_b`/`stem_c`, `cost_scale_TF24`, `root_carbon_per_layer`. `R` and `n` gone from the public namespace. |
| **10b** | Eight dead entities removed; `set_physiology` 14 → 10 arguments; the leaf is now purely intensive; 13 temperature-response parameters made settable. |
| **10c** | The hidden hard-coded atmospheric pressure fixed (`umol_per_mol_to_Pa` derived from `atm_kpa_`). |
| **15** | CI: gcc/clang × Linux/macOS, no R needed. Plus a golden-file regression baseline over 288 operating points, compared bit-exactly. |
| **5** | **Decided: leave XAD as it is.** It arrives via odelia, both packages want the same version, and only forward mode is used so nothing needs linking. No action. |

Changes that alter results or the API live on **`feature/api-cleanup`**, not on
`master`, so `master` stays a drop-in for plant. That branch is 7 commits and is
where 10a, 10b, 10c and the shutdown fix sit.

**Not ours**

| item | what |
|---|---|
| **2** | The shutdown-state leak is **being fixed in plant** (#578, with #577 as the `resize`/`assign` half). We port plant's fix rather than maintaining our own. ⚠️ `feature/api-cleanup` already carries an independent fix — it must be reconciled with plant's, preferring plant's, or the two will diverge. |

**Remaining** — filed as issues; the item column points back into this document.

| issue | item | what | note |
|---|---|---|---|
| [#2](https://github.com/traitecoevo/leaf_cpp/issues/2) | 7b | Extract the soil/root supply path behind an interface | **do first** — #3, the multi-layer λ, and `kmax(h)` all sit on top. Design question settled: **no template needed**, use `std::variant`; measured, see 7b |
| [#3](https://github.com/traitecoevo/leaf_cpp/issues/3) | 7a | Make λ(state) pluggable | needs #2 |
| [#4](https://github.com/traitecoevo/leaf_cpp/issues/4) | 11 | Template `Leaf` on its scalar type | deletes the hand-maintained AD replicas |
| [#5](https://github.com/traitecoevo/leaf_cpp/issues/5) | 6 | R interface (RcppR6) | absorbs the `Control` struct and the dropped-field cleanup |
| [#6](https://github.com/traitecoevo/leaf_cpp/issues/6) | 12 | Calibration vignette, then inversion | needs #4 for trait gradients |
| [#7](https://github.com/traitecoevo/leaf_cpp/issues/7) | 13 | Energy balance, full cut | leaf-to-air VPD is the cheap win; free convection is not worth it |
| [#8](https://github.com/traitecoevo/leaf_cpp/issues/8) | 10a | Signed-vs-magnitude potentials in the type | design question, not a rename |
| [#9](https://github.com/traitecoevo/leaf_cpp/issues/9) | 3 | Finish the plant-side integration | validated already; decisions + one unchecked consumer |
| [#10](https://github.com/traitecoevo/leaf_cpp/issues/10) | 2 | Port plant's shutdown fix | ⚠️ reconcile with our branch fix |
| [#11](https://github.com/traitecoevo/leaf_cpp/issues/11) | 4 | Drop the last R coupling | belongs upstream in odelia |
| [#12](https://github.com/traitecoevo/leaf_cpp/issues/12) | 15 | `R CMD check` integration, C++ API docs | small |
| [#13](https://github.com/traitecoevo/leaf_cpp/issues/13) | 1 | The unexplained 1 ULP | low value; does not affect plant |
| [#1](https://github.com/traitecoevo/leaf_cpp/issues/1) | 14 | Decide the package name | publication framing is item 14 below |

---

## 1. Validate against plant — DONE. The swap is bit-identical.

**The decisive result.** Swapping plant's own leaf for this package changes
**nothing**:

- plant builds clean on `feature/consume-leaf-package` against the installed
  header-only package — zero warnings, no source changes beyond the shim.
- plant's own `tests/testthat/test-leaf.r` passes **unchanged, 218 expectations**.
  That exercises the RcppR6 bindings and the `plant::Leaf` alias, not just the maths.
- A full **SCM regression is bit-identical across 78 of 78 recorded numeric
  nodes**, for both a one-species and a two-species run, including the entire
  collected trajectory — height, mortality, fecundity, storage,
  `net_mass_production_dt`, `competition_effect`, `opt_psi_stem`, `root_mass` — and
  the ODE step sequence. Reproduce with `tests/validate/scm_regression.R` twice and
  `tests/validate/scm_compare.R`.

That last one is the check that mattered, because the SCM integrates the leaf solve
through an adaptive stepper and a discrete node-splitting schedule: a perturbation
far below tolerance could still flip a refinement decision, and that would be
visible. It did not.

**A separate, unresolved 1 ULP.** `tests/validate/compare_with_plant.R` compares a
standalone C++ binary against plant's R-bound leaf over 288 operating points and
finds 585 of 2592 values differing at 1-2 ULP (worst 2.2e-16). The source is
arithmetically identical — a normalised body diff over all 44 shared functions
found only three `size_t` -> `int` loop counters, which cannot change arithmetic,
plus `transpiration_full_integration` (adaptive Simpson by design, off this path).

**The cause is not established, and I am not going to pretend otherwise.** Ruled
out by experiment: the plant version/branch; compiler flags (`-std=c++20` vs
`gnu++20`, `-g`, `-O0`..`-O3`, `-DNDEBUG`, `-fPIC`, `-ffp-contract` on/off);
inlining (`-fno-inline`); and the odelia header version (local checkout vs
installed — these do differ by 97 lines, but only in an adaptive `construct()`
the leaf never calls). An earlier version of this item blamed translation-unit
structure; **the SCM result contradicts that**, since that comparison does change TU
structure and came out bit-identical. So: 1 ULP, cause unknown, confined to the
standalone harness, and below every solver tolerance. Worth a note, not a blocker.

For scale on why 1 ULP is safe here: forcing `-ffp-contract=off` on one side alone
moves the disagreement to **3e-4**, because these nested solvers amplify
perturbations up to about `GSS_tol_abs` (1e-3). Four orders of magnitude separate
"rounding" from "bug", which is what makes the result interpretable.

**Two methodological traps, both hit and both worth recording.** First, the initial
comparison used whatever plant was installed — a Jul 24 build of a *different
branch*, carrying the ATLS thermal-damage layer, while this package came from the
Jul 31 PM branch. It gave identical numbers, because ATLS is default-off and
genuinely bit-identical when off, but that was luck. Build the reference from the
extraction commit; the recipe is in the script header. Second, I read a
"`-ffp-contract=off` on both sides" run as evidence about the source when
`R_MAKEVARS_USER` had silently not applied the flag to plant's build. Check the flag
reached the compile line.

### Still outstanding under this item

- The SCM regression above runs `max_patch_lifetime = 5`, plant's own test scenario.
  A longer run is **not** needed to guard against drift: the two builds are
  bit-identical at every recorded value including the ODE step sequence, and the
  computation is deterministic, so there is no perturbation to accumulate. What a
  longer, drier run would add is **state-space coverage** — hydraulic shutdown, the
  arid corner, tall trees — i.e. branches the 5-year scenario may never enter. That
  matters little on `main`, where those branches are unchanged, and a lot on
  `feature/api-cleanup`, which deliberately changes the shutdown path. So the long
  run belongs to item 2's "quantify before merging", not to this item.
- ~~Run plant's full test suite on the branch~~ **DONE, and it is clean:
  2364 pass, 0 fail, 0 error, 7 skip — IDENTICAL on both builds.**

  Use the right harness. `test_check("plant")` (what `R CMD check` runs via
  `tests/testthat.R`) parents the test environment on the package *namespace*, so
  RcppR6 internals are visible. `test_dir()` parents it on `globalenv()` and
  produces **93 spurious "could not find function" errors** across 25 files, because
  `Parameters`, `Species`, `Node`, `SCM`, `QAG` and friends are unexported. An
  earlier version of this item recorded those 93 as "pre-existing"; they were an
  artefact of my invocation, not a property of plant. Filed as plant #586.

  Also: plant sets `Config/testthat/parallel: true`, and the workers do not inherit
  a non-default `R_LIBS`, so pass `TESTTHAT_PARALLEL=false` when testing an
  installed build in a temporary library. It otherwise fails with an opaque
  `cli_abort` backtrace that says nothing about libraries.

  And verify the library every  run: `/tmp` was cleaned mid-session and one run silently
  fell back to the site-library plant (a different branch). `find.package("plant")`
  at the top of every script is cheap insurance.

- Optionally, track down the 1 ULP. Low value: it does not affect plant, and the
  remaining suspects are narrow.

### Original text, kept for the record

**Blocks everything else.** The conversion from plant was mechanical
(`inline`-ify, hoist includes, swap `util::stop` for a throw, swap `NA_REAL` for
a quiet NaN, split the constants out, retype three loop counters) and none of it
should move a floating-point result. But *should* is not *does*, and R was not
available in the session that did the conversion, so nothing here has been
compared against plant's compiled build.

The regression values currently pinned in `tests/cpp/test_leaf.cpp` were produced
**by this implementation**, not by plant. They guard against future drift; they do
not establish correctness.

Do this:

- Install plant and this package, and run plant's `tests/testthat/test-leaf.r`
  against a plant build that consumes this package (branch
  `feature/consume-leaf-package`, see item 3). It should pass unchanged.
- Compare a fixed grid of operating points — vary `psi_soil` × `PPFD` ×
  `atm_vpd` × number of soil layers — between plant's `Leaf` and this one, and
  require **exact** equality on `opt_psi_stem_`, `root_collar_psi_`,
  `assim_colimited_`, `transpiration_`, `stom_cond_CO2_`, `profit_` and
  `soil_consumption_`. Anything that differs is a conversion bug, not a
  tolerance question.
- Then run a full SCM regression against a stored plant baseline.

Two known and deliberate exceptions to bit-identity:

- `transpiration_full_integration()` now uses adaptive Simpson
  (`leaf/quadrature.hpp`) rather than plant's compiled QAG. It converges on the
  same integral but is not bit-identical. It is a spline-fidelity diagnostic and
  is not on the production path.
- Three loops over `max_soil_layer` were retyped from `size_t` to `int` to clear
  `-Wsign-compare`. Identical for `max_soil_layer >= 0`; the previous form would
  have looped ~1.8e19 times had it ever run before `set_physiology` set the
  value, so the new form is also strictly safer.

Once this is done, replace the pinned values in `test_leaf.cpp` with values
generated from plant and say so in the comment.

## 2. Shutdown-state leak — PLANT'S TO FIX; we port it

**Being addressed in plant** (issues #578 and #577). We do not maintain a separate
fix; when plant's lands we port it.

⚠️ **Reconcile before merging `feature/api-cleanup`.** That branch already carries
an independent fix, written before it was clear plant was taking this on. Two
different fixes to the same defect is exactly how the stem-vs-root vulnerability
mix-up happened elsewhere. Prefer plant's version and drop ours, or diff them
deliberately — do not merge both. Our fix and its measured blast radius (240
mismatches = 48 shutdown rows × 5 flux fields, `psi_stem`/`collar`/`profit`
untouched) are in that branch's commit, which is useful as a cross-check on
whatever plant does.

The diagnosis below is kept because it is the evidence, and because the
"`-R_d_`, not zero" point is easy to get wrong.

`set_shutdown_state()` writes only `root_collar_psi_`, `opt_psi_stem_` and
`profit_`. It does not reset `transpiration_`, `assim_colimited_`,
`stom_cond_CO2_`, `ci_`, `E_up_` or `soil_consumption_`; and `set_physiology()`
does not either, because `setup_clean_leaf()` is called from the constructors
only. So when a leaf hits the hydraulic-shutdown path, those outputs still hold
whatever the *previous* solve on the same object left there.

That matters because plant holds one persistent `Leaf` per `TF24_Strategy`
(`tf24_strategy.h`) and drives every node, every height and every timestep
through it. Demonstrated:

```
psi_soil=  4.0  E= 3.1769e-06  A= 1.7432  soil_cons[0]= 0.000176347  profit= -3.7056
psi_soil= 12.0  E= 3.1769e-06  A= 1.7432  soil_cons[0]= 0.000176347  profit= -8.3846   <- stale
psi_soil= 20.0  E= 3.1769e-06  A= 1.7432  soil_cons[0]= 0.000176347  profit= -8.3846   <- stale
```

`psi_crit` is 5.87 MPa, so the last two rows are shut down and should be
reporting zero water use. Instead they report the 4 MPa solve's fluxes, and
`soil_consumption_` is what feeds plant's patch water balance. On a fresh object
the same hole shows from the other side: the fluxes are never written at all and
stay at the NA sentinel.

Pinned deliberately by `test_shutdown_leaves_stale_state_known_defect()` so that
fixing it is a visible change rather than a silent one.

The fix is to have `set_shutdown_state()` zero the flux outputs (and set `ci_` to
the compensation point, as the other zero-transpiration branches in
`set_leaf_states_rates_from_psi_stem` already do). Two things to settle first:

- **Does it change plant's published results?** It will, wherever shutdown is
  reached — which the dry-margin sweeps do reach. Quantify before merging, and
  land it in plant as its own reviewed change rather than smuggling it in with
  the packaging.
- **Should it be zero or NA?** Zero is right for `soil_consumption_` and
  `transpiration_` (the leaf genuinely uses no water). For `assim_colimited_`,
  note that `profit_` is already set to `-R_d_ - hydraulic_cost_TF(psi_crit)`, so
  the consistent assimilation is `-R_d_`, not zero.

File this against plant as well as here.

## 3. Have plant consume this package

Branch `feature/consume-leaf-package` in plant does this: delete
`inst/include/plant/leaf_model.h` and `src/leaf_model.cpp`, add
`LinkingTo: leaf`, and provide a compatibility shim so that `plant::Leaf` and the
handful of leaf constants plant reads (`kg_per_mol_h2o`) keep resolving. Nothing
in plant's own sources should need to change beyond that shim.

Remaining work on that branch: build it, run the test suite, and decide whether
to keep the `plant::Leaf` alias permanently or migrate plant's ~12,000 lines of
generated RcppR6/RcppExports code to `leaf::Leaf`. The alias is much cheaper and
costs nothing at runtime.

**There is at least one external consumer of plant's leaf headers to check
against.** `Falster-stomatal_analytical_analysis/notes/tf24_closed_form_bench.cpp`
is an Rcpp translation unit that `#include <plant/leaf_model.h>`, does
`using plant::Leaf;`, links against the installed `plant.so`, and reads public
`Leaf` members (`l.b`, `l.c`, `l.beta2`, `l.transpiration(...)`). The shim keeps all
of that working unchanged — which is the main argument for having written a shim
rather than migrating call sites. Checked: it declares locals named `n` and `R`,
which merely shadowed `plant::n` and `plant::R`, so dropping those two names is
safe and in fact removes the shadowing. Build that file against the branch as part
of verifying it.

## 4. Drop the last R coupling

The leaf model is R-free. odelia is not: `odelia/ode_util.hpp` includes
`RcppCommon.h` for `Rcpp::stop`/`Rcpp::warning` and for `Rcpp::as`/`wrap`
declarations on `odelia::util::index`, and `odelia/interpolator.hpp` includes it
for `util::stop`. That is the only reason this package cannot claim to be plain
C++ outright.

`tests/cpp/shim/RcppCommon.h` is a 15-line stand-in that satisfies all of it,
which both keeps the test suite R-free and specifies exactly how small the
problem is. The upstream fix is for odelia to guard its Rcpp bits — e.g. an
`ODELIA_NO_R` switch, or moving `stop`/`warning` to a throw the way
`leaf/util.hpp` does. File against odelia.

## 5. Where XAD comes from — DECIDED: leave it

**Decision: no action.** XAD arrives via odelia, which vendors it; both packages
live in this family and want the same version, and only *forward* mode is used
(`xad::fwd<double>`), which needs no tape and therefore no linking — confirmed by
`nm`: plant's `leaf_model.o` has zero `xad::Tape` symbols. The alternatives below
are recorded in case that ever stops being true.


XAD arrives via odelia, which vendors it. That is fine while both packages live
in this family and want the same version, and it is why this package does not
vendor a second copy. But it means `leaf` depends on an ODE solver in order to
get an autodiff library, which is odd for outside users and will read as odd in a
paper.

Options, roughly in order of preference:

- Leave it. Cheapest, zero divergence risk, and plant needs odelia anyway.
- Split odelia's vendored XAD into its own tiny `xad` R package that both depend
  on. Clean, and helps anyone else in the R world who wants XAD.
- Vendor a second copy here. Standalone, but two copies in one family is a
  version-skew trap.

Note that only *forward* mode is used (`xad::fwd<double>`), which needs no tape
and therefore no linking — confirmed by `nm`: plant's `leaf_model.o` has zero
`xad::Tape` symbols. Any packaging choice can rely on that.

Same question for the spline: `odelia::interpolator` is header-only and already
templated on its scalar type, so it is a good dependency. Only revisit if item 4
stalls.

## 6. Give it an R interface

Right now this is a `LinkingTo`-only package, like `BH`. To be useful to leaf
physiologists — the audience that would otherwise reach for `plantecophys` — it
needs `leaf::Leaf` callable from R.

plant already has the bindings: `inst/RcppR6_classes.yml` describes the `Leaf`
class and RcppR6 generates the glue. Move that definition here and generate
against this package. Then plant's copy can go too.

While doing it, take the chance to give the R side a saner surface than the
C++ one: `set_physiology()` takes fourteen positional arguments, which is
tolerable from a strategy that calls it once and painful from a console.

## 7. Make the components swappable

**This is the argument for the package being a package rather than a file**, and
it is one refactor approached from two directions: make the *cost function*
swappable (7a) and make the *water supply path* swappable (7b). Both serve the
same end — running alternative formulations against identical drivers — and both
are what a fair model comparison requires.

### 7a. The stomatal / optimality formulations

Three formulations already live in `leaf_model.hpp`, but only one of them is a
real citizen:

| formulation | what exists | status |
|---|---|---|
| **TF24** gain-risk | `hydraulic_cost_TF`, `profit_psi_stem_TF`, `optimise_psi_stem_TF`, and the multi-layer `find_root_collar_psi` | the production path |
| **Sperry et al. (2017)** | `hydraulic_cost_Sperry`, `profit_psi_stem_Sperry`, `optimise_psi_stem_Sperry` | second class — `optimise_psi_stem_Sperry` is hardwired to `psi_soil_[0]`, so single-layer only, and nothing routes to it |
| **Medlyn et al. (2011)** | `medlyn_model_gs`, `solve_medlyn_ci_numerical`, `solve_medlyn_ci_analytical` | second class — its own comment says it is "NOT used by the TF24 compute path"; it bypasses the hydraulic solve entirely |

So you cannot presently run the same drivers through Sperry and TF24 and compare,
which is the obvious thing to want.

**But do not build the dispatch around the cost function.** The
`Falster-stomatal_analytical_analysis` project (`atelier/2-research/active/`, draft
manuscript *"The marginal cost of water as a common currency for stomatal
optimality models"*) has already worked out the right abstraction, and it is one
level deeper. Its central result: six widely used models — Cowan-Farquhar,
Medlyn USO, Prentice least-cost, Sperry gain-risk, TF24, Potkay — all maximise a
profit, therefore all satisfy the *same* first-order condition

```
dA/dE = λ
```

and differ **only** in the function λ(state). Given λ, each collapses to the
Medlyn USO functional form, `ci/ca = g1_eff/(g1_eff + sqrt(D))` with
`g1_eff = sqrt(3·Γ*·Patm/(1.6·λ))` — so USO "is not a model but the generic
solution of the family". Verified against plant: solving `dA/dE = λ_TF24`
reproduces TF24's optimum to **0.02% max relative error in ci/ca across 767
operating points**.

That changes the design. What should be pluggable is **λ(state)**, not the whole
optimiser:

```cpp
struct MarginalCost {                       // concept
  double lambda(const Leaf& l, double psi) const;
  double dlambda_dpsi(const Leaf& l, double psi) const;   // analytic
};
```

Six models then become six small functions rather than six solvers, they share
one tested numerical core, and — the point — a comparison is *guaranteed* to be
apples-to-apples, because only λ differs. That is a far stronger claim than "we
implemented them all in one package", and it is the claim the manuscript needs
code to support. The λ table is `main.tex` Table `tab:lambda`; hand-written
`lambda_TF24` and `dlambda_TF24` already exist in that project's
`notes/tf24_common.R` and `notes/tf24_closed_form_bench.cpp`.

Remaining design notes:

- Compile-time dispatch is *plausible* here in a way it was not for 7b, but it is
  still unmeasured — **measure before committing to it.** The relevant difference
  is that the cost core really does inline: `profit_psi_stem_TF`,
  `hydraulic_cost_TF` and `assim_colimited` have no out-of-line symbol in
  `test_golden` at all, so a virtual λ would be introducing an indirect call where
  today there is no call. That is the opposite of the supply path, where the
  measurement (7b) found the call already out-of-line and the dispatch free. Use
  `tests/cpp/bench_solve.cpp`; the recipe is in 7b. Composes naturally with item 11.
- **Prentice least-cost is the one that does not fit cleanly**, and it is not
  Medlyn. Its cost is a *ratio*, `(aE + bV)/A`, not an additive `λE`, so it is only
  a member of the family after a transformation; the project records that
  `A − λE` was "verified numerically" to reproduce its optimal ci, but that
  verification code is not in the repo. Reproduce it before relying on it.
- Medlyn is easy under this scheme, not awkward as previously noted here: it *is*
  the generic solution, with λ constant and fitted. The existing standalone
  `solve_medlyn_ci_*` becomes the λ = const case.
- Keep `optimise_psi_stem_*`'s single-layer forms as the unit-test entry points
  even after the multi-layer versions exist; they are much easier to reason about.
- Beware the **two vulnerability curves**. Stem (`b`, `c`) drives
  `hydraulic_cost_TF`; root (`root_b`, `root_c`) drives uptake. Using the root
  parameters for the cost gave that project a published-draft λ ∝ ψ^3.02 where it
  should have been ψ^0.64. See item 10.

Why this matters beyond tidiness: every R package in this space commits to one
*hydraulically explicit* scheme, or to none (see
[COMPARISON.md](COMPARISON.md)), so none of them can compare the formulations
where the live argument actually is. A package that runs four against identical
drivers, at 4 µs a solve, is a different and more interesting contribution than a
fourth implementation of one of them.

### 7b. The soil-and-root water supply path

**Does the leaf need to know about soil layers?** No — and this is worth acting
on, because the separation is already almost clean.

Measured on the current code. The soil/root transport is 257 of the 1,117 lines of
member-function body (**23%**), plus 13 state members and a block of
`set_physiology`, dominated by `E_from_Soil_to_Root_Collar` (153 lines) and
`dE_from_soil_dpsi_collar` (58). But the *gas-exchange core is entirely
soil-agnostic*. Every one of these never touches soil state:

`transpiration`, `transpiration_to_psi_stem`, `proportion_of_conductivity`,
`psi_stem_to_ci`, `assim_colimited`, `assim_rubisco_limited`,
`assim_electron_limited`, `electron_transport`, `hydraulic_cost_TF`,
`hydraulic_cost_Sperry`, `profit_psi_stem_TF`, `profit_psi_stem_Sperry`,
`set_leaf_states_rates_from_psi_stem`, every temperature and energy-balance
function, and every Medlyn function.

The coupling runs through a single scalar. `find_psi_stem_from_psi_root` is six
lines:

```cpp
E_from_Soil_to_Root_Collar(psi_root, psi_soil);          // soil -> E_up_
double psi_stem = transpiration_to_psi_stem(E_up_, psi_root);
```

So the soil enters the optimisation *only* as a supply function
`E_up = f(P_collar)` — plus its derivative, which
`dE_from_soil_dpsi_collar` already provides. And note that
`optimise_psi_stem_TF` and `optimise_psi_stem_Sperry` already run off
`psi_soil_[0]`, i.e. a single scalar potential: the soil-free path exists.

**So: isolate it, but keep it here — do not push it up into plant.** The reason
not to move it to plant is that the coupling point is not at plant's level. The
golden-section search evaluates the supply function at every candidate collar
potential, inside the inner loop. plant would have to inject a callback into the
leaf optimiser's hot path, which is a worse boundary than the one that exists now.

Extract it instead as a swappable component *within* the package:

```cpp
struct SupplyPath {                              // concept, not a base class
  double uptake(double P_collar) const;          // E_up, kg H2O m-2 leaf s-1
  double duptake_dpsi(double P_collar) const;    // analytic derivative
  double wettest_potential() const;              // for bracketing the solve
};
```

with (at least) two implementations: `MultiLayerRoots` — today's behaviour, with
per-layer vulnerability, root resistance and gravitational head — and
`SinglePotential`, which is one ψ_soil and either infinite or constant
conductance.

#### Does it have to be a template? Measured: no.

**This item previously asserted that `Leaf` had to be templated on the supply
path so the calls would still inline, and that a virtual base or a
`std::function` would put an indirect call inside a loop running ~10³ inner
evaluations. That assertion was wrong, and it was wrong in the same way the
`area_leaf` coupling claim was: by reasoning instead of measuring.**

Two findings, the first of which settles it on its own.

**1. The supply path is not inlined today, so there is no inlining to preserve.**
`E_from_Soil_to_Root_Collar` has 18 out-of-line call sites in `test_golden` at
`-O2` and clang inlines it into none of them — it is ~150 lines and carries
string-building error paths. Nor is it a threshold accident: at
`-mllvm -inline-threshold=2000` (roughly 10× the default) it is *still*
out-of-line at 22 call sites. The spline `eval` it calls per layer is not
inlined either (clang: `cost=405, threshold=225`). The hot path already pays a
call per layer per evaluation.

**2. So the dispatch is free.** `tests/cpp/bench_solve.cpp` times the 288-point
golden grid; the body of `E_from_Soil_to_Root_Collar` was left byte-for-byte
untouched and only the way it is *reached* was varied, with the bench checksum
identical across all four builds to confirm no arithmetic moved:

| dispatch | µs/solve | vs direct |
|---|---|---|
| direct call (today) | 3.52 | — |
| `std::function` | 3.54 | +0.6% |
| virtual base | 3.56 | +1.1% |
| `std::variant` + `std::visit` | 3.61 | +2.6% |

Reproducible to ±0.01 µs at `reps=2000` (±0.5 µs at `reps=40`, which is why an
early 40-rep run looked like noise in both directions). Aggressive inlining is
if anything *worse* — `-O2 -mllvm -inline-threshold=2000` gives 3.66 µs and
`-O3` with the same 3.62 µs, i.e. code bloat costs more than the calls save.

**So use a plain composed class. Do not template `Leaf`.** Templating would turn
`Leaf` into `Leaf<Supply>`, which breaks plant's `plant::Leaf` alias and its
~12,000 lines of generated RcppR6 code (item 3), for a measured benefit of zero.

**Constraint on which runtime mechanism: `Leaf` must stay copyable.** plant's
`make_strategy_ptr(TF24_Strategy s)` takes the strategy **by value**, and
`TF24_Strategy` holds `Leaf leaf;` as a member (`tf24_strategy.h:375`), so a
`std::unique_ptr<SupplyPath>` member would break plant's build. That leaves
`std::variant` (+2.6%, value semantics, no heap — verified `Leaf` stays
copy-constructible and copy-assignable) or a virtual base with a cloning copy
constructor (+1.1%, but boilerplate and a heap allocation per `Leaf` copy).
**Prefer `std::variant`** while there are only two implementations; 2.6% is a
fair price for value semantics, and it is dwarfed by the 6.3×/27× the
closed-form path already offers (item 9).

⚠️ **This result does not transfer to item 7a.** It is specific to the supply
path. The *cost* core is the opposite case: `profit_psi_stem_TF`,
`hydraulic_cost_TF`, `assim_colimited` and `transpiration_to_psi_stem` have **no
out-of-line symbol at all** — they are fully inlined into the golden-section
loop. Making λ virtual could therefore cost real time where making the supply
path virtual does not. Measure 7a separately with the same harness; do not cite
this table for it.

Four reasons this is worth doing rather than merely tidy:

1. **It is what makes model comparison fair.** To compare cost functions you must
   hold the supply side fixed; to compare supply representations you must hold the
   cost function fixed. Right now neither is possible because they are one object.
   And the alternatives 7a wants to add — Medlyn, Prentice least-cost,
   Cowan-Farquhar — are all formulated against a *single* soil water potential, so
   without this the comparison is either unfair or impossible.
2. **The bare-leaf user should not have to build a root system.** Someone coming
   from `plantecophys` has one ψ_soil and no root-mass profile.
   `set_physiology`'s fourteen arguments include three parallel soil vectors;
   `SinglePotential` removes that barrier to entry entirely.
3. **It quarantines the bug-prone part.** Both currently open defects in this code
   — plant #577 (`resize` that should be `assign`) and #578 (the shutdown leak,
   item 2 here) — are in the soil path. It is the most intricate 23% of the model
   and it currently has no independent test surface.
4. **It keeps a genuinely novel capability.** Multi-layer root water uptake with
   per-layer vulnerability curves and gravitational head has no counterpart in
   `plantecophys`, `tealeaves` or `bigleaf` — `plantecophys` explicitly documents
   that soil-to-root conductance is not implemented. This is an asset, not
   baggage; the point is to make it *optional*, not to lose it.

**Independent confirmation that this is the right interface.** The
`Falster-stomatal_analytical_analysis` project derives the multi-layer correction
to the marginal cost of water and verifies it as an identity
(`notes/tf24_multilayer_lambda.R`, central-difference ratio 1.0000):

```
λ_multi = λ_TF24 · [1 + kmax·f(ψ_r)/S]        S = dE_up/dψ_r
```

`S` — the root-network conductance at the collar — is *exactly* `duptake_dpsi` in
the interface above. So the two methods this interface needs are the two the
multi-layer theory needs, which is a good sign the boundary is in the right place.
Note the size of the effect: single-layer λ **always understates** the marginal
cost, by a factor of 2 to 12, and the correction flattens the height scaling of
`g1_eff` from `h^-0.30` to `h^-0.15`. That project records the multi-layer λ as
"the one genuine gap, and the blocker for any plant implementation" — so
implementing it here, behind this interface, unblocks their work.

Do 7b before 7a: the supply interface is the thing the λ functions plug into.

## 8. Report λ and g1_eff as first-class outputs

**Small, and the analytical project asks for it first.** Its
`notes/proposed changes to plant.md` item 1 is "Report `g1_eff` as an aux
variable… Do this one first."

The leaf currently reports `ci_`, `gs`, `A`, `E`, `psi`, `profit_`. It does not
report the marginal cost of water λ, or the equivalent Medlyn slope
`g1_eff = sqrt(3·Γ*·Patm/(1.6·λ))`. Those two are the quantities that make the
model *comparable to the literature*: `g1_eff` is directly comparable to fitted g1
values in the Lin et al. (2015) database and to `plantecophys::fitBB` output, and λ
is the common currency of item 7a.

Both are cheap — λ is analytic given the cost function, and `g1_eff` is one
`sqrt`. Add them as members alongside the existing outputs, and expose them
through whatever R interface item 6 builds.

One unit trap to document while doing it: plant uses **1.67** for the H₂O:CO₂
stomatal diffusion ratio (`H2O_CO2_stom_diff_ratio` in `leaf/constants.hpp`),
whereas Medlyn (2011) and the g1 literature use **1.6**. That is a 2.2% offset in
g1, which matters when comparing against fitted values. Either expose the ratio as
a settable parameter or report `g1_eff` both ways; do not leave it implicit.

## 9. Add the closed-form fast path

**Speed, and it is already written.** The analytical project derived a closed-form
approximation to the TF24 optimum and benchmarked it in C++ against the real
`plant::Leaf` (`notes/tf24_closed_form_bench.cpp`, 226 lines):

| solver | per solve | speedup |
|---|---|---|
| exact `optimise_psi_stem_TF` | 2.611 µs | 1× |
| power law + 1 Newton step | 0.241 µs | **10.8×** |
| explicit form, β₂ = 1/c | 0.056 µs | **47×** |

At β₂ = 1/c the solve is fully explicit and 0.051 of the 0.056 µs is
`set_physiology`, i.e. the leaf solve itself has essentially vanished. That is a
much larger win than anything else on this list — for comparison, header-only
conversion and LTO both measured at zero (see the family memory note), and the
multi-layer case is a *three*-level nest, so the project expects "the prize there
is larger than 10.8×, not smaller".

That project also notes the C++ "reads only public `Leaf` members, so dropping it
in as a `Leaf` method is close to copy-paste". Bring it over as a selectable
solver: exact search / power law + k Newton steps / explicit β₂ = 1/c.

Four things not to get wrong:

- **The Newton step count is 1, deliberately.** k = 2 is *worse* in the tail.
  Do not "improve" it.
- **The guard tests an output** (`ci/ca > 0.5`), so it must be applied post hoc
  with a fallback to the exact solve — you cannot branch on it up front.
- **Report the realised speedup, not the best case.** With a fallback fraction φ
  the realised gain is `1/[φ + (1−φ)/10.8]` — 3.6× at φ = 0.2, not 10.8×.
  Measuring φ on a real water-limited scenario is an open item over there.
- **Smoothness of the argmax is a hard constraint.** plant chose golden-section
  over Brent specifically because its argmax varies smoothly with inputs, which
  the demographic growth-rate gradient depends on (the comment survives in
  `leaf/optimize.hpp`). Any replacement solver must preserve that; the project
  measures it as roughness in dA/dh — 0.0015 closed form against 0.0011 exact.

The larger prize is that the closed form is **analytically differentiable**, so
the demographic gradient could become exact rather than finite-differenced. That
project observes this would make TF24f's whole acclimation-tracking apparatus
redundant. Which makes this item and item 12 the same argument arriving from two
directions.

## 10. Fix the API while it is still cheap

Two halves: get the *names* right (10a) and get the *input set* right (10b). Both
are far cheaper now, with no downstream users, than after item 6 builds an R
interface and item 3 lands in plant.

The governing principle for 10b: **the leaf takes leaf-level quantities only.**
Anything that is whole-plant allometry gets computed on the plant side and passed
in already reduced. `leaf_specific_conductance_max` is the model to copy — plant
computes `kmax = K_s·θ/(h·η_c)` in `tf24_strategy.cpp:377` and hands the leaf a
scalar, so which conductance-versus-height model is in force is plant's business,
not the leaf's. That matters concretely: the open question over whether α = 1 in
`kmax ~ h^-α` is defensible at all (item 14) can be settled without touching this
package.

### 10a. The naming hazards

A new package with no downstream users is the one moment when renames are free.
These are not cosmetic — each has already cost someone real numbers, and the list
is from `notes/proposed changes to plant.md` §7 plus what surfaced during the
extraction.

- **`b` / `c` → `stem_b` / `stem_c`.** This is the expensive one. There are *two*
  Weibull curves: stem (`b`, `c`) driving `hydraulic_cost_TF`, and root
  (`root_b`, `root_c`) driving uptake. The unmarked default is the stem, which is
  not obvious, and the analytical project used root parameters for the stem cost
  and carried λ ∝ ψ^3.02 into a manuscript draft where it should have been ψ^0.64.
  Never leave an unmarked default for a parameter that exists in two versions.
- **`mass_root_prop` → `root_carbon_per_layer`.** It is not a proportion of mass.
- **`g1_TF24` → `cost_scale_TF24`.** It is not a g1 and invites confusion with
  Medlyn's g1 — doubly so once item 8 starts reporting an actual `g1_eff`.
- **Signed-versus-magnitude water potentials belong in the type, not a comment.**
  The convention is currently held together by a comment block above
  `E_from_Soil_to_Root_Collar` and by suffix conventions
  (`psi_soil_` positive, `psi_soil_inverted_` negative). A one-line strong type —
  or at minimum a consistent naming rule enforced in review — removes a whole
  class of sign error. This code already has form here: plant #584 is a dead
  `std::max` clamp caused by a sign slip.
- **`R` and `n` are already gone** from the public namespace (see item 1). Worth
  recording *why* it mattered: the analytical project's
  `tf24_closed_form_bench.cpp` declares locals `const double n = l.c*l.beta2 - 1.0`
  and a Newton residual `R`, both inside a scope where `plant::R` and `plant::n`
  were visible. It compiles only because the locals shadow them. Removing the
  namespace-scope names makes that file strictly safer.

Do this before item 6 builds an R interface, so the R names are right the first
time, and coordinate with plant, since renames cross the shim.

### 10b. Shrink the input set

Audited `set_physiology`'s fourteen arguments, the constructor's nineteen, and the
constants header. Findings, in order of how free they are.

**Already dead — two removed, three blocked on plant.** These are stored and never
read:

| | status |
|---|---|
| `root_mass_` (member) | never even *assigned*. **Removed.** |
| `vcmax_25_to_jmax_25` (constant, 1.67) | zero uses, like the `n` constant. **Removed.** |
| `rho` (`set_physiology` arg) | stored, never read |
| `a_bio` (arg) | stored, never read |
| `sapwood_volume_per_leaf_area` (arg) | stored, never read — and plant *computes* it (`pars.theta·height·η_c`, `tf24_strategy.cpp:383`) purely to throw away |

The last three are still `access: field` in plant's `inst/RcppR6_classes.yml`, and
plant's `test-leaf.r` asserts they start as `NA`, so removing them needs an RcppR6
regeneration and a test edit. Do it with item 6. Removing them cannot change a
result — nothing reads them — and takes `set_physiology` from fourteen arguments
to eleven.

**Move to the plant side.**

- **`area_leaf`.** Its only use is `1.0/area_leaf_`, normalising soil uptake to a
  per-leaf-area basis (`E_from_Soil_to_Root_Collar`,
  `dE_from_soil_dpsi_collar`). If plant passes root carbon already per unit leaf
  area, the leaf becomes *purely intensive* — no extensive quantity anywhere in it,
  which is a much cleaner contract and one sentence to document.
- **The whole root architecture**: root carbon per layer, `soil_depth`,
  `beta_R_H`, `beta_R_V`, and the hard-coded 1/3 : 2/3 split of root carbon into
  vertical and horizontal components. This is a root-system model living inside a
  class called `Leaf`. It is item 7b, and plant already owns the neighbouring
  parameters (`root_mass_carbon_scale`, `rooting_depth_max`, both file-static and
  TODO-flagged in `tf24_strategy.cpp`).

**Constants that are really parameters.** `constexpr` asserts that a value is not
a modelling choice. Several of these are:

- **The six Arrhenius shape parameters** — `vcmax_ha`, `vcmax_H_d`, `vcmax_d_S`,
  `jmax_ha`, `jmax_H_d`, `jmax_d_S`. These are exactly `plantecophys`'s `EaV`,
  `EdVC`, `delsC`, `EaJ`, `EdVJ`, `delsJ`, which are *user parameters* there — and
  whose **defaults it changed at v1.4 after a literature review**. They are also
  precisely what thermal acclimation modifies, so TF24t needs them mutable. This is
  the strongest case in the list.
- **`H2O_CO2_stom_diff_ratio = 1.67`** where the g1 literature uses 1.6 — see
  item 8.
- **The PM energy-balance block**, most of which the header comments already
  concede: `longwave_net_offset = -40` is an explicit placeholder (plant #581);
  `sw_abs_per_par = 2.0` folds in leaf absorptance, which both `tealeaves`
  (`abs_s`, `abs_l`) and `plantecophys` (`LeafAbs`) expose; `latent_heat_vap` is
  fixed at 25 °C where `plantecophys` makes it temperature-dependent;
  `vol_heat_cap_air` is fixed where `plantecophys` derives air density from `Patm`
  and `Tair`; `aerodynamic_resistance_coef`/`_fixed` are boundary-layer
  coefficients that both comparators derive from leaf size and wind.
- The nine Bernacchi kinetic constants (`gamma_*`, `kc_*`, `ko_*`) are a weaker
  case — enzyme kinetics vary less among species — but `bigleaf` exposes them and
  they do get revised.

**Numerical control is currently mixed in with traits.** `GSS_tol_abs`,
`ci_abs_tol`, `ci_niter` and `vulnerability_curve_ncontrol` are *constructor
arguments sitting among the physiological traits*, and `leaf_temp_min`/`_max` and
`integration_tol_` are scattered elsewhere. Collect them into a `Control` struct —
which is plant's own pattern (`plant/control.h`). Tolerances are not traits, and
mixing them means a trait-calibration loop (item 12) has to know which of its
nineteen constructor arguments are not traits.

### 10c. Two latent inconsistencies found during the audit

**`umol_per_mol_to_Pa = 0.1013` silently hard-codes atmospheric pressure.** It
converts µmol mol⁻¹ to a partial pressure in Pa, which requires the total
pressure: `0.1013 = 1e-6 × 101300 Pa`. So the constant *is* P = 101.3 kPa in
disguise — while `atm_kpa_` is a live, settable input used on the conductance side
(`stom_cond_CO2`, `assim_minus_stom_cond_CO2`, the Medlyn coupling). The model is
therefore **internally inconsistent away from sea level**: raise `atm_kpa` to
simulate altitude and the stomatal side responds while Γ*, Kc, Ko, Km and the ci
root-find bounds all keep assuming 101.3 kPa.

This is the same trap `plantecophys` warns about in bold — that setting `Patm`
alone "does not correct for atmospheric pressure effects on photosynthesis rates"
— except there it is documented and here it is not. Fix by deriving the conversion
from `atm_kpa_` rather than hard-coding it. That *will* change results wherever
`atm_kpa != 101.3`, so quantify first; if the driver never varies pressure the
change is a no-op and can land cheaply.

**`kg_to_mol_h2o` and `kg_per_mol_h2o` are not reciprocals.** `1/55.4939 =
0.0180200` against `kg_per_mol_h2o = 0.018015`, a 0.028% discrepancy. The header
comment already admits this is deliberate, "kept at the historical 0.018015 to
preserve results". Worth unifying once there is a bit-identity baseline (item 1) so
the change can be shown to be 0.028% and nothing else.

## 11. Template `Leaf` on its scalar type

This is the strongest technical argument for the package being header-only, and
the reason to have done the split at all.

`leaf_model.hpp` currently carries hand-maintained templated *replicas* of
`assim_colimited` and `hydraulic_cost_TF` in `namespace detail`, existing only so
that forward-mode XAD can differentiate them. Their comment says they "mirror
`Leaf::assim_colimited` and `Leaf::hydraulic_cost_TF` exactly" — which is a
correctness trap every time the real functions change, and exactly the kind of
duplication that silently rots.

A `Leaf<T>` templated on its scalar type deletes the whole class of problem: the
real functions get differentiated directly, and `dprofit_droot_collar_psi` stops
needing to combine hand-written AD with the implicit function theorem and an
analytic spline derivative. odelia's `basic_interpolator<S>` is already templated
on its scalar type, so the spline side already supports it.

Do this *after* item 1, so there is a bit-identity baseline to check against.
Expect the `double` instantiation to be unchanged and verify it.

There is a second, larger payoff, which is item 12: templating on the scalar type
is what turns "we have AD" into "we can calibrate". See below.

## 12. Demonstrate calibration — and then consider inversion

**The claim.** ~4 µs per solve *and* exact derivatives rather than finite
differences makes this an unusually good target for calibration. Gradient-based
optimisers and Hamiltonian samplers want both: many evaluations, and gradients
that are not numerical noise. And this model is close to the worst case for
finite differencing — a golden-section search wrapped around two nested
root-finds, so the objective is only piecewise smooth in its inputs and a
difference quotient picks up the solver tolerances as much as the physics. plant
already hit exactly this: the analytic `dprofit_droot_collar_psi` exists precisely
because the finite-difference gradient was too noisy for TF24f's acclimation
tracking, and #576 turned out to be an AD-versus-finite-difference branch
asymmetry.

**The honest state of it.** The claim is not yet demonstrated, and it is not yet
fully true:

- AD currently differentiates with respect to the **collar potential only**.
  Calibration needs derivatives with respect to **traits** — `vcmax_25`, `b`, `c`,
  `g1_TF24`, `beta2`, the root parameters. Those do not exist. Item 8 is the
  enabler: a `Leaf<T>` seeded on a trait rather than on `psi_stem` gives them
  directly.
- There is no worked example. Until there is, this is a plausible claim about an
  architecture rather than a capability anyone can check.

**So the deliverable is a vignette**, and it should be a real one:

- Take measured or synthetic `A`, `gs` and `psi_leaf` under a few soil-moisture
  and light conditions.
- Fit hydraulic and photosynthetic traits by gradient-based optimisation, using AD
  derivatives.
- Show the thing that makes the argument: **compare against finite differences**
  — same optimiser, same starting point — and report iterations to convergence,
  wall time, and whether the FD version converges at all at realistic solver
  tolerances. If AD does not win convincingly, say so; that is a useful result
  too, and better found in a vignette than in a paper.
- Report the wall time for the whole fit. If it is seconds, that is the headline.

Do this after item 11. Attempting it before will produce a vignette that
finite-differences trait gradients, which is precisely the thing being argued
against.

**Then the bigger question: inversion.** The clearest capability gap against
`plantecophys` is that `fitaci` inverts A-Ci curves for Vcmax/Jmax/Rd and this
package cannot invert anything. Calibrating *hydraulic* traits — P50-type
vulnerability parameters, root conductances — from gas-exchange and water-potential
data would be a genuinely new capability rather than a reimplementation of
`fitaci`, and the AD story is what would make it tractable. Scope it properly
before committing; it is a project, not a task, and it needs someone to decide
what data it is meant to consume.

## 13. Energy balance — the full cut

The Penman-Monteith path (`use_energy_balance_`, default off) is a deliberate
minimal core. Compared with `plantecophys::PhotosynEB` and `tealeaves` it is
missing a lot, but the missing pieces are not equally valuable and one of them is
expensive. In priority order:

1. **Leaf-to-air VPD.** `saturation_vapour_pressure()` and its slope exist and
   are tested but are *not wired in*: stomatal conductance still divides by the
   prescribed air VPD however far Tleaf has risen above Tair, which understates
   the driving gradient exactly when the energy balance matters. Cheap, and the
   pieces are already there. **Do this one.**
2. **Longwave.** Replace the fixed `-40 W m-2` offset with `eps*sigma*Tleaf^4`
   plus an atmospheric term. Makes `Rn` depend on Tleaf, so it needs either the
   radiation-conductance linearisation that plantecophys uses or an iteration.
3. **Free convection.** A Grashof-number term, as in plantecophys and
   tealeaves. **This is the expensive one, and probably not worth it here.** It
   makes the boundary-layer conductance depend on `|Tleaf - Tair|`, which makes
   the balance implicit and forces an inner iteration *inside* a golden-section
   search that already runs ~10³ inner evaluations per solve. At 4 µs per solve
   in a model that calls it millions of times, that is a poor trade. If leaf
   temperature accuracy under low wind ever becomes the priority, prefer a
   one-step correction over a converged inner solve.

## 14. Naming, home, publication

- **Package name.** `leaf` is clear inside this family and too generic outside
  it. Decide before anything is published: `leafhydro`, `hydroleaf` and
  `leafoptim` all say more.
- **Repository home.** DESCRIPTION points at `traitecoevo/leaf`; create it, or
  change the URL.
- **Paper.** Two candidate framings, and they are not the same paper.

  The **software** paper writes itself from [COMPARISON.md](COMPARISON.md): every
  leaf gas-exchange package in R assumes a stomatal conductance model, and this
  one derives stomatal behaviour from hydraulics instead.

  But the more interesting one is already in progress and is not primarily about
  this package: `Falster-stomatal_analytical_analysis` (`atelier/2-research/active/`),
  draft manuscript *"The marginal cost of water as a common currency for stomatal
  optimality models: size dependence, testable contrasts, and a diagnosis"*,
  targeting *New Phytologist* or *PC&E*. Its argument is that six models share
  `dA/dE = λ` and differ only in λ(state), so USO is the generic solution of the
  family rather than a model — which relocates the empirical question onto the
  *shape* of λ. That paper currently has no code artefact beyond TF24-against-plant
  scripts. Items 7a, 8 and 9 would give it one, and a package that runs six λ
  functions through one tested numerical core is much stronger evidence for the
  unification claim than a symbolic derivation alone.

  Sequencing matters here: **that manuscript is the software's first customer, not
  a downstream user.** Its blockers are ours. In particular it records the
  multi-layer λ as "the one genuine gap, and the blocker for any plant
  implementation" — which is item 7b — and its highest-priority open science
  question is whether TF24's `kmax = Ks·θ/(h·η_c)` (Shinozaki's uniform-diameter
  pipe model, with α = 1 in `kmax ~ h^-α`) is defensible at all, given Koçillari
  et al. 2021 find no height trend in leaf-area-specific conductance across 103
  plants. If α is near zero rather than one, "much of the seedling and tall-tree
  diagnosis may be an artefact of the conductance model rather than the cost
  function". That is a modelling decision this package should make configurable
  rather than hard-code — worth folding into item 7b, since `kmax(h)` sits on the
  supply side.

## 15. Housekeeping

- **CI.** GitHub Actions matrix building `tests/cpp` on gcc and clang, Linux and
  macOS. The suite needs no R, so this is fast and catches the portability
  problems (`long double`, FMA contraction) that have bitten plant before —
  see the arm64 note in `hydraulic_cost_Sperry`.
- **Run the C++ suite under `R CMD check`** so `LinkingTo` consumers get told
  when a header breaks.
- **Doxygen or roxygen for the C++ API.** The header comments are unusually
  good — they explain why, not what — and are worth rendering.
