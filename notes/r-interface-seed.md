# Plan: give `leaf` an R interface, as a standalone package (#5)

Implementation plan plus handoff brief. The context was established while landing
#15, and updated after #25/#26 and the plant-side catch-up. Read this first so you
don't re-derive it; the staged plan is at the end.

`PLAN.md` is the reasoning behind every open issue and is the authority; this file
is only the *entry point* for this particular next step, plus the things that are
true right now and not yet written down anywhere else.

## Where the repo is, as of 2026-08-04

`master` carries **everything**. Five PRs landed on 2026-08-03:

| PR | what |
|---|---|
| #17, #18 | the supply path became swappable (`MultiLayerRoots`, `SinglePotential`, enum dispatch) — issue #2 closed |
| #19, #22 | the include graph is R-free: `RcppCommon` shim deleted, odelia pinned `>= 0.2.0` — issue #11 closed |
| #20 | `R CMD check` integration + rendered C++ API docs — issue #12 closed |
| #21 | the 1 ULP resolved — it was R's decimal parser, not the model — issue #13 closed |
| #15 | API cleanup: renames, purely-intensive input set, the pressure fix, the shutdown fix |

then two more on 2026-08-04, which are why this section needed rewriting:

| PR | what |
|---|---|
| #26 | ported plant `develop`'s leaf fixes: four stale-state exits, incl. **plant #577** — see the retraction below |
| #25 | **one representation for water potential: positive magnitudes throughout.** `psi_soil_inverted_` deleted, `root_collar_psi_` → `opt_root_psi_` (positive) |

**The two-branch rule is retired.** There is no longer a `feature/api-cleanup`
where results changes live. Results changes now land on `master` through a PR that
states its measured blast radius against the golden file.

Open issues: **#1** (package name), **#3** (pluggable λ), **#4** (template on
scalar type), **#5** (this one), **#6** (calibration), **#7** (energy balance),
**#9** (plant-side integration), **#24** (the dead clamp, plant #584).

**#8 is answered, not open.** It asked for the signed-vs-magnitude convention in
the type system. #23 built exactly that (`Psi` / `AbsPsi`), it worked, it was
bit-identical, and it was **closed unmerged**: it described a two-convention model
instead of removing one, and no type can cover `dE_from_soil_dpsi_collar`, which is
a derivative and therefore neither type. #25 deleted the second convention instead.
Do not reintroduce the types.

## The goal, and why it is next rather than #3

Make `leaf::Leaf` callable from R, so the package stands on its own for the
audience that would otherwise reach for `plantecophys` — rather than existing only
as plant's `LinkingTo` dependency.

Chosen ahead of #3 (pluggable λ) deliberately. Three reasons, in descending order
of how much they matter:

1. **#5 dissolves hazard 7, which makes every later refactor cheaper.** plant's
   `inst/RcppR6_classes.yml` binds most of `Leaf`'s state with `access: field`, so
   RcppR6 emits `obj_->psi_soil_` as getter *and* setter into plant's generated
   glue. That is why #17 needed coupled commits in two repos, and why #15 drifted
   into conflict. Move the bindings here and plant stops naming those fields at
   all. Doing #3 *after* #5 avoids paying the coupled-review tax on every step of
   it.
2. **It ships work that already exists.** `SinglePotential` is written, wired and
   tested, and is **unreachable from R** — `supply_kind_` is a C++-only API. That
   is the bare-leaf, one-ψ_soil use case which was reason 2 in PLAN 7b for building
   it at all. PLAN 7b-iii stage 2 explicitly deferred wiring it up "with item 6,
   where the whole R interface gets designed, not before".
3. **#6 is blocked on #4 anyway**, so it was never actually next.

## The prerequisite is discharged

PLAN 10a says: do the renames **before** item 6 builds an R interface, so the R
names are right the first time. That was the reason #5 could not start. **#15
merged, so it is done** — `stem_b`, `stem_c`, `cost_scale_TF24`,
`root_carbon_per_leaf_area`, and a 10-argument `set_physiology` are what the R API
should be built against.

**#25 finished the job, and it matters more for an R API than for the C++ one.**
Every ψ the R interface will expose is now a positive magnitude in MPa — one
convention, asserted at the input boundary rather than documented. Two consequences
for the YAML:

- the collar output is **`opt_root_psi_`, positive**. The old `root_collar_psi_` is
  gone by design: a flipped sign under the same name would let an old analysis read
  the wrong value silently, where a rename gives an error.
- `E_from_Soil_to_Root_Collar`, `find_root_psi`, `find_psi_stem_from_psi_root` and
  `dE_from_soil_dpsi_collar` take the soil state as an *argument*, and #25 changed
  what that argument means without changing any signature. They validate it and
  stop on a negative entry. **If you give the R side a friendlier wrapper, keep that
  check** — it is the only thing standing between an old script and a wrong number.

## One decision to make before writing any YAML

**leaf is currently header-only by design, and #5 changes that.** Three documents
assert the old answer and will need revising together, not retrofitting:

- `DESCRIPTION`: "The package ships headers only: there is no compiled code and
  nothing to link against."
- `NAMESPACE`: "no R-level API and no compiled code, so nothing is exported and no
  shared library is loaded."
- `.github/workflows/cpp-tests.yml`: pure C++, no R.

RcppR6 adds `src/`, a shared library, `R/`, and needs an R job in CI. This is
reconcilable with #11 — the *headers* stay R-free and plant keeps `LinkingTo` them,
with the R layer sitting on top — but settle the framing first. Note #11 only just
made the include graph R-free; don't undo that by letting Rcpp back into
`inst/include/`.

## What to build against

plant already has the bindings, and they are the reference implementation: move
`Leaf` out of plant's `inst/RcppR6_classes.yml` and generate against this package.
plant's copy can then go, which is the hazard-7 payoff above.

Two design notes carried from PLAN item 6 and 7b-iii:

- **Give the R side a saner surface than the C++ one.** Ten positional arguments is
  tolerable from a strategy that calls it once and painful from a console.
- **Wire up `supply_kind_` / `single_`** so `SinglePotential` is reachable. PLAN
  7b-iii flags the footgun: a naive setter silently invalidates a configured root
  network. Design it, don't just expose it.
- **Absorb the `Control` struct** (PLAN 10b): `GSS_tol_abs`, `ci_abs_tol`,
  `ci_niter`, `vulnerability_curve_ncontrol`, `integration_tol_`,
  `leaf_temp_min`/`_max`. Deferred *to* this item because it means surgery on a
  19-argument constructor that #5 redesigns anyway.

## Things that will waste your time if you don't know them

These are in `.claude/CLAUDE.md` in full; the three that bite hardest here:

1. **The golden file is the safety net and is bit-exact only on macOS/arm64.** Only
   run `make -C tests/cpp golden` deliberately — running it after an accidental
   change rubber-stamps the change. Magnitudes: ~1e-16 is reassociation, ~1e-4 is a
   real difference. If you need a magnitude, read the summary line, not the FAIL
   lines (that mistake has been made twice).
2. **`test_check("plant")`, not `test_dir()`**, when testing an installed plant, and
   pass `TESTTHAT_PARALLEL=false` in a non-default library. `test_dir` produces 93
   spurious "could not find function" errors that look exactly like a broken build.
3. **Bench interleaved, never sequentially.** Between-process noise is ~±0.1 µs, and
   a sequential A/B got the sign wrong once already.

## ⚠️ RETRACTION: the "two live bugs" section was true for one day

This file used to say plant #577 was live here and that plant had no fix to port.
**Both halves expired on 2026-08-04.** plant merged #585 into `develop`, which fixed
#577 *and* two further stale-state exits, and #26 ported all of it here. #578's fix
had already landed here in #15.

The generalisable part, because this document and PLAN item 2 both got caught the
same way: **a "checked, there is nothing upstream to port" conclusion has a shelf
life of days on an actively developed sibling.** Re-check before acting on it, not
just before writing it down. The check is
`git log <base>...origin/develop -- <the files we forked>`.

One live bug remains: **#24, the dead clamp** (plant #584). `prepare_collar_solve`
takes `std::max(root_crit, -supply_psi_crit())` where it wants
`std::min(root_crit, root_psi_crit)`, so it can never bind. #25 kept it live and
annotated deliberately, because fixing it moves results and deserves its own
measurement rather than riding along with a representation change.

## The plant-side work is DONE, and its survey was wrong in three ways

PLAN item 3's survey of what #15 costs plant (8 hand edits plus a regeneration, one
hand-written C++ break) was accurate as far as it went, and it went too far in
predicting the *risk*. The work landed on plant's `feature/consume-leaf-package` on
2026-08-04. Three things the survey did not predict, each of which cost real time:

1. **A third compile break, in a dependency.** plant needs **odelia at `master`**
   (`d8235d1`), not the `>= 0.2.0` its DESCRIPTION asks for: plant #585 made
   `Patch::ode_rates` non-const while odelia 0.2.0's `r_ode_rates` takes the system
   by `const&`. plant's `develop` does not compile against released odelia either,
   so this is not caused by the swap — but it presents as eight bewildering errors
   inside odelia's headers. Filed as traitecoevo/odelia#48.
2. **The pressure fix is NOT inert for plant, and it dominates everything else.**
   The survey said 10c "should be inert at 101.3 kPa, and plant's tests all use
   that". The *tests* do; `TF24_Environment` sets the **`atm_kpa` driver to 100.5**,
   which nothing had checked. Measured on the one-species SCM scenario, deriving the
   ppm→Pa conversion from `atm_kpa` moves offspring production **+2.4%**, against
   **+0.10%** for the whole rest of the swap combined. TF24's `scientific_version`
   went 4 → 5 for it. This package's golden grid evaluates at 101.3 and is therefore
   blind to it by construction — a good reminder that the golden file bounds *this*
   package's behaviour, not plant's.
3. **"The `area_leaf` division cancels — do not compute it" is right for plant's C++
   and a trap for its tests.** In plant the cancellation is exact
   (`mass_root()` is `pars.a_r1 * area_leaf`, so the per-leaf-area carbon is just
   `root_mass_carbon_scale * pars.a_r1`). But `test-leaf.r` passed absolute carbon
   and `area_leaf` as separate arguments, and dropping the argument without dividing
   leaves a root system 20× too weak at `area_leaf = 0.05`. It compiles, runs, and
   moved the critical-demand collar potential from −0.685 to −2.57 MPa **while
   leaving the zero-uptake collar untouched** — that one is scale-invariant, so of
   two regression guards sitting side by side only one fired.

The stronger verification, worth reusing: `origin/develop` built and tested in a
worktree against the same installed dependencies passes its own suite with zero
failures, which is what licenses attributing every failure on the branch to the
swap rather than to the environment. And the dry scenario gateway
(`PLANT_RUN_SCENARIOS=1`) passes on the branch against develop's blessed baseline —
so despite the +2.4%, every scenario's success/failure classification is unchanged.

## Recommended order from here

1. ~~**#9's plant-side edits**~~ — **done**, see the section above. plant's
   `feature/consume-leaf-package` compiles and passes against `master`.
2. **#5**, this item, with the header-only decision made first.
3. **#4** (template on scalar type), which #6 needs.
4. **#3** (pluggable λ), now cheaper because #5 removed the plant coupling.
   ⚠️ #2's "dispatch is free" result does **not** transfer: the cost core is fully
   inlined where the supply path was already out-of-line, so it needs its own
   measurement with `tests/cpp/bench_solve.cpp`.

## The staged plan for #5

Five stages, each landing on its own. The ordering is chosen so that the two
decisions that are hard to reverse — whether the package stops being header-only,
and what the R-facing names are — come first and cheapest, and so that no stage
leaves `master` unable to build plant.

**Stage 0 — settle header-only, in a docs-only PR.** Nothing else can start until
`DESCRIPTION`, `NAMESPACE` and `.github/workflows/cpp-tests.yml` agree on an answer
(see "One decision to make" above). Write the answer down before writing YAML: the
headers stay R-free and plant keeps `LinkingTo` them, with `src/` and an R job added
alongside. Cheap to land, and it is the thing a reviewer will otherwise argue about
in the middle of a 2000-line diff.

**Stage 1 — the YAML and the generated glue, no design changes.** Copy `Leaf` out of
plant's `inst/RcppR6_classes.yml` verbatim, minus the four fields #15 deleted, with
`root_collar_psi_` → `opt_root_psi_`. Target: `Leaf()` constructs from R and
`set_physiology` / `find_root_collar_psi` / the field getters work. Deliberately
**not** friendlier than plant's version yet — this stage exists to prove the build,
the CI job and `R CMD check` are sound with one moving part, not two.

⚠️ **The C++ suite is the regression baseline for all of this, and it is blind to
the R layer.** Add an R-side test that reproduces two or three golden operating
points through the R API and compares to the recorded `%.17g` values. Without it, a
YAML mistake that mistranslates an argument produces plausible numbers and no
failure.

**Stage 2 — a friendlier surface, as a thin R layer over stage 1.** This is where
the design work goes, and keeping it *above* the generated glue means it can be
rewritten without touching CI or the build:

- **`set_physiology`'s ten positional arguments become a named, defaulted call.**
  Ten is tolerable from a strategy that calls it once and painful from a console.
- **Keep the `psi_soil >= 0` check.** #25 made every psi a positive magnitude and
  asserted it at the input boundary; a friendly wrapper that quietly accepts a
  signed vector would undo the one thing protecting an old script from a wrong
  number. Surface the error, don't smooth it over.
- **Absorb the `Control` struct** (PLAN 10b): `GSS_tol_abs`, `ci_abs_tol`,
  `ci_niter`, `vulnerability_curve_ncontrol`, `integration_tol_`,
  `leaf_temp_min`/`_max`. Deferred *to* this item because it means surgery on a
  19-argument constructor that #5 redesigns anyway.
- **A one-call convenience entry point** — drivers in, operating point out, as a
  data.frame. This is the thing that makes the package legible to someone who would
  otherwise reach for `plantecophys`, and `tests/cpp/test_golden.cpp`'s `solve()` is
  already the shape of it.

**Stage 2b — take resistances, not root carbon.** Coupled with plant, so it wants its
own PR pair, but it belongs here because it is the same surgery on the same argument
list.

`set_physiology` currently takes `root_carbon_per_leaf_area` and calls
`roots_.set_root_network_from_carbon(...)` internally. **The solve never touches root
carbon.** `uptake_impl` and `duptake_dpsi` read exactly two vectors —
`network_.r_R_H_min` and `network_.r_R_V_sum` — plus `grav_head_z_` and
`max_soil_layer`. Everything else the carbon path produces (`c_r_V`, `c_r_H`, `r_R_V`)
is carried only because plant exposes it through RcppR6.

So taking carbon makes the leaf own four things that are not gas exchange:
`beta_R_H`, `beta_R_V`, `dz_`, and the 1/3 : 2/3 horizontal/vertical split. Two of
them are constructor arguments — 2 of the 19 that stage 2 is redesigning anyway.
Taking resistances instead deletes all four from the leaf's surface, and it is exactly
the move `leaf_specific_conductance_max` already makes: plant computes
`kmax = K_s*theta/(h*eta_c)` and hands over a scalar, because which
conductance-versus-height model is in force is not this package's business. Which
root-architecture model is in force is not either. `MultiLayerRoots::set_root_network`
already exists and takes a `RootNetwork` directly; `roots.hpp`'s own comment argues
for this at length, and `set_physiology`'s says the remaining step — hoisting the call
up to plant — "is an API change and belongs with item 10b". It was deferred, not
rejected.

Two objections, both answerable:

- **"It leaves the tested surface."** `roots.hpp` says `root_network_from_carbon`
  stays in the package "so that it remains covered by the golden file — the moment
  this arithmetic crosses the package boundary it leaves the tested surface." Keep the
  function here, as a documented public helper (it already has a value-returning
  overload "for tests and one-off callers", and `test_root_network_from_carbon`
  already exercises it). What moves out is the *call*, not the model. The golden grid
  then calls helper-then-`set_physiology` and covers both.
- **"It costs an allocation per solve."** Only if the caller rebuilds the vectors.
  Building five fresh ones each call measured **+0.074 µs (0.061 → 0.135), about +2%
  of a whole solve**, which is why the in-place `RootNetwork& out` overload exists.
  plant already holds the carbon buffer as a `TF24_Strategy` member and refills it; it
  would hold a `RootNetwork` the same way. Neutral if done that way — and worth
  re-measuring interleaved, per hazard 5, rather than assumed.

The gain for an R user is the real argument: a bare-leaf caller does not have a root
carbon profile, a layer thickness, or an opinion about `beta_R_V`. They have a
resistance, or they have `SinglePotential` (stage 3). Asking them for carbon forces
them through a plant-shaped model to get at a leaf.

**Stage 3 — expose `SinglePotential`.** It is written, wired and tested, and
unreachable from R because `supply_kind_` is a C++-only API. That is the bare-leaf,
one-ψ_soil use case that was reason 2 in PLAN 7b for building it at all, and it is
what makes the item 7a model comparison possible: Medlyn, Prentice least-cost and
Cowan-Farquhar are all formulated against a *single* soil potential.

⚠️ PLAN 7b-iii flags the footgun: **a naive `supply_kind_` setter silently
invalidates a configured root network.** Design the transition, don't just expose
the enum. The obvious shape is two constructors (or two `set_supply_*` calls) that
each leave the object fully configured, with no state in which the tag and the
network disagree.

**Stage 4 — delete plant's copy, and collect the hazard-7 payoff.** Once the
bindings live here, plant's `inst/RcppR6_classes.yml` stops naming `Leaf`'s fields
and hazard 7 dissolves: no future member move needs a coupled commit in two repos.
This is a coupled PR pair like the ones #25/#26 just did, and it is the last stage
because it is the only one that can break plant.

Sequencing note: stages 0, 1, 2 and 3 leave plant untouched, so they can land while
plant's `feature/consume-leaf-package` is still in review. **Stages 2b and 4 are
coupled** and must wait for it to merge, or they will be rebasing against a moving
target. 2b before 4 if both are done, since 4 deletes the bindings that 2b changes.

## One measurement to carry forward

It will come up again the next time something is "supposed to be bit-identical":
**boost's TOMS748 is not sign-symmetric.** #25 predicted the golden file would be bit-identical except for an
exactly-negated collar column. It came out 276/288 exactly negated, with 12 rows off
by 1–3 ULP and 5–11 rows per other column moving at ≤5e-16. Hunted rather than
tolerated, and localised: the flux rewrite *is* exactly antisymmetric (`root_zero_E`,
which depends only on `E_up`, is exactly negated), while `root_crit` — from the same
solver on a reversed bracket — is not. Reproduced independently on a synthetic family
of smooth monotone targets: 44 of 192 bracket midpoints are not the exact negation.
Reversing a bracket's orientation costs you the last two bits, and no amount of
care in the model code will buy them back.
