# Seed: give `leaf` an R interface, as a standalone package (#5)

Handoff brief for starting the R-interface work in a fresh session. Everything
below was established while rebasing and landing #15. Read this first so you don't
re-derive it.

`PLAN.md` is the reasoning behind every open issue and is the authority; this file
is only the *entry point* for this particular next step, plus the things that are
true right now and not yet written down anywhere else.

## Where the repo is, as of 2026-08-03

`master` is `26ab841` and carries **everything**. Five PRs landed in one day:

| PR | what |
|---|---|
| #17, #18 | the supply path became swappable (`MultiLayerRoots`, `SinglePotential`, enum dispatch) — issue #2 closed |
| #19, #22 | the include graph is R-free: `RcppCommon` shim deleted, odelia pinned `>= 0.2.0` — issue #11 closed |
| #20 | `R CMD check` integration + rendered C++ API docs — issue #12 closed |
| #21 | the 1 ULP resolved — it was R's decimal parser, not the model — issue #13 closed |
| #15 | API cleanup: renames, purely-intensive input set, the pressure fix, the shutdown fix |

**The two-branch rule is retired.** There is no longer a `feature/api-cleanup`
where results changes live. Results changes now land on `master` through a PR that
states its measured blast radius against the golden file.

Open issues: **#1** (package name), **#3** (pluggable λ), **#4** (template on
scalar type), **#5** (this one), **#6** (calibration), **#7** (energy balance),
**#8** (signed-vs-magnitude potentials in the type), **#9** (plant-side
integration).

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

## Two live bugs, neither of them fixed

- **plant #577 (`resize` should be `assign`) is live in this package.** `set_physiology`
  calls `soil_consumption_.resize(supply_n_layers(), 0.0)`, but `resize`'s fill
  applies only to newly added elements and the uptake loop writes only up to
  `max_soil_layer` — the deepest *rooted* layer, not the deepest layer. So a solve
  with fewer rooted layers than the last leaves the tail stale. Reproduced:

  ```
  tree     (roots in 3 layers): 0.000522694 0.00026097 0.000172362
  seedling (roots in 1 layer) : 0.000633884 0.00026097 0.000172362
                                            ^^^^^^^^^^^^^^^^^^^^^ carried over
  ```

  **#15's shutdown fix does not cover it** — that runs only on the shutdown path,
  this leaks on every path. One word, but it changes results, so it wants its own
  commit and its own measured blast radius. Andrew measured 33.78% of cohort-time
  records affected on a production run.
- **plant #578 is still open upstream**, even though the fix now lives here. The
  reconciliation went *downstream*: plant had no fix written. Landing both in plant
  is issue #9.

## Where the plant-side work is written up

PLAN item 3 now carries a full survey of what #15 costs plant: 8 hand edits plus a
regeneration, exactly one hand-written C++ break (`tf24_strategy.cpp:445`), and
three corrections worth reading before starting — `rho_`/`a_bio_` are also bound and
removed, `c_r_V_`/`c_r_H_` need no edit, and the `area_leaf` division cancels so it
must not be computed. The real risk there is not the compile errors but the two
places plant compiles clean and the numbers move.

## Recommended order from here

1. **#9's plant-side edits** — plant currently tracks this package's `master` via
   `Remotes:`, so plant is broken against `master` until they land. This is the only
   genuinely urgent item.
2. **#5**, this item, with the header-only decision made first.
3. **#4** (template on scalar type), which #6 needs.
4. **#3** (pluggable λ), now cheaper because #5 removed the plant coupling.
   ⚠️ #2's "dispatch is free" result does **not** transfer: the cost core is fully
   inlined where the supply path was already out-of-line, so it needs its own
   measurement with `tests/cpp/bench_solve.cpp`.
