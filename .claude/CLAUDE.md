# Developer guide — `leaf_cpp`

A header-only C++ leaf gas-exchange model in which stomatal behaviour emerges from
hydraulics. Extracted from the TF24 strategy in
[traitecoevo/plant](https://github.com/traitecoevo/plant).

Read alongside:

- **[README.md](../README.md)** — what it is and how to use it
- **[PLAN.md](../PLAN.md)** — status table, then the reasoning behind every open issue
- **[notes/r-interface-seed.md](../notes/r-interface-seed.md)** — handoff brief for the
  next step (#5, the R interface): where the repo is, what is decided, what bites
- **[COMPARISON.md](../COMPARISON.md)** — how this differs from `plantecophys`, `tealeaves`, `bigleaf`
- **[issues](https://github.com/traitecoevo/leaf_cpp/issues)** — the work queue; PLAN.md is the *why* behind each

Family context lives in [`plant-meta`](https://github.com/traitecoevo/plant-meta).

## Layout

```
inst/include/leaf.hpp          umbrella header — one include is the whole library
inst/include/leaf/
  leaf_model.hpp               the Leaf class: the gas-exchange core, everything inline
  roots.hpp                    MultiLayerRoots (the soil → root-collar water supply)
                               plus root_network_from_carbon, the root architecture
                               model that feeds it resistances
  single_potential.hpp         SinglePotential: the other supply path — one ψ_soil
                               and a constant series resistance
  vulnerability.hpp            the Weibull cumulative-integral builder, shared by both
  potential.hpp                Potential (signed) / Suction (magnitude) — the water
                               potential sign convention in the type (hazard 2)
  constants.hpp                physical constants as inline constexpr
  closed_form.hpp              fast approximate solver, default off, not wired in
  quadrature.hpp               adaptive Simpson (replaced plant's compiled QAG)
  util.hpp                     R-free stop()/sentinels
  uniroot.hpp, optimize.hpp    1-D root finders and optimisers
tests/cpp/                     plain-C++ suite, no R, no framework
tests/cpp/golden/              bit-exact regression baseline, 288 operating points
tests/cpp/bench_solve.cpp      timing harness for the collar solve (hazard 5)
tests/validate/                R scripts comparing against plant (needs R)
```

## Build and test

```sh
make -C tests/cpp            # builds and runs both suites
make -C tests/cpp golden     # regenerate the golden file -- see the warning below
make -C tests/cpp bench      # time the collar solve (not part of `make all`)

# compare the golden file with a tolerance instead of bit-exactly. Correct on a
# platform other than macOS/arm64, wrong as a way to silence a real diff.
make -C tests/cpp GOLDEN_ARGS=--cross-platform
```

`bench` reports min-of-N over the golden grid. Use `reps=2000` (the default)
before believing a small difference: reproducibility is ±0.01 µs there but ±0.5 µs
at `reps=40`, which is wide enough to invent or hide a few-percent effect.

⚠️ **The ±0.01 µs is *within* one process. Between processes it is ~±0.1 µs, so
A/B by running one build and then the other is not good enough.** Measured while
sizing issue #2 stage 1: one sequential pair said the refactor was 4% *faster*;
interleaving the two binaries three times each said it was 1.7% slower, which is
what held up. Build both, keep both, and alternate:

```sh
cp bench_solve bench_before      # ...switch branches, rebuild...
for i in 1 2 3; do ./bench_solve | tail -1; ./bench_before | tail -1; done
```

No R needed. Dependency headers are found via `Rscript` if R is installed, else a
sibling `odelia/` checkout and Homebrew Boost. Override with
`make BH_INC=... ODELIA_INC=...`.

CI runs gcc and clang on Linux and clang on macOS.

### The two R-side entry points

```sh
R CMD check .                # runs the C++ suite too, via tests/cpp.R
doxygen                      # renders the C++ API to docs/html/index.html
```

`tests/cpp.R` is not a duplicate of the CI workflow. It compiles with R's
*configured* compiler (`R CMD config CXX20`) against the *installed* headers,
which is what a `LinkingTo: leaf` consumer gets — so a consumer running their own
`R CMD check` is told when one of our headers stops compiling. Expect exactly one
NOTE, `'LinkingTo' field is unused: package has no 'src' directory`: it is
inherent to a header-only package and is not worth adding compiled code to
silence.

It compiles the two sources directly rather than calling `make`, and
`tests/cpp/Makefile` is `.Rbuildignore`d, because `R CMD check` warns about the
GNU extensions the harness uses. The sanctioned fix, `SystemRequirements: GNU
make`, was tried and reverted — installing this package needs no make, and
`LinkingTo: leaf` consumers would inherit a declaration that is false for them.

Doxygen, not roxygen — roxygen documents R objects and this package has none.
Because every comment here is a plain `//`, which Doxygen ignores,
`tools/doxygen_filter.awk` converts them at render time; **the sources are never
touched**, and CI asserts that every non-comment line survives the filter
byte-for-byte. If you want a literal Doxygen command, write a `///` comment and
the filter will leave it alone. Publishing is off until someone sets the repo
variable `PUBLISH_DOCS=true`.

**"If you add a header, you need do nothing" was too optimistic — `potential.hpp`
tripped two comment-shape rules, and the docs job runs with `WARN_AS_ERROR`, so
either one fails the build.** Both are worth knowing before you write the header,
not after:

- **Don't open a file with a dashed `// ----` banner.** The filter turns the
  file's first comment block into a `\file` block, and a following line of `---`
  makes Markdown read the `\file` line itself as a setext heading — *"'\file'
  command is not allowed in section title"*. Open with a prose summary, as every
  other header does; banners are fine anywhere later in the file.
- **Two qualified names joined by a slash reads as an explicit link request.**
  `std::min/std::max` fails with *"explicit link request to 'max' could not be
  resolved"*; bare `std::variant` and `std::vector<double>` are fine and already
  ship. Write them separately, and backtick them — a code span suppresses
  autolinking.

⚠️ **Both are CI-only.** Doxygen 1.17 (Homebrew) accepts the second happily;
ubuntu's 1.9 does not — the same version skew the filter already documents for a
leading `::`. Reproduce locally with the CI command,
`{ cat Doxyfile; echo "WARN_AS_ERROR = FAIL_ON_WARNINGS"; } | doxygen -`, but
don't read a local pass as CI passing.

## The golden file is the safety net — treat it that way

`tests/cpp/golden/operating_points.tsv` records 288 operating points and is
compared **bit-exactly** (`%.17g` round-trips a double). It is what makes a large
refactor of this code checkable rather than hopeful, and it earned that role
repeatedly: it proved three "surely dead" `set_physiology` arguments really were
dead, confined the shutdown fix to exactly 48 rows × 5 fields, and showed the
`area_leaf` change was 2 ULP.

**Only run `make golden` deliberately.** Running it after an accidental change
rubber-stamps the change. If a diff is intended, regenerate and say so in the
commit message with the measured size of the change.

Guide to magnitudes, measured: **~1e-16 is reassociation, ~1e-4 is a real
difference.** These nested solvers amplify perturbations up to about
`GSS_tol_abs` (1e-3), so there is a four-order-of-magnitude gap between rounding
and bug. Anything in between deserves investigation.

**It is bit-exact only on the platform that generated it — macOS/arm64.** libm's
`exp`/`pow` are not bit-reproducible between glibc on x86-64 and Apple's libm on
arm64, and FMA contraction differs too, so cross-platform bit-equality was never
achievable. CI compares bit-exactly on macOS and with `--cross-platform` elsewhere.

**The size of the cross-platform disagreement is the interesting part, and it is
not one number.** The nine reported fields split into two classes three orders of
magnitude apart:

| field | gcc | clang | why |
|---|---|---|---|
| `profit` | 1.85e-06 | 5.87e-07 | it is the maximum itself — well-conditioned |
| the other eight | 5.53e-04 | 2.73e-04 | evaluated at the **argmax** — sqrt-amplified |

The maximum is *flat*: curvature measured directly at the two worst points gives
k ≈ 1.0 and 0.9 in `profit ≈ p* − k(psi_stem−x*)²`. For a flat maximum an error
`dp` in the profit **value** displaces its **location** by `sqrt(dp/k)`, which is
why the well-conditioned column sits three orders below the other. Checked
pointwise — at those points the residuals imply k = 1.01 and 1.70 against the 1.0
and 0.9 measured from the curvature. It is *not* a global identity: the two column
maxima fall at different operating points, so `sqrt(worst profit)` is not meant to
reproduce `worst argmax`.

Two things follow that matter beyond this file:

- **`profit` is the only reported field that is well-conditioned across platforms.**
  If you need a portable check of the solve, check `profit`, not `opt_psi_stem_`.
- **Eight of the nine fields are pinned to 17 digits but only *determined* to about
  `GSS_tol_abs` (1e-3).** Bit-exactness on one platform is still a good drift
  detector — any code change moves the arbitrary choice inside that window and it
  shows — but that is reproducibility, not determinacy. Don't read a bit-identical
  `opt_psi_stem_` as meaning the argmax is known to 17 digits.

⚠️ **The numbers above were wrong twice, the same way both times.** First this note
said the worst cross-platform difference was 1.7e-15 (13 ULP); then, after that was
caught, 2.1e-07 and 4.5e-04. Both came from the 20 lines `test_golden` prints
before it truncates — reading a truncated failure list as if it were the
distribution. The truncated sample is biased toward whichever operating points
happen to be listed first, which is not a random draw.

The figures in the table are the full-grid maxima CI now prints on every run. The
per-class reporting exists so the number never has to be inferred from the printed
failures again — **if you need a magnitude, read the summary line, not the FAIL
lines.**

**Never regenerate the golden file to make another platform pass.** It just moves
the failure to the platform the file came from.

## Branches

**Single-branch now.** `feature/api-cleanup` merged as #15 on 2026-08-03, and with
it the two-branch rule this section used to state ("results changes live on the
branch, `master` stays a drop-in"). Both halves of that rule are spent: `master` is
no longer a source-level drop-in for plant's leaf, and there is no longer a branch
to quarantine results changes onto.

- **`master`** — everything. **Behaviourally** still validated against plant's own
  leaf as of the pre-#15 state: full suite 0 fail / 0 error, SCM regression
  bit-identical across 78/78 nodes. (The pass *count* is environment-dependent —
  2364 when first recorded, 2431 on a rerun with `NOT_CRAN` unset, on both arms of a
  control. Compare against a control run, never against a remembered number.)
  ⚠️ **#15 deliberately changed results**, so that validation now needs redoing
  against plant — see hazard 7 and PLAN item 9.

What replaces the old rule: **a change that moves results lands on `master` through
a PR that states its measured blast radius against the golden file.** #15 is the
worked example — 240 cells from the shutdown fix, 355 at ≤1.2e-13 from
reassociation, split by cause because the two classes are four orders apart. A PR
that regenerates the golden file without saying what moved and why is the thing to
refuse.

## Hazards, each of which has cost someone real numbers

1. **There are TWO Weibull vulnerability curves.** Stem (`stem_b`, `stem_c`) drives
   `hydraulic_cost_TF`; root (`root_b`, `root_c`) drives uptake. They used to be the
   unmarked `b`/`c` plus `root_b`/`root_c`, and the companion analysis used the root
   parameters for the stem cost and carried λ ∝ ψ^3.02 into a manuscript draft where
   it should have been ψ^0.64. Never leave an unmarked default for a parameter that
   exists in two versions.
2. **Signed versus magnitude water potentials — now in the type, on the supply
   side only.** `leaf/potential.hpp` defines `Potential` (signed, ≤ 0) and
   `Suction` (positive magnitude), and the soil → root-collar path is threaded
   with them: `psi_soil_inverted_` is a `std::vector<Potential>`, and
   `begin_solve` / `uptake` / `uptake_at` / `duptake_dpsi` take `Potential`.
   Mixing the two does not compile, and the tests assert that with the detection
   idiom rather than a comment.

   **Three boundaries to know, because the type stops at each of them:**

   - **The demand side is plain `double` magnitudes on purpose.**
     `transpiration`, `psi_stem_to_ci`, `hydraulic_cost_TF`,
     `proportion_of_conductivity`, `profit_*`, `opt_psi_stem_` and `psi_crit`
     have one convention and no recorded incident. Don't "finish the job" here
     without a reason; it buys `.value` across the whole photosynthesis core.
   - **The R boundary is `double`.** plant binds `psi_soil_`,
     `root_collar_psi_`, `opt_psi_stem_`, `psi_stem` as `access: field, type:
     double` (hazard 7), so those members cannot be retyped, and teaching Rcpp
     these types would reintroduce hazard 9. Every plant-bound field and method
     keeps its signature behind a thin adapter — including `prepare_collar_solve`
     and `profit_at_collar_psi`, which plant calls from C++ rather than R. So
     nothing here needs a coupled plant change.
   - **`dE_from_soil_dpsi_collar` is NOT covered.** It differentiates with
     respect to the *signed* potential, so it returns a negative number where a
     conductance is wanted — that produced a negative λ before it was caught.
     A derivative is neither a `Potential` nor a `Suction`, so
     `marginal_cost_water_multilayer` still carries the negation and the comment.
     Closing it needs a conductance type on the supply contract.

   ⚠️ **plant #584 is LIVE here and deliberately left that way.** The clamp in
   `prepare_collar_solve` — `std::max(-root_crit, -root_psi_crit)` — compares a
   magnitude against a signed potential, so it can never bind. Reproduced at
   `psi_soil = 5.90` with `psi_crit = 5.91988`: the bracket runs to 5.906974,
   past `root_psi_crit = 5.870283`. Fixing it moves results, so it is annotated
   in place with the intended line quoted, and filed as **issue #24**. It now
   only compiles by reaching past the types into `.value`, which is the point.
   Don't tidy that line without the blast-radius measurement.

   One more thing the conversion surfaced and did not change: `E_column` parks a
   *magnitude* in `root_collar_psi_` mid-solve as scratch, in a member documented
   as signed in every branch. Nothing reads it before the next write.
3. **Argmax smoothness is a hard constraint, not a preference.** plant chose
   golden-section over Brent because its argmax varies *smoothly* with inputs, and
   the demographic growth-rate gradient depends on that. Any change to the solver
   must re-measure it. The comment in `optimize.hpp` explains why.
4. **The leaf is purely intensive.** Every input is per unit leaf area or a
   dimensionless/intensive driver. Nothing scales with plant size — whole-plant
   allometry (`kmax(h)`, root carbon totals) is computed on the plant side and passed
   in already reduced. Keep it that way; it is a one-sentence contract.
5. **Hot-path discipline — but measure it, because the intuition has been wrong.**
   `find_root_collar_psi` runs ~10³ inner evaluations per solve, and plant calls it
   millions of times. λ and `g1_eff` are *accessors*, not stored state, for this
   reason, and don't add a `pow()` to `set_leaf_states_rates_from_psi_stem`.

   What is *not* true is the blanket "prefer compile-time dispatch, a virtual call
   in the solve would cost". Measured (`make -C tests/cpp bench`), the code splits
   into two halves that answer oppositely:

   - **The soil/supply path is already out-of-line.** `E_from_Soil_to_Root_Collar`
     is inlined into none of its 18 call sites at `-O2`, nor at
     `-inline-threshold=2000`; neither is the per-layer spline `eval`. So making it
     virtual costs **+1.1%**, `std::function` +0.6%, `std::variant` +2.6% — noise.
     Issue #2 needs no template. Forcing more inlining is *worse* (3.66 vs 3.52 µs).
   - **The cost core is fully inlined.** `profit_psi_stem_TF`, `hydraulic_cost_TF`,
     `assim_colimited` and `transpiration_to_psi_stem` have no out-of-line symbol at
     all. Issue #3's pluggable λ lands here, so it is a genuinely different question
     and needs its own measurement.

   Before arguing from inlining, check: `nm -C test_golden | grep <fn>`. No symbol
   means inlined; a symbol plus `bl` call sites in `objdump -d` means it is not.
6. **There are TWO supply paths, chosen by `Leaf::supply_kind_`.** `MultiLayerRoots`
   (default) and `SinglePotential`. Both are held as members and selected by an
   enum, which measured **free** where `std::variant` cost +1.0% — a predictable
   branch in front of an already-out-of-line call disappears into it. Don't
   "tidy" this into a variant without re-measuring; PLAN 7b-iii stage 2 has the
   numbers and the reproduction recipe, including the trap that a
   compile-time-known tag folds the branch away and reports a false zero.

   When adding a third path, the contract is five methods — `begin_solve`,
   `uptake`/`uptake_at`, `duptake_dpsi`, the collar-potential bound, and the layer
   count — plus one rule: **`duptake_dpsi` returning NaN means "fall back to
   finite differences", so returning 0 or throwing silently degrades TF24f's
   acclimation gradient.** A path with no branch kinks simply never returns it.
7. **Moving a public member is a plant API break, because RcppR6 binds fields by
   name.** plant's `inst/RcppR6_classes.yml` lists most of `Leaf`'s state as
   `access: field`, and the generator emits `obj_->psi_soil_` — a getter *and* a
   setter — straight into `src/RcppR6.cpp`. Relocating a member into a sub-object
   therefore breaks plant's generated glue, not just its own sources.

   It is cheap to fix, which is the part worth knowing: for `access: field` the
   template pastes `name_cpp` verbatim after `->`, so a dotted path works and the
   R-side name (taken from the YAML key, not `name_cpp`) does not move:

   ```yaml
   psi_soil_:  {type: "std::vector<double>", access: field, name_cpp: "roots_.psi_soil_"}
   r_R_H_min:  {type: "std::vector<double>", access: field, name_cpp: "roots_.network_.r_R_H_min"}
   ```

   The path can be nested arbitrarily, as the second line shows — so a member can
   be moved as deep as the design wants without the R API noticing.

   Check `git grep 'access: field' inst/RcppR6_classes.yml` in plant before
   moving anything public, and land the two changes together — plant tracks this
   package's **master** via `Remotes:`, so a merge here is what breaks it.
8. **`soil_consumption_` is not cleared between solves (plant #577), and it is
   still live here.** `set_physiology` calls
   `soil_consumption_.resize(supply_n_layers(), 0.0)`, but `resize`'s fill applies
   only to newly added elements, and the uptake loop writes only up to
   `max_soil_layer` — the deepest *rooted* layer, not the deepest layer. So a solve
   with fewer rooted layers than the last one leaves the tail holding the previous
   plant's values, and plant bills them to the patch water balance. Reproduced
   here: a 3-layer tree then a 1-layer seedling on the same `Leaf` gives the
   seedling bit-identical layers 1–2. **The shutdown fix does not cover this** —
   that one runs only on the shutdown path, this leaks on every path. One-word fix
   (`assign`), but it changes results; see PLAN item 2.
9. **No Rcpp in the leaf, and none in the include graph either.** `util.hpp` throws
   `std::runtime_error` instead of `Rcpp::stop`, and uses a quiet NaN instead of
   `NA_REAL`. odelia's `ode_util.hpp` was the last R touchpoint; it was fixed upstream
   (traitecoevo/odelia#44), so `tests/cpp/shim/` is gone and the suite compiles against
   the real headers with no R anywhere (issue #11, closed). Don't reintroduce Rcpp —
   and note the tests will no longer tell you if you do it indirectly by leaning on a
   header that pulls R in, because there is no longer a shim to notice.

## Validating against plant

`tests/validate/` holds the harnesses. Read the header of
`compare_with_plant.R` before using them — it records what has been ruled out.

Three things that will waste your time otherwise:

- **Build the reference plant from the commit this package was extracted from**, not
  whatever is installed. The installed plant may be a different branch; during this
  work it was, carrying the ATLS thermal-damage layer, and gave identical numbers
  only because ATLS is default-off.
- **Print `find.package("plant")` in every script and check it.** `/tmp` was cleaned
  mid-session once and a run silently fell back to the site-library plant while
  reporting plausible numbers.
- **Use `test_check("plant")`, not `test_dir()`**, to test an installed plant.
  `test_check` parents the test environment on the package namespace so RcppR6
  internals are visible; `test_dir` doesn't, and produces 93 spurious
  "could not find function" errors that look exactly like a broken build. Also pass
  `TESTTHAT_PARALLEL=false` when the package is in a non-default library, since the
  parallel workers don't inherit `R_LIBS`. Filed as traitecoevo/plant#586.

## Related work

- **plant** — `feature/consume-leaf-package` consumes this package via a compatibility
  shim aliasing `plant::Leaf`. Issue #9.
- **odelia** — supplies the spline interpolator and the vendored XAD. Only *forward*
  mode is used, which needs no tape and so no linking.
- **The companion manuscript** — `Falster-stomatal_analytical_analysis` in atelier,
  *"The marginal cost of water as a common currency for stomatal optimality models"*.
  **It is this package's first customer, not a downstream user**: its blockers are
  ours (the multi-layer λ is issue #2), and its central result — that six optimality
  models differ only in λ(state) — is what issue #3 is built on.

## Writing commits here

The commit log is doing real work in this repo: several commits record *why* a
change is safe, what its measured blast radius was, and which explanations turned
out to be wrong. Two of the more useful ones retract an earlier claim. Keep that up —
state the measurement, not the intention, and if an earlier commit was wrong say so
plainly rather than quietly correcting it.
