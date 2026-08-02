# Developer guide — `leaf_cpp`

A header-only C++ leaf gas-exchange model in which stomatal behaviour emerges from
hydraulics. Extracted from the TF24 strategy in
[traitecoevo/plant](https://github.com/traitecoevo/plant).

Read alongside:

- **[README.md](../README.md)** — what it is and how to use it
- **[PLAN.md](../PLAN.md)** — status table, then the reasoning behind every open issue
- **[COMPARISON.md](../COMPARISON.md)** — how this differs from `plantecophys`, `tealeaves`, `bigleaf`
- **[issues](https://github.com/traitecoevo/leaf_cpp/issues)** — the work queue; PLAN.md is the *why* behind each

Family context lives in [`plant-meta`](https://github.com/traitecoevo/plant-meta).

## Layout

```
inst/include/leaf.hpp          umbrella header — one include is the whole library
inst/include/leaf/
  leaf_model.hpp               the Leaf class: ~2000 lines, everything inline
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
```

`bench` reports min-of-N over the golden grid. Use `reps=2000` (the default)
before believing a small difference: reproducibility is ±0.01 µs there but ±0.5 µs
at `reps=40`, which is wide enough to invent or hide a few-percent effect.

No R needed. Dependency headers are found via `Rscript` if R is installed, else a
sibling `odelia/` checkout and Homebrew Boost. Override with
`make BH_INC=... ODELIA_INC=...`.

CI runs gcc and clang on Linux and clang on macOS.

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

## Branches

- **`master`** — stays a drop-in replacement for plant's own leaf. Validated: plant's
  full suite is 2364 pass / 0 fail / 0 error, and an SCM regression is bit-identical
  across 78/78 nodes.
- **`feature/api-cleanup`** — everything that changes results or the API: the renames,
  the shrunken input set, the pressure fix, and an independent shutdown fix that
  **must be reconciled with plant's** (issue #10).

Keep them separate. `master` being a drop-in is what lets plant adopt this without
a coupled review.

## Hazards, each of which has cost someone real numbers

1. **There are TWO Weibull vulnerability curves.** Stem (`stem_b`, `stem_c`) drives
   `hydraulic_cost_TF`; root (`root_b`, `root_c`) drives uptake. They used to be the
   unmarked `b`/`c` plus `root_b`/`root_c`, and the companion analysis used the root
   parameters for the stem cost and carried λ ∝ ψ^3.02 into a manuscript draft where
   it should have been ψ^0.64. Never leave an unmarked default for a parameter that
   exists in two versions.
2. **Signed versus magnitude water potentials.** `psi_soil_` holds positive
   magnitudes; `psi_soil_inverted_` holds signed (negative) potentials. The
   convention is held together by comments, not types (issue #8). `dE_from_soil_dpsi_collar`
   differentiates with respect to the *signed* potential, so it returns a negative
   number where a conductance is wanted — that one produced a negative λ before it
   was caught. plant #584 is a dead clamp from the same class of slip.
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
6. **No Rcpp in the leaf.** `util.hpp` throws `std::runtime_error` instead of
   `Rcpp::stop`, and uses a quiet NaN instead of `NA_REAL`. The only R touchpoint left
   in the include graph is odelia's `ode_util.hpp`; `tests/cpp/shim/RcppCommon.h` is a
   15-line stand-in that both keeps the tests R-free and specifies exactly what needs
   fixing upstream (issue #11). Don't reintroduce Rcpp.

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
