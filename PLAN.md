# Next steps

## Status, 2026-08-03

Tracked as issues in [traitecoevo/leaf_cpp](https://github.com/traitecoevo/leaf_cpp/issues);
this file keeps the reasoning behind each one.

**Done**

| item | what |
|---|---|
| **1** | Validated against plant. The swap is bit-identical: plant's full suite 0 fail / 0 error on both builds, `test-leaf.r` 218 expectations, SCM regression bit-identical across 78/78 nodes. Harnesses live in `tests/validate/`. **And the 1 ULP is resolved** (#21): there was no offending function — all 2352 finite values are bit-identical to plant. The 585 decomposed exactly into 345 from R's decimal parser (not correctly rounded, ~18% of inputs off by 1 ULP) and 240 from the shutdown NA sentinel, which #15 has now fixed. The disagreement was in the measuring instrument. |
| **8** | λ and `g1_eff` reported as outputs, λ verified against finite-difference `dA/dE`, and the multi-layer λ identity implemented and verified. |
| **9** | Closed-form fast path (`leaf/closed_form.hpp`): 6.3× with one Newton step, 27× for the explicit β₂=1/c form. Default off; accuracy characterised. |
| **10a** | Renames done: `stem_b`/`stem_c`, `cost_scale_TF24`, `root_carbon_per_leaf_area`. `R` and `n` gone from the public namespace. |
| **10b** | Eight dead entities removed; `set_physiology` 14 → 10 arguments; the leaf is now purely intensive; 13 temperature-response parameters made settable. |
| **10c** | The hidden hard-coded atmospheric pressure fixed (`umol_per_mol_to_Pa` derived from `atm_kpa_`). |
| **15** | CI: gcc/clang × Linux/macOS, no R needed, plus `R CMD check` and rendered C++ API docs (#20). And a golden-file regression baseline over 288 operating points, compared bit-exactly on macOS/arm64 and with `--cross-platform` elsewhere — see the note under item 1 on why bit-exact cannot be a cross-platform gate, and on the flat-optimum amplification that sets the tolerances. |
| **5** | **Decided: leave XAD as it is.** It arrives via odelia, both packages want the same version, and only forward mode is used so nothing needs linking. No action. |
| **7b** | The supply path is swappable. All four stages of 7b-iii merged (#17 → `10115e1`, #18 → `cfd5dcf`): `MultiLayerRoots`, the resistance interface, `SinglePotential`, and an enum-tag dispatch that measured **free**. Issue #2 closed; #3 and the multi-layer λ are unblocked. |
| **2** | The shutdown-state leak, **fixed and merged** (#15) — and plant **#577 too**, together with two more stale-state exits ported from plant `develop` (#585). See item 2; the reconciliation with plant reversed direction twice. |
| **4** | **The include graph is R-free.** The `RcppCommon` shim is deleted (#19) and odelia is pinned at `>= 0.2.0`, its first release with an R-free core (#22). ⚠️ Restated by item 6a once the R layer landed: the guarantee is *directional* — nothing reachable from `<leaf.hpp>` includes Rcpp — not "no Rcpp under `inst/include/`", which the generated `RcppR6_*.hpp` now break. What enforces it is `cpp-tests.yml` building on runners with no R. |
| **10a tail** | **Signed-vs-magnitude potentials — issue #8, closed.** Answered by **#25**: remove one convention rather than type both. #23 typed them (`Psi`/`AbsPsi`), worked, was bit-identical, and was **closed unmerged** — it described the two-convention model instead of removing one, and no type covers `dE_from_soil_dpsi_collar`, which is a derivative and therefore neither. Do not reintroduce the types. |
| **3** | **The plant-side integration — issue #9, closed.** plant's `feature/consume-leaf-package` (plant #591) compiles and passes against this package's `master`. Its survey was accurate about the work and wrong about the risk in three ways, all recorded under item 3. |
| **6d stage 3** | **`SinglePotential` is reachable from R — issue #32.** `leaf_supply_single()` / `leaf_supply_multilayer()`. The 7b-iii footgun is designed out: no settable tag at any level, both entry points reconfigure completely, and the fields are bound read-only. |
| **6** | **The R interface — issue #5.** `leaf::Leaf` is callable from R: generated RcppR6 bindings tied back to the golden file bit-exactly, then a hand-written surface over them — `leaf_solve()` (drivers in, operating point out, vectorised), `leaf_traits()` / `leaf_control()` splitting the constructor's 19 arguments, `leaf_model()` / `set_drivers()` / `operating_point()`, and `vignette("leaf")`. λ and `g1_eff` are exposed for the first time. The model stays R-free and gained a CMake package, so it is still linkable from C++ or Python. **#32** landed on top of it. **#33** and **#34** are split out and coupled with plant #591. Reasoning in item 6; #31 was found on the way. |

**Everything above is on `master`.** `feature/api-cleanup` merged as
[#15](https://github.com/traitecoevo/leaf_cpp/pull/15) (`26ab841`) on 2026-08-03,
which retires the two-branch split this section used to describe. There is no
longer a branch where "changes that alter results" live, and no rule to keep them
apart: results changes now land on `master` through a PR with their blast radius
measured, which is what the golden file is for.

**What #15 moved**, measured cell-by-cell over the 288-point golden grid and split
by cause, because the two classes are four orders of magnitude apart and should
never be quoted as one number:

| cause | cells | rows | worst relative change |
|---|---|---|---|
| the shutdown fix (`nan` → a real value) | 240 | 48 | n/a — was unset |
| ppm-to-Pa + the `area_leaf` reassociation | 355 | 124 | **1.2e-13** (`profit`) |

240 is exactly 48 shutdown rows × 5 flux fields (`ci`, `assim`, `transpiration`,
`gc`, `e_up`); `psi_stem`, `collar` and `profit` are untouched, so the structural
outputs were already right. The second row is rounding on this guide's own scale
(~1e-16 is reassociation, ~1e-4 is a real difference).

Three things from that work are worth keeping, because they are properties of the
code rather than of the merge:

- **`area_leaf` is out of the *supply contract*, not just out of `set_physiology`.**
  `uptake` / `uptake_at` / `duptake_dpsi` no longer take it, and `inv_area_leaf` is
  gone from `uptake_impl`. So hazard 4 ("the leaf is purely intensive") is now true
  of the interface a *third* supply path would implement, not merely of the public
  entry point — it cannot reintroduce an extensive quantity without changing the
  interface. Consequence: `SinglePotential::resistance_` is **per unit leaf area**,
  documented at the member.
- **Why dropping `area_leaf` was safe is a homogeneity property, not a
  coincidence**, and it is worth being able to state, because "it must cancel
  somewhere" was the first reviewer reaction. Every term in `r_R` is *inversely*
  linear in root carbon — including the cumulative `r_R_V_sum`, which is a sum of
  `1/C[j]` terms — so `r_R` is homogeneous of degree −1 in the carbon vector: scale
  carbon by `1/A` and `r_R` scales by exactly `A`, reproducing the explicit `1/A`
  the old code applied and nothing else. That holds for every layer and every
  operating point, which is what makes it a property rather than luck.
  `beta_R_H`/`beta_R_V` keep their numerical values and are now read as resistances
  per unit leaf area. The residual **2 ULP** is pure reassociation:
  `(dpsi*(1/A))/r_R` versus `dpsi/(r_R*A)`.
- **The rebase preserved arithmetic, and that was checkable rather than hoped.**
  The golden file came through the rebase **byte-identical**, while `master` had in
  the meantime moved the entire supply path into `roots.hpp` behind an enum
  dispatch. Cost: nothing — 3.53 µs/solve on both arms, interleaved ×4 at
  reps=2000.

**Remaining** — filed as issues; the item column points back into this document.

| issue | item | what | note |
|---|---|---|---|
| [#3](https://github.com/traitecoevo/leaf_cpp/issues/3) | 7a | Make λ(state) pluggable | **unblocked** — #2 is done. Measure the dispatch separately: the cost core is fully inlined where the supply path was not, so 7b's "dispatch is free" result does **not** transfer |
| [#4](https://github.com/traitecoevo/leaf_cpp/issues/4) | 11, **11a** | Template `Leaf` on its scalar type | **reordered, and 11a is DONE** (#35 + #36, 2026-08-04). Golden section missed `dprofit = 0` by 6.2e-04, which for hydraulic traits made the trait response **smooth, plausible and sign-inverted** (`root_b`: −2.6e-03 against a true +2.6e-04, arbitrated derivative-free). Replaced by a root-find: residual improved on 240/240 rows, **24.5% faster**, hazard 3 ~1000× smoother. **Still open:** the replica drift (5.53e-4) and `Leaf<T>` for the trait partials, both smaller than originally scoped |
| [#33](https://github.com/traitecoevo/leaf_cpp/issues/33) | 6d | Take resistances, not root carbon | #5 stage 2b, **blocked on plant #591**. Before #34 |
| [#34](https://github.com/traitecoevo/leaf_cpp/issues/34) | 6d | Delete plant's `Leaf` bindings | #5 stage 4, **blocked on plant #591**. The hazard-7 payoff, and the only stage that can break plant |
| [#6](https://github.com/traitecoevo/leaf_cpp/issues/6) | 12 | Calibration vignette, then inversion | **UNBLOCKED** — the gate was 11a, not all of #4, and 11a landed (#35 + #36). A plain central difference now gives ~4 correct digits at any relative step from 1e-8 to 1e-2, so a fit can start. AD is still the vignette's headline result (AD *against* FD), so #4's remainder still improves this rather than gating it |
| [#7](https://github.com/traitecoevo/leaf_cpp/issues/7) | 13 | Energy balance, full cut | leaf-to-air VPD is the cheap win; free convection is not worth it |
| [#28](https://github.com/traitecoevo/leaf_cpp/issues/28) | 13 | Temperature-dependent outgoing longwave in the Penman-Monteith Rn | from plant #581 / #567 review |
| [#31](https://github.com/traitecoevo/leaf_cpp/issues/31) | 6 | `profit_psi_stem_TF` returns a plausible number below `psi_upstream` | found writing the stage 2 vignette. Negative conductance, and profit is **discontinuous** at the boundary. Unreachable from plant's solve, which is why it survived — the first thing an R user does is plot the profit function |
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

**The 1 ULP is RESOLVED, and it was never in the model — issue #13.**
`tests/validate/compare_with_plant.R` reported 585 of 2592 values differing at
1-2 ULP with "cause unknown". The 585 decomposes exactly:

| count | cause |
|---|---|
| 345 | R's decimal parser reading the golden file |
| 240 | the shutdown-state NA sentinel (48 rows × 5 flux fields) |
| **585** | |

**R's string-to-double conversion is not correctly rounded.** `as.numeric`,
`scan` and `read.delim` all share it, and it returns a double one ULP off the
correctly rounded value for about 18% of inputs — `"26.550866314209998"` parses
to `0x1.a8d0593240001p+4` where the correct nearest double is `0x1.a8d059324p+4`.
The golden file is written by C++ at full `%.17g` and read into R; plant's values
are computed in-process and never touch a string. So the parser perturbed one
side only, and the script blamed the two implementations for it.

Reading through `tests/validate/tsv_to_hex.c` — parse with the C library's
`strtod`, re-emit as hex, which R reads exactly — the comparison is now
**bit-identical across all 2352 finite values, worst relative difference 0.0**.
`tests/validate/compare_primitives.R` confirms it independently by calling the
underlying functions directly: 329 values, all bit-identical, from `arrh_curve`
up through the spline-backed `transpiration` and the two that iterate.

This also explains why every hypothesis was ruled out with nothing left standing:
plant version, compiler flags, `-ffp-contract`, inlining, odelia header version,
translation-unit structure — none was ever involved. And it is why the SCM
regression was bit-identical: that comparison never round-trips through a text
file parsed by R.

The other 240 are not arithmetic. `set_shutdown_state` assigns
`root_collar_psi_`, `opt_psi_stem_` and `profit_` and leaves
ci/assim/transpiration/gc/e_up untouched, so golden carries the NA sentinel where
plant carries a number. That is the shutdown-state leak (item 2, plant #578), and
48 × 5 is exactly the blast radius recorded for that fix. They are now reported
in their own column instead of being counted as mismatches.

**Expect zero from now on.** A finite difference is real: the nested solvers
amplify perturbations up to about `GSS_tol_abs` (1e-3), so genuine arithmetic
differences arrive around 1e-4, four orders of magnitude above rounding.

**A second, larger instance of the same effect: the golden file is not portable,
and the amplification is what makes it interesting.** CI's first run to reach
Linux (2026-08-03 — the workflow had been watching `branches: [main]` on a repo
whose default branch is `master`, so it had never executed) found **1761 of 2592
values differing under g++ and 1800 under clang++**. `exp`/`pow` are not
bit-reproducible between glibc on x86-64 and Apple's libm on arm64, and FMA
contraction differs, so a bit-exact cross-platform comparison was never
achievable; the file's real job only ever needed bit-exactness on *one* platform.

The magnitudes split cleanly into two classes:

| field | gcc | clang |
|---|---|---|
| `profit` | 1.85e-06 | 5.87e-07 |
| `psi_stem`, `collar`, `ci`, `assim`, `transpiration`, `gc`, `e_up`, `uptake` | 5.53e-04 | 2.73e-04 |

**That split is the flat-optimum amplification, measured.** `find_root_collar_psi`
maximises profit over the collar potential, and the maximum is flat — curvature
measured directly at the two worst operating points gives k ≈ 1.0 and 0.9 in
`profit ≈ p* − k(psi_stem − x*)²`. For a flat maximum an error `dp` in the profit
*value* displaces the *argmax* by `sqrt(dp/k)`, which is why the well-conditioned
row sits three orders below the other. Checked pointwise: at those points the
residuals imply k = 1.01 and 1.70, against 1.0 and 0.9 from the curvature itself.
Note it is not a global identity — the row maxima fall at different operating
points, so `sqrt(worst profit)` is not meant to reproduce `worst argmax`.

So `test_golden` takes `--cross-platform`, with per-class tolerances (profit 1e-5,
argmax-derived 5e-3 — 5.4× and 9.0× headroom on the measured worst). The profit
tolerance is deliberately not loosened further: 1e-4 is the scale at which a real
change shows, so a profit gate approaching it would gate nothing. CI runs
bit-exact on macOS/arm64 and `--cross-platform` elsewhere. Both modes print the
worst value per class. **Regenerating the golden file to make a second platform
pass is the wrong move** — it just relocates the failure.

Two implications worth carrying:

- **`profit` is the only reported field that is well-conditioned across platforms.**
  For a portable check of the solve, compare `profit`, not `opt_psi_stem_`.
- **Eight of the nine fields are pinned to 17 digits but determined only to about
  `GSS_tol_abs` (1e-3).** Bit-exactness on one platform remains a sound drift
  detector, but it is reproducibility of an arbitrary choice inside the solver's
  tolerance window, not determinacy of the argmax.

⚠️ **These numbers were wrong twice, the same way both times.** This paragraph
first put the worst difference at 1.7e-15 (13 ULP) and called the whole thing
reassociation; once that was caught it said 2.1e-07 and 4.5e-04. Both readings
took the 20 lines `test_golden` prints before truncating as though they were the
distribution. That sample is biased toward whichever points are listed first, and
here those were the well-conditioned ones. The figures in the table are now the
full-grid maxima that CI prints on every run — **for a magnitude, read the summary
line, not the FAIL lines.**

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
  mattered little while the shutdown path was unchanged, and a lot once #15 changed
  it deliberately. **#15 has now merged, so this is outstanding against `master`
  and belongs to item 9** (landing the change in plant), not to this item: the
  shutdown fix moves dry-margin water balance, and a 5-year scenario may never
  enter the branch that moved.
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

- ~~Optionally, track down the 1 ULP~~ **DONE (#13), and the answer was that there
  was nothing to track down in the model: it was R's decimal parser reading the
  golden file. See the resolution above.**

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

## 2. Stale state between solves — all four exits now fixed here

**This section has been wrong twice, in opposite directions, and the history is
worth keeping because the mistake is repeatable.**

1. It first said plant owned the fix and we would port plant's when it landed.
2. Reconciled 2026-08-03: **checked, and plant #578/#577 were both still open with
   no fix written** — no commit on any plant branch touched `set_shutdown_state`,
   and `soil_consumption_.resize` was still `resize`. The earlier reading had
   mistaken *"Andrew has filed it"* for *"Andrew is fixing it"*. So our fix stood
   and plant was downstream of it (#15).
3. **Retracted 2026-08-04: statement 2 has expired.** plant merged
   [#585](https://github.com/traitecoevo/plant/pull/585) into `develop`, which
   fixes #577 (`.assign`) *and* two further stale-state exits, and bumps TF24's
   `scientific_version` to 4 for the shutdown half. So plant did acquire fixes to
   port, five days after we concluded it never would. All of them are now in this
   package (branch `fix/port-plant-develop-leaf-fixes`), each as its own commit
   with its own measurement.

**The lesson, since this is the second time the direction flipped: a "checked, and
there is nothing upstream" conclusion has a shelf life of days on an actively
developed sibling.** Re-check before acting on it, not just before writing it
down. `git log <base>...origin/develop -- <the files we forked>` is the check.

The shutdown half of the fix was merged here as #15 (`26ab841`), before plant had
one. Two things made that safe rather than a land-grab:

- **The semantics were already agreed on plant #578, in writing.** The one comment
  on that issue is Daniel's, and it specifies exactly what our fix implements:
  zero for `transpiration_` / `stom_cond_CO2_` / `soil_consumption_`, the
  compensation point `gamma_ * umol_per_mol_to_Pa` for `ci_`, and **`-R_d_`, not
  zero,** for `assim_colimited_`. So there is no divergence to resolve — the fix
  *is* the position recorded on the issue.
- **It lands on precisely the three exits the issue names, and no others.**
  `set_shutdown_state` is called from exactly three places in
  `prepare_collar_solve`, which are #578's three rows. The fourth early exit
  (`assim_max_ < 0`) does **not** call it — it sets its own operating point and
  calls `E_from_Soil_to_Root_Collar` — and #578 says that one is already correct.
  Fixing inside `set_shutdown_state` therefore cannot touch it. That is the check
  worth repeating if the fix is ever moved.

  ⚠️ **The second bullet's premise — that #578 says the `assim_max_ < 0` exit "is
  already correct" — turned out to be wrong**, and plant #585 fixed it. It sets
  `profit_` but leaves `transpiration_`, `stom_cond_CO2_` and `assim_colimited_`
  holding the previous solve's values, which made `profit_ == assim_colimited_ -
  hydraulic_cost_TF()` false in that one branch (measured: 7.79 against a reported
  -1.47). Now ported. The *scoping* argument still holds — fixing inside
  `set_shutdown_state` genuinely cannot reach that exit — but "cannot reach it" was
  read as "does not need to".

### plant #577 — the `resize`/`assign` half. **Now fixed here** (ported from plant #585)

Was `soil_consumption_.resize(supply_n_layers(), 0.0)` where it wants `.assign`.
`supply_n_layers()` returns `soil_number_of_depths_` (every layer), while the
uptake loop runs to `max_soil_layer` (the deepest *rooted* layer), and `resize`'s
fill argument applies only to *newly added* elements — so layers a plant has no
roots in are never written and never cleared, and plant bills all of them to the
patch water balance. Andrew's protocol, run against this package:

```
tree     (roots in 3 layers): 0.000522694 0.00026097 0.000172362
seedling (roots in 1 layer) : 0.000633884 0.00026097 0.000172362
                                          ^^^^^^^^^^ ^^^^^^^^^^ bit-identical
```

Andrew measured **33.78%** of cohort-time records on a production run carrying at
least one stale layer, and 25.73% carrying three of five.

**Our shutdown fix did not cover this**, and it is worth being clear why, because
the two look like one bug: `std::fill` in `set_shutdown_state` runs only on the
shutdown path, whereas #577 leaks on *every* path, shutdown or not.

**How it was measured, and the trap in measuring it.** The golden file is
**bit-identical** after the fix, all 288 points, and that says nothing either way:
`test_golden` constructs a fresh `Leaf` per grid point on purpose, so
`soil_consumption_` is always empty when `set_physiology` runs and `resize` fills
all of it. The bug needs a *reused* `Leaf`. The measurement is a new test that
solves a three-layer-rooted tree and then a top-layer-only seedling on one `Leaf`
with the same three soil layers — only the rooted depth shrinks, which is why
resizing never noticed. Reverting the one word makes it report 1.31e-04 and
4.75e-05 kg H2O m^-2 s^-1 of phantom uptake in layers 2 and 3. That is the ~1e-4
"real difference" scale, four orders above reassociation.

**Still to do, and it is plant-side:** close #578/#577 on plant once
`feature/consume-leaf-package` merges — #585 already fixed both in plant's own leaf,
so the issues are stale rather than open. Tracked under issue #9.

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

**The branch is built and passing (plant #591, open), and its survey — the "8
hand edits plus a regeneration" below — was accurate about the *work* and wrong
about the *risk* in three ways. Each cost real time, and none was predicted:**

1. **A third compile break, in a dependency, presenting as eight errors inside
   somebody else's headers.** plant needed odelia at `master`, not the
   `>= 0.2.0` its DESCRIPTION asked for: plant #585 made `Patch::ode_rates`
   non-const while odelia 0.2.0's `r_ode_rates` took the system by `const&`.
   plant's `develop` did not compile against released odelia either, so the swap
   did not cause it — but nothing in plant named which odelia it wanted, so it
   arrived as template instantiation errors pointing nowhere near the cause.
   Filed traitecoevo/odelia#48; **resolved** by odelia 0.2.1 (#49) and leaf 0.1.0
   (#29), and plant now pins both ways — a `LinkingTo (>= x)` floor, which R
   checks at install, plus `Remotes: ...@sha`, which fixes what gets fetched.
   Neither subsumes the other, and the version has to be bumped *first* or the
   floor means nothing: leaf sat at 0.0.1 across #15, #24, #25 and #26, so
   `>= 0.0.1` was satisfied by a version plant cannot compile against.
2. **The pressure fix is NOT inert for plant, and it dominates everything else.**
   The survey said 10c "should be inert at 101.3 kPa, and plant's tests all use
   that". The *tests* do; `TF24_Environment` sets the **`atm_kpa` driver to
   100.5**, which nothing had checked. Measured on the one-species SCM scenario,
   deriving the ppm→Pa conversion from `atm_kpa` moves offspring production
   **+2.4%**, against **+0.10%** for the entire rest of the swap combined. TF24's
   `scientific_version` went 4 → 5 for it. **This package's golden grid evaluates
   at 101.3 and is blind to it by construction** — a standing reminder that the
   golden file bounds *this* package's behaviour, not plant's.
3. **"The `area_leaf` division cancels — do not compute it" is right for plant's
   C++ and a trap for its tests.** In plant the cancellation is exact
   (`mass_root()` is `pars.a_r1 * area_leaf`, so per-leaf-area carbon is just
   `root_mass_carbon_scale * pars.a_r1`). But `test-leaf.r` passed absolute
   carbon and `area_leaf` as separate arguments, and dropping the argument
   without dividing leaves a root system 20× too weak at `area_leaf = 0.05`. It
   compiles, runs, and moved the critical-demand collar potential from −0.685 to
   −2.57 MPa **while leaving the zero-uptake collar untouched** — that one is
   scale-invariant, so of two regression guards sitting side by side only one
   fired.

Two things from that work worth reusing:

- **The verification that licenses attributing failures to the swap.**
  `origin/develop`, built and tested in a worktree against the *same* installed
  dependencies, passes its own suite with zero failures. Without that control,
  every failure on the branch is ambiguous between the swap and the environment.
  The dry scenario gateway (`PLANT_RUN_SCENARIOS=1`) also passes against
  develop's blessed baseline, so despite the +2.4% every scenario's
  success/failure classification is unchanged.
- **"Checked — there is nothing upstream to port" has a shelf life of days on an
  actively developed sibling.** Both this document and the handoff notes were
  caught by it within 24 hours: plant merged #585 into `develop`, which fixed
  #577 *and* two further stale-state exits, and #26 ported all of it here. The
  check is `git log <base>...origin/develop -- <the files we forked>`, and it has
  to be re-run before *acting* on the conclusion, not just before writing it
  down.

Branch `feature/consume-leaf-package` in plant does this: delete
`inst/include/plant/leaf_model.h` and `src/leaf_model.cpp`, add
`LinkingTo: leaf`, and provide a compatibility shim so that `plant::Leaf` and the
handful of leaf constants plant reads (`kg_per_mol_h2o`) keep resolving. Nothing
in plant's own sources should need to change beyond that shim.

Remaining work on that branch: build it, run the test suite, and decide whether
to keep the `plant::Leaf` alias permanently or migrate plant's ~12,000 lines of
generated RcppR6/RcppExports code to `leaf::Leaf`. The alias is much cheaper and
costs nothing at runtime.

### What #15 costs plant — surveyed, and smaller than expected

Surveyed against `feature/consume-leaf-package` after #15 merged. **plant's
hand-written C++ never touches a renamed or removed member; there is exactly one
hard break.** All 59 `active:` entries in the YAML's `Leaf` block are
`access: field` — there is no `access: member` anywhere in it.

Mechanical, compile-forced — 8 hand edits plus a regeneration:

| file:line | edit |
|---|---|
| `inst/RcppR6_classes.yml:64,68,69,73` | delete the `sapwood_volume_per_leaf_area_`, `area_leaf_`, `rho_`, `a_bio_` bindings |
| `inst/RcppR6_classes.yml:118` | `set_physiology` 14 args → 10 |
| `src/tf24_strategy.cpp:426` | build per-leaf-area carbon directly |
| `src/tf24_strategy.cpp:445` | the new 10-argument call (the *only* hand-written C++ break) |
| `src/tf24_strategy.cpp:383` | delete the now-dead `sapwood_volume_per_leaf_area` |
| — | `make RcppR6` → 5 generated files + `man/Leaf.Rd` |

Three things the survey corrected, each of which would have cost time:

- **`rho_` and `a_bio_` are also bound, and also removed.** They were on nobody's
  list. On `master` all three of `sapwood_volume_per_leaf_area_`, `rho_` and
  `a_bio_` were pure dead stores — assigned, never read.
- **`c_r_V_` / `c_r_H_` need no edit at all.** They already route through
  `name_cpp: "roots_.network_.c_r_V"`, and that struct is byte-identical across the
  supply-path refactor. This is hazard 7's dotted-path trick already paying off.
- **The `area_leaf` division cancels exactly — do not compute it.** `mass_root()` is
  strictly linear in `area_leaf` (`pars.a_r1 * area_leaf`), so the vector is
  `root_mass_carbon_scale * pars.a_r1 * (prev_q - q)` with `area_leaf_` gone
  entirely. Multiplying by `area_leaf_` and dividing back is algebraically identical
  but **not** bit-identical, and would jitter golden values for nothing.

**`src/tf24_strategy.cpp:46` (`leaf.soil_consumption_[soil_layer] * area_leaf_`)
stays correct**, which is worth stating because it looks like it should not:
`soil_consumption_` keeps its per-unit-leaf-area meaning under the change.

High-volume but dull: `tests/testthat/test-leaf.r` has ~38 `set_physiology` calls,
~26 `Leaf(...)` calls, ~10 removed-field assertions, and one exact error-string
match (`mass_root_prop` → `root_carbon_per_leaf_area`); `test-pm-leaf-demo.R:47-51`
has one more. A helper wrapping `set_physiology` would collapse most of it — the
file is currently begging for one.

**Where the actual risk is: plant compiles clean and the numbers move.** The loud
failures are all mechanical. These two are not:

- **The ppm-to-Pa change must be confirmed against plant's reference outputs, not
  taken on trust.** It should be inert at 101.3 kPa, and plant's tests all use that,
  but any scenario off sea level genuinely moves. Check `inst/scenarios/`.
- **The shutdown fix changes dry-margin water balance, and that shift is the fix.**
  Any plant regression baking in the old stale numbers will fail and needs
  regenerating *with a note saying why*. Run a **long, dry** SCM scenario — the
  5-year test scenario may never enter the shutdown branch at all, which is the
  outstanding item from item 1.

Two optional and deliberately deferred: renaming plant's own
`pars.b`/`pars.c`/`pars.g1_TF24` (`tf24_strategy.h:76,77,81`) — plant has precisely
the ambiguity hazard 1 is about, with the stem pair sitting 300 lines from
`root_b`/`root_c`, but renaming is a **user-facing break** to plant's strategy
parameter names, propagating to hyperpar columns, `test-strategy-tf24*.R` and
`inst/scenarios/*.csv`; and exposing the thirteen newly-settable
temperature-response constants, which is the whole point of having made them
settable (TF24t acclimation) but is purely additive.

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

## 4. Drop the last R coupling — DONE

**Fixed upstream and here.** odelia's `ode_util.hpp` included `RcppCommon.h` for
`Rcpp::stop`/`Rcpp::warning` and for `Rcpp::as`/`wrap` declarations on
`odelia::util::index`, which made it the one thing standing between this package
and plain C++. Filed as traitecoevo/odelia#43, fixed in traitecoevo/odelia#44.

odelia took the "move `stop`/`warning` to a throw" route rather than an
`ODELIA_NO_R` guard: a guard would have left R as the default, so plain-C++
consumers would still have needed to know to define the macro. `util::stop` now
throws `std::runtime_error` exactly the way `leaf/util.hpp` does, `util::warning`
writes to stderr, and the `util::index` / `as` / `wrap` block turned out to be
dead code — declared, never defined, no callers — and was deleted.

`tests/cpp/shim/RcppCommon.h` is gone. The suite builds against odelia's real
headers with nothing standing in for R, and reproduces all 288 golden points
bit-identically, so this was a clean removal rather than a swap. Note that
odelia now runs a `tests/standalone/` check in CI on a runner with no R, which
is what keeps this from silently regressing on us.

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

Three decisions and a staged plan follow: 6a settles whether the package stops
being header-only, 6b what to build the bindings against, 6c whether to use a
generator at all, and 6d is the five stages themselves and where they have got
to.

### 6a. "Header-only" — DECIDED: two layers, and the headers keep the guarantee

Three documents asserted that this package ships no compiled code
(`DESCRIPTION`'s Description field, `NAMESPACE`'s comment, and by omission
`.github/workflows/cpp-tests.yml`, which builds with no R present). RcppR6 adds
`src/`, `R/` and a shared library, so one of those two facts has to give.

**The decision: the package gains a compiled R layer, and `inst/include/` keeps
every guarantee it has today.** Two layers with a one-way dependency:

| layer | contents | who consumes it | rule |
|---|---|---|---|
| `inst/include/` | the model | `LinkingTo: leaf` (plant), and anyone embedding it in C++ | plain C++. No Rcpp, no `R.h`, no R at all. Never includes anything from `src/` |
| `src/`, `R/` | RcppR6-generated glue plus a hand-written R layer | R users, via `library(leaf)` | may use Rcpp freely. Includes *downward* into `inst/include/` and never the reverse |

This is what makes the change compatible with item 4 (#11), which only just got
the include graph R-free. **#11's guarantee is about the include graph reachable
from `leaf.hpp`, not about the tarball containing no `src/`** — those were the
same statement while there was no R layer, and they stop being the same
statement here. The one that matters to a consumer is the first: plant compiles
`leaf.hpp` into its own translation units, and what it must not inherit is
Rcpp, not the existence of a `.so` it never links against.

**The invariant that replaces "header-only", and how it stays true.**
`cpp-tests.yml` builds and runs the whole C++ suite with **no R installed on the
runner** — that is the assertion, and it is already written and already green.
It keeps working unchanged, because it compiles `tests/cpp/*.cpp` against
`inst/include/` and never looks at `src/`. So the day someone reaches for
`Rcpp::stop` inside a header, three CI jobs go red before review. Adding the R
job *alongside* it rather than replacing it is the whole point; a single
`R CMD check` job would compile the headers with R present and could not tell
the difference.

**What this costs, stated rather than discovered later:**

- The `R CMD check` NOTE `'LinkingTo' field is unused: package has no 'src'
  directory` disappears. `.claude/CLAUDE.md` says to expect exactly that NOTE;
  that sentence becomes wrong and must move with the code.
- `DESCRIPTION` gains `Imports: R6, Rcpp`, `LinkingTo: Rcpp`, and loses the
  sentence "The package ships headers only: there is no compiled code and
  nothing to link against." `NAMESPACE` gains `useDynLib` and stops being a
  comment block.
- Installing the package now needs a compiler. It always effectively did — a
  `LinkingTo: leaf` consumer compiles these headers, and BH/odelia/Rcpp all
  require a toolchain — but "needs a compiler" moves from the consumer's build
  to ours, and a binary-only R installation can no longer install it from
  source. `tests/cpp.R` already handles the no-toolchain case by skipping; the
  install itself cannot.
- `tests/cpp/` and its Makefile stay exactly as they are. The C++ suite remains
  the regression baseline, and stage 1 adds an R-side test *in addition to* it
  rather than porting it.

**The alternative that was rejected: a separate `leafr` package** holding the R
layer, leaving this one pure. It preserves the current framing exactly and it is
the wrong trade. Two repos to keep in step for one model, and the family has
just finished paying for that failure mode in the other direction — see the
pinning discipline under item 3, where a sibling's version not moving across
four merges made `>= 0.0.1` meaningless. It also puts a choice in front of
precisely the audience item 6 exists for: someone who wants what `plantecophys`
gives them should type `install.packages` once, not first work out which of two
packages is the one with the functions in it.

### 6b. What to build the bindings against — plant's, moved rather than rewritten

plant already has them, and they are the reference implementation: `Leaf` moves
out of plant's `inst/RcppR6_classes.yml` and is generated against this package.
That is deliberate rather than lazy. plant's version is the surface its
`test-leaf.r` has been exercising for years, so translating it means a failure in
that stage is a translation error and nothing else — the redesign happens in a
separate, later stage where it cannot be confused with one. plant's copy then
goes, which is the hazard-7 payoff (stage 4).

Two design notes carried from item 7b-iii:

- **Wire up `supply_kind_` / `single_`** so `SinglePotential` is reachable. 7b-iii
  flags the footgun: a naive setter silently invalidates a configured root
  network. Design it, don't just expose it. (Stage 3.)
- **Absorb the `Control` struct** (item 10b): `GSS_tol_abs`, `ci_abs_tol`,
  `ci_niter`, `vulnerability_curve_ncontrol`, `integration_tol_`,
  `leaf_temp_min`/`_max`. Deferred *to* this item because it means surgery on a
  19-argument constructor that item 6 redesigns anyway. (Stage 2.)

### 6c. RcppR6 versus odelia's hand-written bindings — DECIDED: both, in layers

odelia solves the same problem differently, and it is worth saying why this
package does not simply copy it. odelia is **also** Rcpp + R6 — the difference is
that its glue is hand-written: `// [[Rcpp::export]]` free functions taking an
`Rcpp::XPtr`, wrapped in an R6 class that holds the pointer. No generator.

The two are not actually in competition here, because **stage 2 writes a
hand-written R layer either way.** Named arguments, defaults, a `Control` object
and a one-call entry point are design work, and no generator produces them. So
the only thing in question is the ~90 low-level accessors underneath: 60 field
getters and setters and 30 method forwards.

**Generate those.** They contain no design content, and hand-writing them is
90 opportunities to transpose an argument — which is precisely the error class
that produces plausible numbers and a green C++ suite, and precisely why the
R-side golden test had to be written. Generation also keeps stage 4 a *move*
rather than a rewrite: plant is an RcppR6 package, plant's `Leaf` block is the
same YAML dialect, and plant's `test-leaf.r` asserts the R-side names that
`access: field` preserves. And the dotted `name_cpp` trick — `roots_.psi_soil_`
reaching a moved member without the R name changing — is what makes hazard 7
cheap; it has no hand-written equivalent that is not just as much code.

**The real cost is that RcppR6 is not on CRAN**, and that is worth defusing
rather than living with. It is defused by **not declaring it at all**: the
generated files are committed, so nobody installing or checking this package
needs RcppR6 — only a developer regenerating does. It is therefore absent from
`Suggests` and from `Remotes`, and installed explicitly by the one CI job that
regenerates and diffs the output. That job exists because a committed generated
file can go stale silently, which is the one genuine downside of generating.

**If this package ever does need to drop RcppR6** — CRAN submission is the
plausible trigger — odelia's approach is the migration target, and the migration
is cheap *because* of the layering above: the public R API lives in the
hand-written layer, so the generator underneath can be replaced without the
surface moving. That is a reason to build the layer even while RcppR6 stays.

### 6d. The staged plan — 0, 1 and 2 are done; 3 is next; 2b and 4 are blocked

Five stages, each landing on its own, ordered so that the two hard-to-reverse
decisions — whether the package stops being header-only, and what the R-facing
names are — come first and cheapest, and so that no stage leaves `master` unable
to build plant.

**Stage 0 — settle header-only. DONE**, as 6a above. Written down before any
YAML, which is the point: it is what a reviewer would otherwise argue about in
the middle of a 2000-line diff.

**Stage 1 — the YAML and the generated glue, no design changes. DONE.** `Leaf`
came out of plant's YAML minus the four fields #15 deleted, with
`root_collar_psi_` → `opt_root_psi_` and the constructor's R-side names corrected
to `stem_c` / `stem_b` / `cost_scale_TF24`. λ and `g1_eff` were added, since they
are what #5 asks for and were pure additions of existing C++ accessors.
`R CMD check` is Status OK with zero NOTEs.

⚠️ **The C++ suite is the regression baseline and is blind to the R layer**, so a
mistranslated argument gives a green suite and plausible R numbers.
`tests/testthat/test-golden.R` ties four golden operating points back through the
R API, compared bit-exactly. Four things that were not in the plan:

1. **RcppR6 forces the `.h` / `.hpp` split.** It hardwires `#include <leaf.h>` in
   what it generates and `<leaf/RcppR6_pre.hpp>` in the umbrella it expects, so
   the R-facing umbrella *must* be `leaf.h` and the model's had to stay
   `leaf.hpp`. A file extension now carries the one-way-dependency rule, which is
   more weight than an extension should bear — hence the banner in `leaf.h`.
2. **The generated `RcppR6_*.hpp` live in `inst/include/leaf/` and do include
   Rcpp**, so "no Rcpp under `inst/include/`" is no longer the invariant. The
   invariant is *directional*; hazard 9 in the developer guide was rewritten.
3. **The R test's expected values must be C99 hex floats.** Decimal `%.17g`
   strings would fail for ~18% of them against a model that is exactly right —
   issue #13's parser problem, arriving from a new direction.
4. **The golden file's bit-exactness is conditional on the optimisation level**,
   not only the platform: identical at `-O1`/`-O2`/`-O3`, off by 3.47e-15 at
   `-O0`, which does not contract `a*b + c` into an FMA. Found because CMake's
   default build type is flagless.

**Stage 2 — a friendlier surface, as a hand-written R layer over stage 1. DONE.**
Kept *above* the generated glue so it can be rewritten without touching CI or the
build — which is also what makes 6c's generator replaceable.

- `leaf_traits()` / `leaf_control()` split the constructor's 19 arguments, four of
  which are tolerances. A test asserts the two partition it exactly.
- `leaf_model()`, `set_drivers()`, `operating_point()` — named and defaulted.
- `leaf_solve()` — drivers in, operating point out as a data.frame, vectorised.
- `vignette("leaf")`.

Two decisions worth keeping: **the `psi_soil >= 0` rejection is surfaced, not
smoothed over** — a friendly wrapper is exactly where someone would `abs()` it,
and that check is all that stands between a pre-#25 script and a plausible wrong
number. And **the leaf-temperature clamp is NOT in `leaf_control()`** though 10b
listed it: it is a guard keeping the Arrhenius block finite on the energy-balance
path, not a tolerance anyone tunes, and making it settable would let a caller get
NaNs back with no indication why.

Writing the vignette also found **#31**: `profit_psi_stem_TF` returns a plausible
number for `psi_stem < psi_upstream`, built on a negative conductance, with a
discontinuity of 1.58 units at the boundary. Unreachable from plant's solve,
which is why it survived; the first thing an R user does is plot the profit
function.

**Stage 2b — take resistances, not root carbon. BLOCKED on plant #591**, since it
is coupled. It belongs here because it is the same surgery on the same argument
list.

`set_physiology` takes `root_carbon_per_leaf_area` and calls
`roots_.set_root_network_from_carbon(...)` internally. **The solve never touches
root carbon.** `uptake_impl` and `duptake_dpsi` read exactly two vectors —
`network_.r_R_H_min` and `network_.r_R_V_sum` — plus `grav_head_z_` and
`max_soil_layer`. Everything else the carbon path produces (`c_r_V`, `c_r_H`,
`r_R_V`) is carried only because plant exposes it through RcppR6.

So taking carbon makes the leaf own four things that are not gas exchange:
`beta_R_H`, `beta_R_V`, `dz_`, and the 1/3 : 2/3 horizontal/vertical split. Taking
resistances instead deletes all four from the leaf's surface, and it is exactly
the move `leaf_specific_conductance_max` already makes: plant computes
`kmax = K_s*theta/(h*eta_c)` and hands over a scalar, because which
conductance-versus-height model is in force is not this package's business. Which
root-architecture model is in force is not either. `MultiLayerRoots::set_root_network`
already exists and takes a `RootNetwork` directly.

Two objections, both answerable:

- **"It leaves the tested surface."** `roots.hpp` keeps `root_network_from_carbon`
  here "so that it remains covered by the golden file — the moment this arithmetic
  crosses the package boundary it leaves the tested surface." Keep the function, as
  a documented public helper; what moves out is the *call*, not the model. The
  golden grid then calls helper-then-`set_physiology` and covers both.
- **"It costs an allocation per solve."** Only if the caller rebuilds the vectors.
  Building five fresh ones each call measured **+0.074 µs (0.061 → 0.135), about
  +2% of a whole solve**, which is why the in-place `RootNetwork& out` overload
  exists. plant already holds the carbon buffer as a `TF24_Strategy` member and
  refills it; it would hold a `RootNetwork` the same way. Neutral done that way —
  and worth re-measuring interleaved, per hazard 5, rather than assumed.

The gain for an R user is the real argument, and stage 2 sharpened it: a bare-leaf
caller has no root carbon profile, no layer thickness and no opinion about
`beta_R_V`, so `set_drivers()` currently invents a default for them and says so in
its own documentation. That is the smell. They have a resistance, or they have
`SinglePotential`.

**Stage 3 — expose `SinglePotential`. DONE (#32).** `leaf_supply_single(resistance,
gravity_head)` and `leaf_supply_multilayer()`, chosen through `leaf_model(supply =)`
and `leaf_solve(supply =)`. That is the bare-leaf, one-ψ_soil use case that was
reason 2 in 7b for building it, and it is what makes the 7a model comparison
possible: Medlyn, Prentice least-cost and Cowan-Farquhar are all formulated
against a *single* soil potential.

7b-iii's footgun — **a naive `supply_kind_` setter silently invalidates a
configured root network** — is designed out rather than documented. There is no
settable tag at any level: C++ gained `set_supply_multilayer()` /
`set_supply_single()`, each of which reconfigures the object completely and
clears the solved state, and `supply_kind` / `single_resistance_` /
`single_gravity_head_` are bound **read-only**. There is no intermediate state in
which the tag and the supply disagree, and a test asserts the assignments fail.

Two smaller things it forced, both worth keeping:

- `setup_clean_leaf()` now clears **both** supply paths' soil state, not just the
  active one. A Leaf switched between paths is precisely the case where the
  inactive path's stale state could come back — hazard 8, from a new direction.
- `set_drivers()` **refuses** `soil_depth` / `root_carbon_per_leaf_area` on the
  single path rather than ignoring them. Silently ignoring a profile someone took
  the trouble to pass is how a plausible wrong number gets made.

**Stage 4 — delete plant's copy, and collect the hazard-7 payoff. BLOCKED on
plant #591.** Once the bindings live here, plant's `inst/RcppR6_classes.yml` stops
naming `Leaf`'s fields and hazard 7 dissolves: no future member move needs a
coupled commit in two repos. A coupled PR pair, and last because it is the only
stage that can break plant.

**Sequencing.** Stages 0–3 leave plant untouched and can land while plant #591 is
in review. 2b and 4 are coupled and must wait for it to merge or they rebase
against a moving target. **2b before 4**, since 4 deletes the bindings 2b changes.

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

⚠️ **This preference has since been overturned — see 7b-iii stage 2.** Two things
this table does not say: `std::variant` is the *slowest* of the three options
here, and the numbers were taken against the pre-stage-1 structure. Re-measured
after stage 1, with a real second alternative, an **enum tag + `switch` is free**
(−0.8%, i.e. at or below a direct call) while `std::variant` costs +1.0% — a
third of the 2.6% below, but not nothing. The *conclusion* that survives from
this section is the constraint — `Leaf` must stay copyable, and no template —
not the ranking.

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

### 7b-i. What actually moves — read this before starting

Inventoried against the code so the refactor starts from a list rather than a
guess.

**State that moves into `MultiLayerRoots`** (line numbers as of `25b6599`):

| group | members |
|---|---|
| root traits and curves | `root_b`, `root_c`, `root_psi_crit`, `root_vuln_from_psi`, `root_vuln_integral_from_psi` |
| soil geometry | `soil_number_of_depths_`, `max_soil_layer`, `soil_depth_`, `z_soil_mid_`, `use_precomputed_z_soil_mid_`, `dz_`, `grav_head_z_` |
| soil state | `psi_soil_` (positive magnitudes; `psi_soil_inverted_` deleted in #25) |
| resistance network | `r_R_H_min`, `r_R_V`, `r_R_V_sum` |
| per-solve cache | `root_vuln_integral_soil_` |
| outputs | `soil_consumption_`, and the `E_up_` it accumulates |

**Functions that move:** `E_from_Soil_to_Root_Collar`, `dE_from_soil_dpsi_collar`,
`setup_root_vulnerability`, the root-network block of `set_physiology`
(`r_R_H_min` / `r_R_V` / `r_R_V_sum` and the `z_soil_mid_` fallback), and the
soil-side cache build at the top of `prepare_collar_solve`.
`build_cumulative_vulnerability_integral` is **shared with the stem** — leave it
on `Leaf` or make it a free function; do not drag it across.

**Functions that stay on `Leaf`,** because they span both sides: `E_column`,
`E_column_zero`, `find_root_psi`, `find_psi_stem_from_psi_root`,
`prepare_collar_solve`, `find_root_collar_psi`.

### 7b-ii. Four things the proposed interface does not yet account for

The `SupplyPath` sketch above is close but incomplete. Each of these will bite
mid-refactor if it is not decided up front.

1. **`uptake()` is not a pure function, and plant writes back into its output.**
   It fills `soil_consumption_[i]` per layer, which feeds plant's patch water
   balance. Checked on plant's side, and it is worse than read-only: after solving
   at each height, `tf24_strategy.cpp:505-508` **assigns** the crown-integrated
   value straight back into `leaf.soil_consumption_[a]` (and `leaf.E_up_`,
   `leaf.profit_`, …), reading it at `:46` as `leaf.soil_consumption_[soil_layer]`.

   So `soil_consumption_` is not simply an output the supply path owns — it is a
   buffer plant reaches into by name and overwrites. **Recommendation: leave
   `soil_consumption_` and `E_up_` as members of `Leaf`** and have the supply
   path write into a buffer handed to it by reference. Moving them into
   `MultiLayerRoots` breaks plant's access path for no benefit, and would also
   entangle the crown write-back with the supply path's per-solve scratch state.

   Note the deliberate unit split while touching this: `E_up_` is kg H₂O m⁻² s⁻¹,
   `soil_consumption_[i]` is **mol**, converted downstream in plant.
2. **`wettest_potential()` needs a per-solve entry point.** Today the soil-side
   caches (`psi_soil_inverted_`, `root_vuln_integral_soil_`) and the wettest
   layer are built in one pass at the top of `prepare_collar_solve`. Give the
   concept a `begin_solve()` that does both and returns the wettest potential —
   otherwise the cache is either lost or rebuilt per call, and it is a measured
   hot-path optimisation (it collapses ~2 spline evals per layer to ~1).
3. **The cache fast path is selected by pointer identity.**
   `E_from_Soil_to_Root_Collar` tests `&psi_soil == &psi_soil_inverted_` to decide
   whether the cache is valid. Once that vector lives in another object the test
   still compiles and silently changes meaning. This is the single most likely
   place to lose bit-identity — the golden file will catch it, but know where to
   look.
4. **`duptake_dpsi` returning NaN is a contract, not a failure.**
   `dE_from_soil_dpsi_collar` deliberately returns NaN at the three branch kinks
   (equal potentials, gravity balance, the ψ=0 split) so the caller falls back to
   finite differences. Preserve that, and document it on the concept — a naive
   implementation that throws or returns 0 would silently degrade TF24f's
   acclimation gradient.

### 7b-iii. Suggested staging

Each stage is checkable bit-exactly against the golden file, which is the point
of doing it in stages.

1. **Pure move, no interface — DONE.** `MultiLayerRoots` as a plain member held by
   value, `Leaf` forwarding to it. No variant, no virtual. Golden bit-identical
   over all 288 points, 105/105 unit checks. What the stage actually taught, in
   descending order of how much it will matter later:

   * **A fifth trap, and the expensive one: moving a public member breaks
     plant's *generated* RcppR6 glue.** 7b-i's inventory was written against this
     package only. plant binds eleven of the moved fields with `access: field`,
     which emits `obj_->psi_soil_` as both getter and setter into
     `src/RcppR6.cpp`. The fix is one YAML line each
     (`name_cpp: "roots_.psi_soil_"` — the template pastes it verbatim after
     `->`, and the R-side name comes from the YAML key so it does not move), plus
     two lines in `tf24_strategy.cpp` for the `z_soil_mid_` write. Cheap, but it
     means **stage 4 is not optional and not last**: plant tracks this package's
     `master` via `Remotes:`, so the merge is what breaks it. Land them together.
   * **Trap 3 (pointer identity) was survivable, and is now half-defused.**
     `&psi_soil == &psi_soil_inverted_` keeps its exact meaning after the move,
     because the member and the test moved together and callers still pass the
     member by reference. The hot path now goes through `roots_.uptake()`, which
     takes an explicit cache flag; the identity test survives only in
     `uptake_at()`, the arbitrary-vector entry point, so that the R-facing
     `Leaf::E_from_Soil_to_Root_Collar` cannot change behaviour. It disappears
     for good in stage 2, when `E_column` / `find_root_psi` stop threading a
     `psi_soil` vector through at all.
   * **Trap 1 (`soil_consumption_` / `E_up_`) held.** They stayed on `Leaf` and
     are handed to `uptake` by reference, as recommended. No friction.
   * **It costs 1.7%** (3.53 vs 3.47 µs/solve, interleaved ×3 at reps=2000),
     with no dispatch added yet. Not attributed: `uptake_impl` is still
     out-of-line, as `E_from_Soil_to_Root_Collar` was, but its argument list went
     from 2 (+`this`) to 6, and it is called ~10³ times per solve. That is the
     leading candidate, not a measured cause. One hypothesis was tested and
     **refuted**: accumulating `E_up` into a local to break the assumed aliasing
     with the `soil_consumption` buffer changed nothing, so it was reverted
     rather than kept on a story it did not earn.
   * **Budget note for stage 2.** The +2.6% predicted for `std::variant` was
     measured against the *old* structure. If it stacks, stage 2 lands around
     +4.3% total. That is still inside the band the closed-form path (6.3×/27×)
     dwarfs, but say the number rather than discovering it.

   **Stage 4 was done alongside, not deferred**, for the reason above. plant's
   `feature/consume-leaf-package` carries the eleven YAML lines and the two in
   `tf24_strategy.cpp` (`3efe9c47`), plus five re-pointed paths once the supply
   path started taking resistances (`bbd47a36`). Run **with a control**, both from
   a clean `src/`:

   | build | result |
   |---|---|
   | leaf `master` + plant unchanged | 2431 pass / 0 fail / 7 skip |
   | this branch + plant's YAML change | 2431 pass / 0 fail / 7 skip |

   Identical, including the skip list. Note this is **2431, not the 2364 recorded
   under item 1** — the control shows that gap is environmental (`NOT_CRAN`
   unset, so seven tests skip) and predates this work. The earlier figure should
   not be read as a target to match.

   The control was worth its cost for a second reason: `R CMD INSTALL` does not
   clean, and a header-only `LinkingTo` dependency changing underneath a stale
   `.o` moves no `.cpp` timestamp, so nothing rebuilds. That produced a load
   failure naming a field accessor nobody had touched. Written up under "Build &
   regeneration workflow" in plant's `agents.md`, with the `nm` check that
   identifies it.
2. **Introduce the concept — DONE**, together with stage 3. Golden bit-identical
   over all 288 points; 134 checks; the dispatch measured **free** (3.537 vs 3.563
   µs/solve against the pre-dispatch build, at or below it in every round).

   `Leaf` holds **both** alternatives as members and selects with a
   `SupplyKind` enum. Default `MultiLayer`, so every existing caller — plant
   included — keeps today's behaviour without knowing this exists.

   **The three R-facing signatures were saved, and cheaply.** The obstacle was
   that `E_column`, `find_root_psi` and `find_psi_stem_from_psi_root` thread
   "the current soil state, signed" as a `std::vector<double>`, and three of those
   are RcppR6-exposed. Rather than removing the threading — an API change that
   would have forced this onto `feature/api-cleanup` — `SinglePotential` simply
   carries a one-element vector of its own, and `Leaf::supply_psi_soil_inverted()`
   returns whichever is active. One duplicated double, and the entire solve stayed
   supply-agnostic without a signature moving. **This also means stage 4 is a
   no-op here**: no plant-visible name moved, so plant needs no YAML change.

   Dispatch is needed at only five points, which is the measure of how well the
   boundary was drawn: `begin_solve`, `uptake`, `duptake_dpsi`, the collar-potential
   bound (`root_psi_crit` for roots, `psi_crit` for a constant-conductance path),
   and the layer count that sizes `soil_consumption_`.

   **Still open, deliberately.** Two things:

   * The `psi_soil` threading survives, so the pointer-identity cache test in
     `uptake_at` survives with it. Removing both is the API change described
     above and belongs with 10b. Nothing forces it now that the concept works
     without it.
   * **`supply_kind_` is not reachable from R.** Switching supply paths is a C++
     API only; plant's YAML exposes neither it nor `single_`. That is the right
     default — plant always wants the multi-layer path, and exposing a setter
     that silently invalidates a configured root network is a footgun — but it
     means the bare-leaf, one-ψ_soil use case that motivated `SinglePotential`
     (reason 2 in 7b) is still unreachable for an R user. Wire it up with item 6,
     where the whole R interface gets designed, not before.

   **The supply path now takes resistances, not root carbon — done ahead of this
   stage, because the concept depends on it.** `MultiLayerRoots::set_root_network`
   takes a `RootNetwork` (`r_R_H_min`, `r_R_V_sum`, plus three diagnostics), and
   the carbon → resistance map is the free function `root_network_from_carbon`.

   The reason is that the solve reads **exactly two** of those vectors, plus
   `grav_head_z_` and `max_soil_layer`. Nothing in `uptake` or `duptake_dpsi`
   touches root carbon, the 1/3 : 2/3 split, `dz`, or either `beta_R_*` — those
   are inputs to a *root architecture* model that happens to run just before.
   This resolves a contradiction that had been sitting in this document: **10b
   lists "the whole root architecture" under "move to the plant side", while 7b
   above says "do not push it up into plant".** Both are right, and they meet at
   the resistance interface — the architecture is plant's, the transport solve is
   ours. It is the same split `leaf_specific_conductance_max` already makes, which
   10b names as "the model to copy".

   Three things learned doing it:

   * **The map stays in this package, as a free function.** Moving the arithmetic
     into `tf24_strategy.cpp` would land it beside `root_mass_carbon_scale` and
     `rooting_depth_max` — already TODO-flagged as magic numbers — and, worse,
     take it outside the golden file's reach. It is now covered *twice*: through
     `set_physiology` as before, and by a direct unit test, which is the payoff
     for it being a free function at all.
   * **Building a fresh `RootNetwork` per call cost +0.074 µs on
     `set_physiology`** (0.061 → 0.135 µs, ~+2% of a whole solve): five vector
     allocations replacing in-place resizes that had been reusing capacity.
     `std::move` did not help — the allocation, not the copy, was the cost. Fixed
     by having `MultiLayerRoots` hold the `RootNetwork` and the map fill it in
     place; the value-returning overload remains for tests and one-off callers.
     Measured back to 0.056–0.061 µs, and the solve itself unchanged.
   * **A zero-carbon layer gets `r_R_H_min = 0`** — *zero* horizontal resistance,
     i.e. infinite soil-to-root conductance where there are no roots. Backwards.
     Preserved rather than silently changed: it looks unreachable from plant
     today (`max_soil_layer` truncates at the last non-zero layer, and plant's
     `Q()` distribution will not give an exact interior zero). Passing resistances
     from plant forces the question to be answered explicitly, which is a reason
     to finish the job.

   **The other half is deferred on purpose.** Having *plant* call
   `root_network_from_carbon` and pass resistances through `set_physiology`
   changes that signature, and signature changes belonged with 10b, which reworked
   `set_physiology` (14 → 10 args, `root_carbon_per_leaf_area`) and **merged in
   #15**. So the reason for deferring is spent: this is now free to do, and doing it
   takes
   `beta_R_H`/`beta_R_V` out of `Leaf`'s constructor and makes `soil_depth`
   droppable, since `dz` and `grav_head_z_` are its only remaining consumers.
   **Which dispatch mechanism — MEASURED against the current code, and the answer
   is an enum tag.** `std::variant` was never the fastest; it was the *slowest* of
   the three options in #14, chosen for copyability. Re-measured here against the
   post-stage-1 structure, with `SinglePotential` as a genuine second alternative
   and all three arms executing the same `MultiLayerRoots` code (identical bench
   checksum, so only the *reach* varied), 10 interleaved rounds at reps=2000:

   | option | measured | vs direct | verdict |
   |---|---|---|---|
   | direct call (today) | 3.507 µs | — | baseline, not an option once there are two paths |
   | **enum tag + `switch`** | **3.478 µs** | **−0.8%** | **free.** Consistently at or below the direct call across all 10 rounds |
   | `std::variant` + `std::visit` | 3.542 µs | +1.0% | real, but a third of the +2.6% the old table claimed |
   | template `Leaf<Supply>` | — | — | still rejected: breaks the `plant::Leaf` alias and ~12,000 lines of RcppR6, for a benefit the numbers above say is under 1% |
   | `std::function` | — | — | measured +0.6% in #14, but a poor fit: this is a stateful *four-method* interface, not one callback |
   | virtual base | — | — | +1.1% in #14, plus clone boilerplate and a heap allocation per `Leaf` copy |

   **So: hold both alternatives as members and switch on a tag.** `SinglePotential`
   is four doubles, so carrying it alongside `MultiLayerRoots` costs ~32 bytes per
   `Leaf` and nothing at all in time. It is trivially copyable, needs no `visit`,
   no clone, and no union.

   *Why* it is free is worth knowing, because it is a durable property rather than
   a lucky number: `uptake_impl` is out-of-line and stays that way (hazard 5), the
   tag never changes within a solve so the branch predicts perfectly, and a
   predictable conditional in front of an out-of-line call disappears into it.
   `std::variant` cannot match that because `std::visit` emits **two out-of-line
   `__dispatch` thunks** — an indirect jump, which predicts worse than a direct
   branch. That is visible in the binary: `nm -C bench | grep __dispatch`.

   ⚠️ **Reproducing this: the tag must be runtime-unknowable.** Set it from a
   constant and the branch constant-folds, every arm measures zero, and it looks
   like a free lunch that is not one. The measurement above set it from `argc`.
   Verify before believing any arm: `nm -C bench | grep -c SinglePotential`
   should be 0 for the direct arm and non-zero for the others — if the dispatch
   folded, the unused path's symbol disappears entirely.

   **The one argument for `std::variant` that survives** is not speed: a variant
   makes the invalid state unrepresentable, whereas holding both members means
   both objects always exist and only one is meaningful. That is worth 1.0% only
   if the invalid state is a real risk. With one tag written in one place, it is
   not — but revisit if a third and fourth supply path arrive, where holding all
   of them stops being sensible.

   **The binding constraint is copyability, not speed.** plant's
   `make_strategy_ptr(TF24_Strategy s)` takes the strategy **by value** and
   `TF24_Strategy` holds a `Leaf` member, so a `std::unique_ptr<SupplyPath>`
   member breaks plant's build outright. Everything above except the virtual base
   satisfies that for free.

   **And the whole table is stale.** It was measured against the *pre-stage-1*
   structure, where the supply call was a direct member call taking two arguments.
   It now goes through `roots_.` into a six-argument out-of-line `uptake_impl`.
   Whether +2.6% stacks on stage 1's +1.7% or partly absorbs into it is an open
   question — **re-measure the shortlist against the current code before
   choosing**, using the interleaved method (hazard note under "Build and test";
   a sequential A/B got the sign wrong once already).

   **Two interface details this stage has to settle.**

   * **"Per-layer output access" is doing real work in that sentence.**
     `MultiLayerRoots` writes `soil_consumption_[i]` for every rooted layer;
     `SinglePotential` has one layer or none. plant reads that buffer *by index*
     and integrates it into the patch water balance, so the concept must express
     how many entries an implementation writes — it is not a detail that can be
     left to each implementation.
   * **Retiring the pointer-identity test is an API change.** Four `Leaf`
     functions still thread a `psi_soil` vector down to `uptake_at`: `E_column`,
     `E_column_zero`, `find_root_psi`, `find_psi_stem_from_psi_root`. Stopping
     that — so the supply object owns the current soil state and callers just say
     `uptake()` — is what finally kills the address-identity cache test. But
     `find_root_psi`, `find_psi_stem_from_psi_root` and
     `E_from_Soil_to_Root_Collar` are all RcppR6-exposed *with* that argument, so
     it changes three R-facing signatures. Results do not move, but the R API
     does. Decide deliberately which branch that belongs on rather than
     discovering it mid-refactor; see the branch rule in the developer guide.

3. **Add `SinglePotential` — DONE**, before stage 2 as the ordering warning
   required. `leaf/single_potential.hpp`: one ψ_soil, a constant series
   resistance, the same four-method contract. Golden unaffected because it is not
   the default. Nineteen checks of its own, including that its analytic
   `duptake_dpsi` matches a central difference, that a zero resistance throws
   rather than returning an infinity, and — through a whole `Leaf` — that drier
   soil and a higher series resistance both cost carbon, which is the contract
   that makes the two paths comparable at all.

   It has **no branch kinks**, so it never returns the NaN that `MultiLayerRoots`
   uses to request a finite-difference fallback. That NaN is a *contract*
   (7b-ii item 4), not an expectation: a path that never needs it is conforming,
   not broken.
4. **Check plant — DONE for every stage.** Stage 1 needed eleven YAML lines and
   two in `tf24_strategy.cpp`; the resistance change re-pointed five paths; stages
   2 and 3 needed **nothing at all**, because they only *added* members and moved
   no plant-visible name. Repeat it for any later change that moves one, since
   `master` is no longer a source-level drop-in.

Stage 1 was the only one that should have been able to break anything. That was
half wrong and half right: it *did* break plant's generated glue (trap 5), which
nobody predicted — but stage 2, which looked like it would have to change three
R-facing signatures, turned out not to need to. Letting `SinglePotential` carry
its own one-element vector of signed potentials kept the whole solve
supply-agnostic without touching a signature, and cost one duplicated double.

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
- ~~**Signed-versus-magnitude water potentials belong in the type, not a comment.**~~
  **Both halves of this were wrong, and the record is worth keeping.** The premise
  ("the convention is held together by a comment block") was right; the prescribed
  fix was not. #23 built the strong types, and they worked — bit-identical, and
  they found a live bug. The review criticism that stuck was that they *describe* a
  two-convention model instead of removing one, and that no type can cover
  `dE_from_soil_dpsi_collar`, which is neither `Psi` nor `AbsPsi`. #25 deleted the
  second convention instead: one representation, positive magnitudes, asserted at
  the input boundary. `psi_soil_inverted_` is gone, so there is nothing left to
  type. The lesson generalises: reach for a type to enforce an invariant you
  actually need, not to police one you could remove.
- **`R` and `n` are already gone** from the public namespace (see item 1). Worth
  recording *why* it mattered: the analytical project's
  `tf24_closed_form_bench.cpp` declares locals `const double n = l.c*l.beta2 - 1.0`
  and a Newton residual `R`, both inside a scope where `plant::R` and `plant::n`
  were visible. It compiles only because the locals shadow them. Removing the
  namespace-scope names makes that file strictly safer.

Do this before item 6 builds an R interface, so the R names are right the first
time, and coordinate with plant, since renames cross the shim.

### 10b. Shrink the input set

Done: the five dead items (`root_mass_`, `vcmax_25_to_jmax_25`, plus `rho`,
`a_bio`, `sapwood_volume_per_leaf_area`, taking `set_physiology` from 14 arguments
to 11), three further dead constants found later (`gamma_c`, `kc_c`, `ko_c`), and
the constants-that-are-parameters: thirteen temperature-response parameters are now
settable members, including `rd_to_vcmax_ratio_`, which had been the bare literal
`0.015` inline.

Still open, and each for a stated reason:

- ~~**`area_leaf`**~~ -- DONE. And my earlier claim that it was coupled to 7b was
  wrong: `r_R = beta/c_r` is exactly linear in root carbon, so `E_i` depends on
  `root_carbon / area_leaf` and nothing else, and passing the ratio is a four-line
  change that needs no interface. beta_R_H and beta_R_V are unchanged -- the scaling
  cancels. The leaf is now purely intensive: no extensive quantity anywhere in it.

  **The rebase onto 7b's supply interface extended this, and improved it.** 7b had
  meanwhile threaded `area_leaf` through the supply contract itself --
  `uptake`, `uptake_at`, `duptake_dpsi` all took it. Since the carbon arriving is
  already divided by it, the parameter is now gone from all three, and
  `inv_area_leaf` is gone from `uptake_impl`. So "purely intensive" is now true of
  the *supply contract* too, not just of `set_physiology`, and a third supply path
  cannot reintroduce an extensive quantity without changing the interface. One
  consequence to know: `SinglePotential::resistance_` is **per unit leaf area**,
  documented at the member.
- **The root architecture** -- that is item 7b itself.
- **A `Control` struct** for `GSS_tol_abs`, `ci_abs_tol`, `ci_niter`,
  `vulnerability_curve_ncontrol`, `integration_tol_`, `leaf_temp_min`/`_max`.
  Deferred deliberately: it means surgery on a 19-argument constructor, and item 6
  is going to redesign that argument list anyway when it builds the R interface.
  Doing it twice would be wasteful and doing it now would make item 6 harder.
- **The PM energy-balance constants** (`longwave_net_offset`, `sw_abs_per_par`,
  `latent_heat_vap`, `vol_heat_cap_air`, `aerodynamic_resistance_*`) -- these
  belong with item 13, which is going to revisit that whole block; making them
  settable before deciding what the block *is* would fix the wrong interface.


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

### 11a. Sharpen the outer solve FIRST — DONE, and it reordered the item

**Both PRs landed 2026-08-04: #35 (the sentinel) and the collar root-find.** The
decision below was written first and survived implementation; three things came out
of building it that the write-up did not predict, and they are recorded at the end
of this sub-item under "What implementation changed". Headline results:

| | golden section | root-find |
|---|---|---|
| `|dprofit|` at the returned collar, 198 interior rows | median 7.8e-04 | **median 5.6e-15** |
| rows with an improved residual | — | **240 / 240**, none worse |
| distinct argmax values over 11 trait steps | 6 / 11 | **11 / 11** |
| argmax second differences in a trait | 3.9e-04 (= the step: noise) | **3.4e-07** |
| µs/solve, interleaved at reps=2000 | 3.51 | **2.65 (24.5% faster)** |
| worst golden-file change | — | **1.5e-03** (2160 of 2592 cells, 240 rows) |



**Decision: replace the golden-section maximisation in `find_root_collar_psi` with
a safeguarded root-find on `dprofit_droot_collar_psi == 0`, as its own PR, before
any `Leaf<T>` work.** Everything below is measured through the #30 R interface at
this package's defaults, and the prototype was built in R (`uniroot` on the bound
derivative) so the conclusion did not have to be taken on trust.

Settled in writing before code, the way 6a settled header-only. Issue #4's four
comments sized the work; this sub-item records what survived checking. **Two of
its claims are corrected below and one is strengthened** — see the last two
paragraphs of the first section.

#### What was verified

The two load-bearing measurements from #4 comment 4 reproduce exactly. At
`psi_soil = 2.0`, `PPFD = 900`, `atm_vpd = 2.0`, one layer: `dprofit` at the
golden-section answer is **−6.219675e-04** against a local slope of **−8.9566
MPa⁻¹**, so the returned collar sits **6.94e-05 MPa** from the true stationary
point. Both figures are #4's to every digit quoted there.

Extended to the whole 288-point golden grid, which #4 asked for and did not have:

| | |
|---|---|
| feasible rows (GSS actually runs) | **240** of 288; the 48 at `psi_soil = 6` shut down |
| `dprofit` monotone decreasing across the bracket | **240 / 240** |
| `∂²profit/∂ψ²` at the answer — the IFT denominator | **−1.56 to −61.4**, never near zero |
| offset from the stationary point, interior rows | median **1.01e-04 MPa**, max **2.70e-04** |

So monotonicity is a property of the grid, not of one operating point, and the
IFT denominator is well conditioned everywhere — #4 comment 3's point 2 asked for
exactly this re-measurement and it holds.

**Correction 1 — the argmax is at a BOUNDARY on 42 of the 240 feasible rows**
(24 pinned at the wet end `root_zero_E`, 18 at the dry end
`min(root_crit, supply_psi_crit())`), all at `psi_soil` of 3 or 4. #4 assumed an
interior stationary point throughout. `dprofit` does not cross zero on those
rows, so **the root-find must be a constrained solve**: check the sign at both
endpoints and return the endpoint when there is no interior crossing. That is
what any safeguarded bracketing method does when handed a non-bracketing
interval, so it is a branch rather than a difficulty — but it must be written,
and it must be written *deliberately*, because the failure mode is returning a
point strictly inside a bracket whose maximum is on the edge. The remaining 198
interior rows all bracket correctly (`dprofit > 0` above the wet end, `< 0` at
the dry end).

**Correction 2 — the `0.0` sentinel fires at the WET BRACKET ENDPOINT, which
makes fixing it a hard prerequisite rather than defence in depth.** #4 flagged
the sentinel as a trap; it is worse placed than that suggests. At exactly
`bound_a = root_zero_E` uptake is zero by construction, so `psi >= psi_stem` and
`dprofit_droot_collar_psi` returns its `0.0`. A bracketing solver evaluates that
endpoint *first*, to check the bracket brackets, and reads `0.0` as "the root is
here" — returning the zero-transpiration point as the optimum. At the default
point profit there is **−1.897** against **2.516** at the true optimum. The
saving grace is that the region is narrow: bisecting on it, the sentinel extends
at most **3.46e-07 MPa** into the bracket (max over the grid; median 1.22e-08),
never more than 6.2e-07 of the bracket width. So the fix is cheap — move the
sentinel out of the searched range, to NaN or a separate feasibility flag — and
**it must land before or with the solver change, not after.**

**Strengthened — the staircase is not the worst of it. For hydraulic traits the
current solver returns a smooth, plausible, SIGN-INVERTED trait response.** This
is the finding that most changes the item, and it is why #4's instruction to
check a hydraulic trait was the right one.

#4 comment 3 measured `vcmax_25` and found the argmax piecewise constant, so a
finite difference is exactly zero below `h ≈ 0.1` and the composite `dA/dθ`
silently drops its second term. That reproduces: 6 distinct `psi_stem` values
across 11 samples, treads ~9e-4 ≈ `GSS_tol_abs`, FD exactly `0.000000e+00` for
relative `h ≤ 1e-4`, and `dA/dvcmax_25` by FD at `h = 1e-3` is **8.21e-03**
against a true **1.72e-02** — #4's number, and it is 52% low rather than merely
"biased".

But traits in the hydraulic path move the collar *bracket*, so the finite
difference does cross treads, and what it lands on is not the derivative:

| trait | naive FD (rel h 1e-6) | resolved truth | |
|---|---|---|---|
| `vcmax_25` | 0.000000e+00 | +4.008e-03 | term dropped entirely |
| `cost_scale_TF24` | 0.000000e+00 | −1.765e-01 | term dropped entirely |
| `stem_b` | +4.996e-01 | +9.257e-01 | 46% low |
| `root_b` | **−2.324e-03** | **+4.749e-03** | **wrong sign** |
| `beta_R_H` | **+1.443e-05** | **−3.510e-05** | **wrong sign** |

`d psi_stem/d θ`, central difference, default operating point. "Resolved truth"
is a least-squares slope over ±5% at n = 41, and it is corroborated below.

For `root_b` the golden-section collar decreases smoothly with the trait
(slope **−2.639e-03**, second differences 4.7e-08 — it looks *clean*) while the
truth increases (**+2.587e-04**). Wrong sign, and an order of magnitude too
large. Arbitrated three ways, which matters because a wrong gradient is
plausible: a 20001-point scan of `profit` with a parabolic refinement and **no
derivative anywhere in it** agrees with the root-find to **3–7e-10 MPa** and
gives slope **+2.5868e-04** against the root-find's **+2.5868e-04**; `profit` at
the root-find collar is `>=` the golden-section value at every sample; and
`|dprofit|` at the root-find answer is **~1e-14** against **6.2e-04** at the
golden-section answer.

That is a materially worse failure than a staircase, and it is the strongest
argument in this whole item. A staircase announces itself — a zero gradient, or
visible treads. A smooth wrong sign does not: a gradient-based calibration handed
`root_b` walks it the wrong way and either diverges or parks confidently on the
wrong value, with nothing in the trace to say so. It also cannot be repaired by
any amount of finite-difference step tuning or by differentiating through the
iterations, **because the error is in the solved argmax itself rather than in how
it is differentiated.** Only sharpening the solve fixes it. plant #406 predicted
this class of problem in the abstract — "a noisy finite-difference gradient can
push `x` the *wrong* way, not just slow it"; this is a measured instance, and it
is not noisy.

#### The prototype: it works, and it settles hazard 3 favourably

`uniroot` on `dprofit == 0` to `tol = 1e-14`, wet endpoint nudged past the
sentinel, substituted for the golden section. Compared against the current solver
at the same operating point:

| | golden section | root-find |
|---|---|---|
| distinct collar values over 11 samples at rel step 1e-3 (`vcmax_25`) | **6 / 11** | **11 / 11** |
| collar second differences, `vcmax_25` | 3.94e-04 (≈ the step: noise) | **3.44e-07** |
| collar second differences, `cost_scale_TF24` | 5.45e-04 (≈ the step: noise) | **6.75e-07** |
| FD of `d psi_stem/d vcmax_25`, rel h from **1e-8 to 1e-2** | 0 until 1e-3, then drifting | **stable to 7 digits** |
| FD of `d A/d root_b`, same range | 8.62e-03 (64% low) | **stable to 7 digits** |

The staircase is gone at source. Finite differences work across eight decades of
step size — which is also the statement that differentiating *through* the
iterations would now work, since both fail for the same structural reason and it
has been removed.

**Hazard 3 is improved, not threatened, and this is now measured rather than
argued.** The guide says plant chose golden section over Brent because the argmax
must vary *smoothly* with inputs, and that any solver change must re-measure it.
Re-measured: smoothness of the argmax in the traits improves by **~800–1000×**
(second differences 3.9e-04 → 3.4e-07). #4 predicted the direction and was right.
The honest caveat is that hazard 3's actual concern is smoothness in **plant
state** — height, light — feeding the demographic growth-rate gradient, and that
is a plant-side measurement this package cannot make while plant #591 blocks
end-to-end validation. The mechanism is the same one, so the same direction is
expected; it is not the same measurement, and the PR should say so rather than
claim the hazard discharged.

One new non-smoothness does appear and should be recorded rather than discovered
later: where a row transitions from an interior optimum to a pinned one, the
argmax has a genuine **kink** — a constrained optimum coming off its constraint.
That is a property of the constrained problem, not of the solver, and the
golden-section solver has it too. It is not a regression, but it is a real
non-differentiable set in trait space that a calibration can walk onto, and the
42 pinned rows are all at `psi_soil` of 3 to 4, i.e. the dry end where a fit is
most likely to wander.

#### Cost: unmeasured, and not to be claimed until it is

Indicative evaluation counts over the 111 interior rows sampled: a root-find to
~1e-14 takes a median of **12** `dprofit` evaluations (min 10, max 37) against
**18** profit evaluations for golden section to `GSS_tol_abs` — ratio 0.67. But a
`dprofit` evaluation costs strictly more than a profit evaluation (the same
`find_psi_stem_from_psi_root` and `ci` root-find, plus two AD sweeps, two spline
derivatives, and a possible finite-difference fallback), so **an eval-count ratio
of 0.67 is not a speedup and must not be quoted as one.** #4 comment 4's "it is
probably faster" is plausible and unverified.

The real number needs `make bench` at `reps=2000`, **interleaved**, per hazard 5
— and per the guide's own warning that a sequential A/B got the sign wrong once
while sizing #2. Two builds kept side by side, alternated three times. Do it
after the C++ lands and before the PR claims anything about speed.

#### What this means for `Leaf<T>`, and for the three AD questions

Templating is still wanted, but the reason narrows and the ordering changes.

- **The outer solve no longer needs the implicit function theorem.** At a root of
  `dprofit`, `∂profit/∂ψ = 0` holds to ~1e-14 rather than to 6e-04, so the
  envelope theorem applies cleanly and `dψ*/dθ` is one well-conditioned division
  (denominator −1.56 to −61.4 across the grid). #4 comment 3 designed the IFT
  composition as the *only* correct option; it becomes the cheap and obvious one.
- **The replica drift is a live defect and is independent of all of this.**
  `detail::assim_colimited_ad` no longer mirrors `Leaf::assim_colimited` (5.53e-4
  on the golden grid, per #4 comment 2), so the AD derivative and the forward
  model are not currently derivatives of the same function. Fix it regardless of
  what happens to the solver.
- **Item 12 may not be blocked as hard as PLAN says.** After the solver change a
  plain central difference gives ~4 correct digits at any relative step from 1e-8
  to 1e-2, which a gradient-based fit can use. That does not retire the AD
  argument — the vignette's headline result is AD *against* FD, and AD is exact —
  but "item 12 cannot start until #4 lands" should be re-read as "cannot start
  until 11a lands", which is a much smaller gate.

On the three questions raised when this was handed over:

- **Derivatives across many parameters.** Forward mode costs one sweep per
  parameter, which for the 10–15 traits a calibration fits is 10–15× a solve.
  Vendored XAD already supports **vector forward mode** — `xad::fwd<T, N>` seeds
  N directions and reads N sensitivities from one sweep, N fixed at compile time
  (`XAD/Interface.hpp:59`, `XAD/Literals.hpp:49`). **Nothing in leaf, plant or
  odelia uses N > 1 today**, and there is a trap: `ExprTraits<FReal<Scalar,N>>::
  vector_size` is hard-coded to 1 while the `FRealDirect` variant reports N, so
  `xad::fwdd` is the safer carrier if `fwd` misbehaves at N > 1. Prefer this over
  reverse mode here: PLAN item 5's settled decision is that only forward mode is
  used so **nothing needs linking**, and a tape in leaf would spend that. Reverse
  mode is the right tool at plant scale, where `ff16_production_kernel.h` already
  demonstrates trait gradients through a templated kernel in one adjoint sweep.
- **Not carrying nested solvers through the tape.** This is plant #537's rule,
  and it names `dprofit_droot_collar_psi` a deliberate **cut-point** — hard-coded
  to `xad::fwd<double>`, taking and returning plain `double`, therefore opaque to
  an outer AD pass and uncompilable against a future `Leaf<T>`. Its design rule:
  prefer derivative code templated on the scalar type, and where a hand
  derivative is kept for speed, register it as a custom adjoint / IFT checkpoint
  the outer tape can traverse, with a templated reference path alongside. That
  rule, not the replica deletion, is the real content of "template `Leaf`" — and
  the solver change **reduces** the number of cut-points needed, because the
  outer optimisation stops being one. odelia already relies on the composition
  property that makes this safe: forward tangents need no tape, so they never
  contend with an outer adjoint (`odelia/ode_jacobian.hpp:11`).
- **TF24f as a way to avoid nested solvers — and the thing nobody has written
  down.** TF24f makes the collar potential a sixth ODE state with
  `dψ/dt = k_acclim · ∂profit/∂ψ` (`plant/src/tf24f_strategy.cpp:26-38`), chosen
  as gradient ascent precisely because the relaxation form `dx/dt = (x*−x)/τ`
  would need the argmax explicitly — the search being deleted (plant #525). So
  **TF24f's fixed point is exactly `dprofit == 0`, which is exactly what 11a
  proposes to solve.** They are one equation approached from two directions:
  TF24f relaxes onto it over demographic time, the root-find lands on it at once.

  The consequence is a correctness argument for 11a that is not in #4 at all.
  plant #529 notes that TF24f's `k → ∞` limit is the quasi-steady state
  `∂profit/∂ψ = 0`, i.e. "solving the optimisation that TF24 already does each
  step". That is not currently true: TF24's golden-section answer misses
  `dprofit = 0` by 6.2e-04, so TF24 and TF24f's fixed point are **different
  points**, ~7e-05 MPa apart, and no choice of `k` makes them agree. After 11a
  they are the same point. TF24f does *not*, however, help item 12 — calibrating
  against steady-state gas exchange needs the fixed point, which is the
  root-find. TF24f buys the demographic path, not the leaf-level fit.

  Worth noting for whoever picks up plant: TF24f keeps the two bracket root-finds
  (`root_crit`, `root_zero_E`) and the `ci` root-find on every step, and it calls
  `find_root_collar_psi` once per individual at birth to seed the state
  (`plant/src/tf24f_strategy.cpp:51`). So it deletes the golden-section layer, not
  the nesting.

#### Sequencing

1. **The sentinel fix**, on its own — `dprofit_droot_collar_psi`'s `0.0` out of
   the searched range. Separable, small, and golden-bit-identical (no golden
   column comes from the gradient; hazard 8 says the same of the last three
   fixes). Land it first.
2. **The solver change**, as a deliberate results-moving PR with its blast radius
   measured against the golden file and split by cause, the way #15 was done. It
   moves every output by roughly the argmax correction, ~1e-04 MPa on the collar,
   which is at this guide's "real difference" threshold — and it is a
   **correctness improvement**, checkable by the residual: `|dprofit|` at the
   returned point goes from 6.2e-04 to ~1e-14, and `profit` does not decrease at
   any sampled point. Include the constrained-optimum branch, the interleaved
   bench, and the honest note that plant-side smoothness is unmeasured.
3. **Then reassess `Leaf<T>`** against what is left: the replica drift fix, the
   trait partials, and #537's cut-point rule. The scope should be smaller than
   #4 comment 3 assumed.

⚠️ **The golden file says nothing about whether any of the derivative work is
right** — no golden column comes from a gradient. Do not read a bit-identical
golden run as evidence here. The checks that do bite are the three arbiters used
above: a derivative-free scan of the maximum, `profit` not decreasing, and
`|dprofit|` at the returned point.

#### What implementation changed — three things the write-up above got wrong

Recorded because each cost a debugging cycle and two of them contradict what is
written above.

1. **Returning the raw bound on a pinned row is WRONG, because of #31.** The plan
   said "return the endpoint when there is no interior crossing", which is what any
   safeguarded method does. It moved 22 golden rows' profit **down by 1.44**. The
   cause is #31: below `psi_upstream` the profit algebra runs on a negative
   conductance, so profit is **discontinuous** across the feasibility boundary
   rather than merely steep — at `psi_soil = 4`, `vpd = 2`, 5 layers it is −4.696
   just inside and −6.136 at `bound_a`. So #4's warning to "keep the search inside
   the bracket" applies to the *answer* as well as to the search. The fix is to
   return the stepped-inside point, which also beat golden section on profit at both
   worst rows. Consequence worth carrying: a wet-pinned answer is determined to the
   step-in scale (~1e-6 of the bracket width), not to `collar_root_tol`.
2. **`profit` is the wrong instrument for checking this change**, and the write-up
   above proposes it. It is the maximum, so it is flat, and its own floor is the
   nested `ci` root-find's 1e-7 tolerance: two of 288 rows end ~6e-7 *lower* in
   profit while their residual improves by ten orders of magnitude. **Check
   `|dprofit|` at the returned point** — 240/240 improved, none worse. A test
   asserting "profit never decreases" would have encoded the noise floor.
3. **The speedup is real and larger than "probably faster" — 24.5%**, 2.65 against
   3.51 µs/solve, identical across four interleaved pairs at reps=2000. The
   eval-count estimate above (12 vs 18, "not a speedup") was too pessimistic
   because it priced a `dprofit` evaluation without noticing the other half of the
   change: the gradient was split into a wrapper plus `dprofit_at_collar_psi`, so
   the solve seats the supply caches **once** rather than per evaluation — the same
   saving #530 made for the finite-difference path.

⚠️ **And one process note.** R does not track header dependencies, so
`R CMD INSTALL` after editing `inst/include/` silently reuses a stale
`src/RcppR6.o` and the R layer keeps running the OLD model. Two rounds of R-side
diagnostics were measured against the previous solver before this was spotted, and
the numbers looked plausible throughout. **`rm -f src/*.o src/*.so` before
reinstalling**, and sanity-check one value against the C++ suite.

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

**Do this after item 11 — confirmed as a decision, 2026-08-04, not just a
preference.** Attempting it before produces a vignette that finite-differences
trait gradients, which is precisely the thing being argued against. There is no
partial version worth doing either: a calibration built on FD gradients would
have to be thrown away rather than upgraded, because the headline result is the
comparison between the two. **Item 11 (#4) is therefore the gate on this item,
and nothing here should start until AD trait gradients exist.**

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
- **Run the C++ suite under `R CMD check`** — DONE (#12). `tests/cpp.R` compiles
  and runs `tests/cpp` using R's *configured* compiler (`R CMD config CXX20`)
  against the *installed* headers, which is what a `LinkingTo` consumer actually
  gets. `R CMD check` is clean: one NOTE, `'LinkingTo' field is unused: package
  has no 'src' directory`, which is inherent to a header-only package and cannot
  be fixed without adding compiled code we do not want.

  It compiles the two sources **directly** rather than calling `make`, and that
  is deliberate. `R CMD check` scans every Makefile in the tarball and warns
  about GNU extensions, which `tests/cpp/Makefile` uses freely — `?=`,
  `$(shell)`, `$(wildcard)`, `ifeq`. The sanctioned way to silence that is
  `SystemRequirements: GNU make`, and it was declared briefly before being
  reverted: it would be a **false statement about the package**. Installing leaf
  needs no make whatsoever — there is no compiled code — and because
  `SystemRequirements` is package-level metadata, every `LinkingTo: leaf`
  consumer would inherit a declared dependency that is untrue for them. Only the
  developer harness needs make. So the Makefile stays for developers and is kept
  out of the tarball via `.Rbuildignore`, which costs two compiler invocations in
  `tests/cpp.R` and no lies in the metadata.
- **Doxygen for the C++ API** — DONE (#12). Not roxygen, and not a close call:
  roxygen documents R objects and this package has none. When item 6 gives it an
  R interface, roxygen becomes right for *that* surface and the two coexist.

  The obstacle was that every comment here is a plain `//`, which Doxygen does
  not read as documentation, and converting ten headers to `///` would be a
  large diff through exactly the files #15 was rewriting at the time.
  `tools/doxygen_filter.awk` does it at render time instead — comments only,
  sources untouched. It also promotes each file's opening block to `\file`,
  wraps indented tables and examples in `\verbatim` so Markdown cannot reflow
  them into a paragraph, and escapes the characters Doxygen would otherwise read
  as commands. CI renders with warnings promoted to errors, and asserts that
  every non-comment line survives the filter byte-for-byte — because a filter
  bug would otherwise show wrong signatures with nothing to indicate it.

  Publishing is opt-in: set the repo variable `PUBLISH_DOCS=true` and Pages'
  source to "GitHub Actions". Until then every run uploads the site as the
  `api-docs` artifact.
