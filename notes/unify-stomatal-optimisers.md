# Unify the stomatal optimisers: one Sperry, a measured stem curve, gradients for every model

**Branch `refactor/unified-psi-stem-solver`. Status as of 2026-08-22, rewritten from the working plan that lived outside the repo in `~/.claude/plans/`.** It lived there for the first half of this work and was lost sight of for the second half, which is why it is now version-controlled beside the code it describes. `PLAN.md` remains the per-issue status document; this is the plan for one branch.

## Where this started, and the premise it refuted

The original question was whether `optimise_psi_stem_Sperry` / `_TF` / `_ProfitMax` could collapse into one function differing only in cost and benefit.

**Refuted, measured.** `optimise_psi_stem_TF` and `find_root_collar_psi` are *not* two implementations of one model: they optimise different variables over different supply topologies, disagree on 20 of 30 driver rows, and the gap does **not** close as the root resistance → 0 (1.67 MPa at r = 1e-2). The code always said so — "non-**root-based** profit optimisation methods". As that resistance vanishes the collar loses its freedom, `[root_zero_E, root_crit]` collapses, and the collar solve correctly reports `determined` rather than optimising. So the decision variable is a structural axis in its own right, and the entry points cannot collapse across it.

What *could* collapse, and now has, is everything on one side of that axis: the single-layer optimisers share one body.

Also measured: **#119 made the single-layer solvers exact** (0 of 30 rows short against a 20001-point reference). So a first-order-condition-only route rests solely on the gradient requirement — sufficient, but narrower than assumed, because a scan argmax is piecewise-constant in the traits and cannot feed a trait gradient.

## Decisions settled

| question | decision | state |
|---|---|---|
| The prescribed-λ "Sperry" arm | **Deleted.** ProfitMax *is* Sperry, and Sperry's point is that the cost scaling is emergent; a fixed λ contradicts the model rather than varying it. Keeps the arm computing `λ* = \|A\|max / k_span` each step, and frees `lambda_` to mean Cowan-Farquhar's λ alone | ✅ landed |
| Stem curve parameterisation | **(P50, c)** | ✅ landed |
| `psi_crit` | **Derived, not a parameter** — for every model, from one named critical-fraction constant that Sperry's `k_crit` also reads | ✅ landed |
| Model-scoped parameter names | **Prefixed with the model's code name**, initials-plus-year where the paper has no name: `TF24_*`, `CF77_*`, `JS22_*`, `CMax_*` | ✅ landed |
| λ in the gradient enumeration | **Appended one slot**, as `CF77_lambda_` rather than the `par_lambda` the original plan named -- the rename landed first | ✅ landed |
| A parameter absent from the active model | **Refuses, naming the model** and which of the two axes -- supply path or active model -- ruled it out. Distinguishes "structurally not in this objective" (refuse) from "in the objective but inactive at this operating point" (zero + `status`) | ✅ landed |
| ProfitMax's derived λ | **Not an available parameter** for it, and it has no gradient route at all: its normaliser comes from a scan whose argmax is not differentiable. `vcmax_25` and `jmax_25` stay partials at fixed λ | ✅ landed |

## The model space, as it now stands

Seven single-layer entry points plus the collar solve. Five subtract a cost; two multiply.

| code name | paper | cost `C(ψ)` | parameters |
|---|---|---|---|
| `TF24` | Towers & Falster (2024) | `γ(1−f)^β` | `TF24_cost_scale`, `TF24_beta2` |
| `CF77` | Cowan & Farquhar (1977) | `λE` | `CF77_lambda_` (an input, no default) |
| `ProfitMax` | Sperry et al. (2017) | `[k(ψ_s)−k(ψ)]/[k(ψ_s)−k_crit]` | none — λ* is emergent |
| `JS22` | Joshi & Stocker (2022), hydraulic term only | `γ_J(Δψ)²` | `JS22_gamma` |
| `CMax` | Wolf et al. (2016), as Anderegg et al. (2018) parameterised it | `½a(ψ²−ψ_s²)+b(ψ−ψ_s)` | `CMax_a`, `CMax_b` |
| `SOX` | Eller (2018, 2020) | product: `A·g`, `g = (f−0.05)/0.95` | none |
| `JW26` | Jones et al. (2026) | product: `A·g`, `g = 1−ψ/ψ_crit` | none |

The last four were added after the original plan was written, so every count in it was stale. The trait vector is **15** and `gradient_par_names()` is **18** -- twelve traits, the three the new curves brought, and three non-traits (`leaf_specific_conductance_max`, `resistance`, `CF77_lambda_`) -- against the 15 the old plan projected.

⚠️ **`JS22`, `CMax` and `JW26` are the objectives, not the models.** `JS22` omits Joshi's `α·J_max` capacity term and the joint optimisation over `J_max` that goes with it. `JW26` uses a derived `psi_crit` and our saturating supply where the paper supplies its own `Pcrit` and a supply linear in the potential. Deriving `psi_crit` is deliberate — it is what puts `JW26` and `SOX` on the same two anchors, so the comparison between them is structural rather than a fitted match — but it means neither reproduces its paper's numbers.

## What landed

**Coverage and instrumentation.** `psi_stem_optima.tsv` now records **4608 rows** over solver × topology × energy balance × thermal cost × drivers, in two passes: a fresh `Leaf` per row, and one reused `Leaf` with outputs poisoned so stale state is visible. `operating_points.tsv` (576) and `primitives.tsv` (544) are unchanged throughout. `bench_solve` prices every arm.

**One optimiser body.** The six single-layer curves that share a shape — `TF`, `CF77`, `JS22`, `CMax`, `SOX`, `JW26` — were six near-copies, 61–67% line-identical. They are now one templated body over three `if constexpr` tables (`check_cost_parameters`, `profit_psi_stem_for`, `lambda_for`), each written once with every arm explicit and the last asserting. `util::maximise_over_closed_interval` is called from one place where it was called from six. `evaluate_psi_stem` shares the same tables. Net −66 lines, all 4608 golden rows bit-identical.

`ProfitMax` stays outside deliberately: it seeds `|A|max` and the conductance span before searching and writes its own normalised members, so its body is a different shape rather than a differently parameterised one.

**The (P50, c) reparameterisation.** `stem_b`/`root_b` are derived read-only accessors; `psi_crit`/`root_psi_crit` are derived from one named critical-fraction constant. Two checks were deleted rather than maintained: `check_psi_crit_domain` is provably vacuous (a derived P95 is inside P99 for any `c > 0`, since `ln 20 < ln 100`), and the consistency check has nothing left to compare.

It could not be bit-identical, and not for the reason assumed: the defaults were round in `stem_b` and `psi_crit`, not in P50, so today's `psi_crit` was P95 rounded to seven figures and `f(psi_crit)` was `0.049999968739194864` rather than 0.05. Adopting round `P50 = 3.4` shifts `stem_b` by +3.13e-08 and `psi_crit` by −4.65e-08, which propagates over the 576-point grid to a worst **4.39e-05 on assimilation** — a ~1000× amplification of a 5e-08 parameter shift, biologically nothing but ⚠️ **just under this project's own 1e-4 "real difference" threshold rather than far below it.** Approved on that basis and the reference values updated.

**The three chain rules are asserted** (`test-gradient.R:82`): `d/dP50` at fixed `c` equals `(ln 2)^{-1/c}` times the b-gradient — on interior rows only, because a derived `psi_crit` adds a term that is exactly zero at an interior optimum and order 1 at a pinned one; and `d/dc` at fixed `P50` equals `d/dc|_b + (db/dc)·d/db` with `db/dc = −0.36651·b/c²`.

**`dprofit_dpsi_stem<CostCurve>`** exists for the additive curves, verified against a central difference of the objective itself. ⚠️ It refuses the two product curves at compile time: it computes `dA/dψ − dC/dψ`, and a product's derivative is `(dA/dψ)·g + A·g'`.

**A maths vignette**, `vignettes/the-models.Rmd`, covering the objective, the cost table, λ as the common currency, the gates, the degenerate states, and appendices on (P50, c), gradients and numerical practice. Every number in prose is produced by the chunk beside it.

**Two bugs found and fixed, both worth carrying:**

`bench_gradient` segfaulted for three trait-count changes. The timing loop ran to a literal `13` over an 11-element `kBase`/`kNames`, handing a garbage pointer to `printf("%s")`. It hid twice over: the read is only fatal at some address layouts, so it ran clean under `lldb` (which disables ASLR) and crashed when run directly; and a piped stdout is block-buffered, so the report died in the buffer and it *looked* like a crash during initialisation. `-fsanitize=address` named it in one line. **Reach for ASan before reading a stack trace.**

The product objectives reported a shut-down leaf as profitable. Past `psi_crit` the reduction factor goes negative, and the objective is `A·g` with `A = −R_d < 0` at closure — so two negatives gave **+0.0125 at ψ_soil 6.0 and +0.0633 at 7.0**, against TF24's −8.48 and −8.85, growing as the soil dried. Both factors now clamp at zero in the NaN-propagating form (`g < 0 ? 0 : g`, so a NaN is returned rather than silently zeroed).

## Open work

Items 4, 5 and 6 have landed. What is left is the batch gradient route, the comment sweep, and the deferred first-order-condition polish -- which item 5 promoted from optional to load-bearing.

### Item 4 — λ into the gradient, and per-model availability — LANDED

Append `par_lambda` (17 → 18). Extend `.gradient_available_pars(single)` → `(single, model)` at `R/gradient.R:903`, which is already the single home for this rule and is consumed by both entry points and `print.leaf_batch`. Absent slots carry `NA` and are never read, exactly as `resistance` does today.

**Partly landed**: the two positional tail writes in `.gradient_theta_matrix` are now name-based. They were `length(par_names)` and `length(par_names) - 1L`, correct only while `kmax` and `resistance` are the last two columns in that order — so appending λ would have written `kmax` into the new slot and `resistance` into `kmax`'s, with no length change to notice.

### Item 5 — differentiate any additive curve's optimum — LANDED, with one measured limitation

`leaf_gradient()` gains `model`. `"collar"` (the default, bit-identical) is the TF24 cost over the root-collar potential; `"TF24"`, `"CF77"`, `"JS22"` and `"CMax"` differentiate the corresponding **stem** optimum. Checked against differencing the solve on all four: agreement 5e-05 to 2e-03, which bounds the finite difference's error rather than the composite's. `CF77_lambda_` is differentiable, which is what item 4 appended it for.

⚠️ **A stem route reports `status == "pinned"` at a perfectly interior optimum, so `method = "auto"` differences the solve.** Measured: the collar solve root-finds `dprofit == 0` and leaves stationarity at **1.15e-15**; the stem entry points scan a 64-point grid and refine the winning cell, leaving **2.2e-08 to 1.2e-06**. `stationarity_tol`'s 1e-08 default was measured in the empty band between the *collar* route's own two populations and does not transfer. The premise is untested on a stem route rather than failing there — `method = "ift"` gives numbers that agree with differencing the solve.

**So the deferred FOC item is now load-bearing rather than optional.** A first-order-condition polish inside the winning cell would make the stem optima stationary and let `"auto"` use the composite. That is the same multi-start-FOC work the last section defers, and this is the concrete reason to do it. The test asserts today's classification, so the polish will fail it and say so.

⚠️ **What remains for the calibration**: `leaf_gradient_batch()` — the C++ composite in `gradient.hpp` — still hardwires `find_root_collar_psi()`. The R route is done and the batch route is not, and the batch route is the one a fit uses. It must stay bit-identical with R, which is what makes it the larger half.

⚠️ **`SOX` and `JW26` have no route.** They maximise `A * g(psi)`, so the derivative is `(dA/dpsi)*g + A*g'` rather than `dA/dpsi - dC/dpsi`. Refused by name. The plan's warning about the log link putting the residual in log units is still ahead, not behind: the first thing to bite was the *solver*, not the units.

Keep the active-set classification and the envelope theorem untouched: they depend only on stationarity, which is the payoff of the FOC route.

### Item 6 — the (P50, c) converter — LANDED; the downstream copies are not retired

`psi_at_plc()`, `weibull_s50()` and `weibull_p50_c()` are exported, round-tripped through every admissible pair, and `psi_at_plc` is asserted against the model's own derived `psi_crit` and `stem_b`.

⚠️ **`S50` paired with a quantile other than P50 does not determine the curve.** `S50(c) = K*c*L^(1/c)/psi_f` has an analytic minimum at `c = ln L`, so it is U-shaped and any slope above that minimum is matched by **two** curves. Found by the round trip failing with "no root in (0, 20]" — both bracket ends were positive because the root sat in a dip. Each branch is searched from the analytic turning point; the upper one is returned because measured angiosperm stems sit at `c = 1.8-2.6` while the lower root lands at 0.57-0.88, and the other solution is reported rather than discarded.

⚠️ **The downstream copies are NOT retired, so the derivation still exists twice.** Published curves arrive as almost anything except `(P50, c)` — a `P88`, a `P95`, a `P50` plus a slope — and all of it inverts in closed form. For a remaining conductivity fraction `f`:

$$P_f = P_{50}\big(\log_2(1/f)\big)^{1/c} \quad\Longleftrightarrow\quad P_{50} = P_f\big(\log_2(1/f)\big)^{-1/c}$$

with `log₂(1/f)` = 1, 3.0589, 4.3219, 6.6439 for `f` = 0.50, 0.12, 0.05, 0.01. From any two quantiles both parameters follow, `c = ln(L₂/L₁)/ln(P_{f₂}/P_{f₁})`. And the slope at P50, the other commonly published form, is `S50 = 34.657·c/|P50|`.

Ship one exported function taking any two of `{P50, P88, P95, P99, Px at a given x, S50, c}` and returning `(P50, c)`, plus the forward `psi_at_plc()`. Do **not** make S50 internal — it has P50 built into it by construction. Defaults and guards: `c = 2.0` for angiosperm stems, guard `c > 0`, warn below ~1 (usually a measurement artefact rather than a trait).

⚠️ **The derivation currently exists twice.** `leaf_calibration_test` has four hand-written `psi_at_plc` copies (`01-leaf-predict.R:62` and `:154`, `03-analysis-steps.R:221`, `05-exact-gradient.R:99`, `08-unmeasured-traits.R:286`) plus its Jacobian rows and a `stop()` refusing to free `stem_b` and `stem_c` separately *because* `psi_crit` depended on both. Item 3 said to retire them in the same change and that did not happen. Retiring them is part of this item, or the model and its consumer disagree the first time either is edited.

⚠️ **Sign convention.** Published `P50` is negative (−3.95 MPa); this package is positive magnitudes throughout, so `P50` here means `|P50|`.

### Item 7 — the comment sweep

Code, `README.md`, `.claude/CLAUDE.md`. A separate code-free change. `leaf_model.hpp` is now ~4400 lines of which over half are comments.

**The test**: would this help a reader who never saw the previous version? If it only makes sense as a diff against something gone, it goes.

**Remove** — change narrative ("used to", "no longer", "before #X", struck-through text); issue, PR and plan identifiers; decision provenance ("tried and reverted", "closed unmerged"); self-narrative ("it cost me", "this entry said otherwise"); dated claims.

**Keep, rewritten in the present tense as a rule rather than a story** — contracts; hazards as standing rules rather than incident reports; measurements that inform a present decision (keep the number, drop the campaign); the reason for a non-obvious choice where it constrains future edits.

⚠️ **Over-deletion is the risk.** Much of that mass is load-bearing hazard documentation with measurements attached — the two Weibull curves, the accumulator-width trap, the argmax-smoothness constraint. Those survive; only their history goes. Where a comment's only content is history, check first: if "X used to be Y and that broke Z" is the only record of Z, the replacement is "X must be Y because Z", not silence.

## Deferred, with a number attached

**The interior second hump.** A FOC-only route knowingly regresses an arm that currently measures exact. The fix must be **multi-start FOC, not a grid** — several sub-brackets, each an exact stationary point, compared on profit — because a grid argmax cannot feed a gradient.

The number motivating it **cannot be measured until FOC-only exists**: it is a property of the replacement, not of the current code. All arms are exact on all 1920 rows of the current sweep (5 ψ_soil × 4 T × 3 VPD × 2 PPFD × 4 gate combinations, against a 20001-point reference of each arm's own objective) — zero rows short, worst 0.000e+00. The 3.05e-01 figure sometimes quoted was measured with **endpoints only**, i.e. the half-fix before the scan existed. So: implement FOC-only behind a flag, measure its shortfall against the same reference, and only then decide whether multi-start is needed.

Deleting the prescribed-λ arm removed the most *structurally* multimodal cost curve: its marginal cost was `λ·kmax·|f'|`, and `|f'| → 0` at **both** ends for `c > 1`. Keep a synthetic peaked curve in the probe if that stress case is wanted later.

## Live hazards

1. **`n_pars` is a compile-time constant** in fixed-size stack arrays and in test/bench sites that initialise `theta` with a literal list. An aggregate initialiser shorter than `n_pars` is legal C++ and zero-fills, so adding a parameter shifts `kmax` and `resistance` down a slot and drops `resistance` off the end **with no diagnostic**. This bit twice in one session — `test_leaf.cpp`'s batch `theta` (symptom: three observations failed to solve) and `test_golden.cpp`'s `by_solver_threw[4]` indexed by an enum that grew to 8. Both now carry `static_assert`s tying the count to one named constant. Count against `n_pars` when touching either.
2. **The trait vector is bound positionally in four places** — C++ `apply()`, R's `.gradient_setter`, the batch route, and R's derived enumeration. Two fail loudly on a length change; the positional tail writes failed silently until item 4 named them.
3. **`operating_points.tsv` validated plant's 78/78 SCM nodes.** Anything here must leave it bit-identical or move it deliberately with an explanation. The collar path never calls `optimise_psi_stem_*`, so a movement there means a leak into shared code.
4. **A derived `psi_crit` makes the P50 identity conditional** on the optimum being interior. Asserting it unconditionally fails on pinned rows for a correct implementation.
5. **The Choat 2003 `c` values are not trustworthy** (13–14 points each, already 25–30% PLC by 0.8 MPa — a native-embolism artefact). The defensible claim rests on the three Blackman angiosperms measured the same way.
6. **plant is not updated.** Four parameter names moved (`beta2`, `cost_scale_TF24`, `lambda_`, and `lambda_analytical_`, which is deleted here but still bound in plant's `inst/RcppR6_classes.yml:79`). plant tracks this package's master via `Remotes:`, so merging breaks it until they land together.

## Verification

```sh
cd /Users/z2209343/GitHub/plant-family/phylloptim
make -C tests/cpp clean && make -C tests/cpp     # test_leaf, test_golden, test_primitives
make -C tests/cpp bench                          # CI builds these; `all` does not
rm -f src/*.o src/*.so && R CMD check .          # tests/cpp.R against INSTALLED headers
Rscript tools/fingerprint.R
```

`rm -f src/*.o src/*.so` is not optional — R does not track header dependencies, and a stale `.o` reads as a real bug or as nondeterminism between processes.

Anything that moves recorded numbers needs `gradient_golden.tsv` regenerated on macOS/arm64 from an *installed* build (`tools/gradient_golden.R` enforces both) **and** `leaf_behaviour_fingerprint()`'s digest bumped, or `test-fingerprint.R` fails. That golden file's `single-potential` case is the only coverage of the two non-trait parameters.

⚠️ **When finite-differencing anything computed through a root-find, sweep `h` and take the floor of the V.** The noise floor is solver-tolerance/`h`, so relative error goes as 1/`h` for six decades before truncation takes over: measured 1.3e-01 → 1.9e-07 for `h` from 1e-8 to 1e-3, then back to 1.3e-04 at 1e-2. A "small" `h = 1e-6` sits three orders up the noise side and reads as a wrong derivative.

⚠️ **Check that an assertion was actually exercised.** A first-order-condition test whose λ was 100× off-scale put all nine rows on a bound, so the identity was never evaluated and the test passed. Print the count of rows that reached the assertion. λ's scale is set by the leaf — `marginal_cost_water()` runs 9e4–3e5 at the defaults.
