# Next steps

**What this file is.** The reasoning behind each open issue in
[traitecoevo/phylloptim](https://github.com/traitecoevo/phylloptim/issues), plus the
closed decisions that still constrain the code. The issues are the work queue; this
is the *why* behind each.

**What it is no longer.** It used to carry the full retrospective for every
completed item — 3142 lines, most of it the order in which things came to be
understood. That belongs in the merged PR, which squash-merge copies into permanent
history and which the issue links to. Where a closed item survives below, it is
compressed to the decision, the measurement something else points at, and the ⚠️
that would cost someone time.

⚠️ **Two other files carry knowledge this one deliberately does not.**
[`.claude/CLAUDE.md`](.claude/CLAUDE.md) holds the hazards, the golden-file
magnitudes and the build discipline — where the two overlapped the guide was
consistently the *more current*, so the duplicate is gone from here. `NEWS.md` holds
the user-visible history.

## Status, 2026-08-05

### Remaining

| issue | item | what | note |
|---|---|---|---|
| [#3](https://github.com/traitecoevo/phylloptim/issues/3) | 7a | Make λ(state) pluggable | **unblocked** — #2 is done. ⚠️ Measure the dispatch separately: the cost core is fully inlined where the supply path was not, so 7b's "dispatch is free" does **not** transfer |
| [#4](https://github.com/traitecoevo/phylloptim/issues/4) | 11e, 11f | Trait gradients — stage 2, and `stem_c`'s rebuild | Stages 1 and 3 landed ([#42](https://github.com/traitecoevo/phylloptim/pull/42), [#46](https://github.com/traitecoevo/phylloptim/pull/46)). ⚠️ Read 11e before quoting any speedup |
| [#6](https://github.com/traitecoevo/phylloptim/issues/6) | 12 | Real-data calibration, then inversion | **The one with the most downstream** — it specifies what is left of #4, and doing it found #38, #40, #41 and #52. ⚠️ The *synthetic* vignette is done ([#53](https://github.com/traitecoevo/phylloptim/pull/53)) and **AD lost**: read 12a before repeating the comparison |
| [#52](https://github.com/traitecoevo/phylloptim/issues/52) | 12 | `leaf_gradient()` rebuilds a `Leaf` per call, with no way to pass one in | From #53. A third of a per-observation gradient's cost, and why item 12's "reuse one object" advice cannot be followed as written |
| [#7](https://github.com/traitecoevo/phylloptim/issues/7) | 13 | Energy balance, in priority order | Leaf-to-air VPD is the cheap win and is **not wired in**; free convection is not worth it |
| [#28](https://github.com/traitecoevo/phylloptim/issues/28) | 13 | Temperature-dependent outgoing longwave in the Penman-Monteith Rn | from plant #581 / #567 review |
| [#31](https://github.com/traitecoevo/phylloptim/issues/31) | 31 | `profit_psi_stem_TF` returns a plausible number below `psi_upstream` | Found writing the vignette. Profit is **discontinuous by 1.58** at the boundary |
| [#33](https://github.com/traitecoevo/phylloptim/issues/33) | 6d | Take resistances, not root carbon | **blocked on plant #591**. Before #34 |
| [#34](https://github.com/traitecoevo/phylloptim/issues/34) | 6d | Delete plant's `Leaf` bindings | **blocked on plant #591**. The hazard-7 payoff, and the only stage that can break plant |
| [#38](https://github.com/traitecoevo/phylloptim/issues/38) | 38 | `psi_crit` above the curve's domain fails with an error naming neither | Bounds how far `stem_b` can move in one step; see 11f |
| [#40](https://github.com/traitecoevo/phylloptim/issues/40) | 40 | Cross-validate the Medlyn path against `plantecophys` | From #6. Ties to #3, which needs Medlyn first-class |
| [#41](https://github.com/traitecoevo/phylloptim/issues/41) | 41 | `R_d` cannot be set from R: `rd_to_vcmax_ratio_` is unbound | From #6. **74% Rd overstatement** for *Q. ilex* against Sabot's data |
| [#49](https://github.com/traitecoevo/phylloptim/issues/49) | 14 | Is `kmax ~ h^-1` defensible? Koçillari et al. 2021 find no height trend | Split out of item 14, which had it filed nowhere |
| [#50](https://github.com/traitecoevo/phylloptim/issues/50) | 8 | `H2O_CO2_stom_diff_ratio` is 1.67, the g1 literature uses 1.6 | **2.2% offset** in every reported `g1_eff` |
| [#51](https://github.com/traitecoevo/phylloptim/issues/51) | 10c | `kg_to_mol_h2o` and `kg_per_mol_h2o` are not reciprocals | **0.028%**. The bit-identity baseline it was waiting for now exists |
| [#55](https://github.com/traitecoevo/phylloptim/issues/55) | 6d | The temperature cache silently invalidates a trait sweep | #42 documents the trap; nothing prevents it. A bare field write still compiles |
| [#56](https://github.com/traitecoevo/phylloptim/issues/56) | 8, 10c | `gs` is reported only to CO2, and the hydraulic path mixes kg and mol bases | Coupled to #50 and #51 — settle the constants first |
| [#57](https://github.com/traitecoevo/phylloptim/issues/57) | 12 | A failed solve throws, so a fit needs a per-row `tryCatch` | Wants an opt-in per-row status and a condition class, not a silent NA (#39) |
| [#58](https://github.com/traitecoevo/phylloptim/issues/58) | 15 | A consumer cannot tell when the model's behaviour changed | Results moved three times on `0.1.0`; bump the version, and consider a fingerprint |

Not filed as issues, and small enough to stay here: the **`Control` struct** and the
**constants that are really parameters** (item 10b); the **longer, drier SCM run**
(item 1); the **publication** half of item 14.

### Done

All on `master`. One line each; the PR holds the reasoning.

| item | what | detail |
|---|---|---|
| **1** | Validated against plant — the swap was bit-identical, and the 1 ULP was R's decimal parser rather than either model | [#21](https://github.com/traitecoevo/phylloptim/pull/21), issue #13 |
| **2** | Stale state between solves: all four exits fixed, here and in plant | [#15](https://github.com/traitecoevo/phylloptim/pull/15), [#26](https://github.com/traitecoevo/phylloptim/pull/26), plant #577/#585 |
| **3** | plant consumes this package: `feature/consume-leaf-package` compiles and passes | plant #591 (open), issue #9 |
| **4** | The include graph is R-free, and the guarantee is *directional* | [#19](https://github.com/traitecoevo/phylloptim/pull/19), [#22](https://github.com/traitecoevo/phylloptim/pull/22) |
| **5** | XAD stays as it is — forward mode only, so nothing links | below |
| **6**, 6a–6c, 6d stages 0–3 | The R interface: bindings, a usable surface, the single-potential path | [#30](https://github.com/traitecoevo/phylloptim/pull/30), issues #5, #32 |
| **7b**, 7b-i–iii | The supply path is swappable, and the dispatch measured free | [#17](https://github.com/traitecoevo/phylloptim/pull/17), [#18](https://github.com/traitecoevo/phylloptim/pull/18), issue #2 |
| **8** | λ and `g1_eff` reported, and the multi-layer λ identity verified | ⚠️ tail open as #50 |
| **9** | Closed-form fast path, default off, accuracy characterised | below — ⚠️ φ still unmeasured |
| **10a** | Renames: `stem_b`/`stem_c`, `cost_scale_TF24`; `R` and `n` gone | [#15](https://github.com/traitecoevo/phylloptim/pull/15) |
| **10a tail** | One representation for ψ — positive magnitudes throughout | [#27](https://github.com/traitecoevo/phylloptim/pull/27), issues #8, #23, #25 |
| **10b** | Eight dead entities gone, `set_physiology` 14 → 10 args, the leaf purely intensive | [#15](https://github.com/traitecoevo/phylloptim/pull/15) — ⚠️ tail open, and one claim overstated: see #41 |
| **10c** | The hidden hard-coded atmospheric pressure fixed | [#15](https://github.com/traitecoevo/phylloptim/pull/15) — ⚠️ second half open as #51 |
| **11** | `Leaf<T>` — **closed, superseded, not going to be built** | below |
| **11a, 11b** | The collar solve solves its own first-order condition; golden section gone | [#36](https://github.com/traitecoevo/phylloptim/pull/36) |
| **11c, 11d, 11e** | Trait gradients by the implicit function theorem, with an active-set guard | [#42](https://github.com/traitecoevo/phylloptim/pull/42) |
| **11e tail** | `leaf_solve()` 16× faster; the two non-trait parameters differentiable | [#43](https://github.com/traitecoevo/phylloptim/pull/43), [#45](https://github.com/traitecoevo/phylloptim/pull/45), issues #39, #44 |
| **11f** | A gradient in `stem_b` no longer rebuilds the spline — 24.5× | [#46](https://github.com/traitecoevo/phylloptim/pull/46) |
| **14** | The package is named `phylloptim` | [#47](https://github.com/traitecoevo/phylloptim/pull/47), issue #1 |
| **15** | CI on gcc/clang × Linux/macOS with no R, plus `R CMD check` and rendered C++ docs | [#20](https://github.com/traitecoevo/phylloptim/pull/20) |

---

# Remaining work

## 7a. The stomatal / optimality formulations — #3

The package carries three formulations and only one is usable:

| formulation | state |
|---|---|
| TF24 hydraulic gain-risk | production; the whole solve is built on it |
| Sperry et al. (2017) cost | present but **hardwired to `psi_soil_[0]`**, and nothing routes to it |
| Medlyn et al. (2011) USO | present but **bypasses the hydraulic solve altogether** |

So you cannot run the same drivers through two of them and compare, which is the
obvious thing to want. The goal is each as a first-class member, alongside Prentice
et al. (2014) least-cost and Cowan-Farquhar.

**What should be pluggable is λ, not the cost function.** All of these maximise a
profit, so all satisfy `dA/dE = λ` and differ **only** in λ(state). Given λ, each
collapses to

```
ci/ca = g1_eff / (g1_eff + sqrt(D)),   g1_eff = sqrt(3·Γ*·Patm / (1.6·λ))
```

Verified numerically: **0.02% maximum relative error in ci/ca across 767 operating
points**. So six models become six small functions over one tested numerical core,
and the comparison is apples-to-apples by construction rather than by bookkeeping.
This is the companion manuscript's central result (item 14); provenance is
`main.tex` table `tab:lambda` and `notes/tf24_common.R`.

The interface is two methods — λ at a state, and an analytic `dlambda_dpsi` —
because the collar solve needs the derivative, and a finite difference there would
put solver noise into the first-order condition.

⚠️ **The dispatch cost is unmeasured here, and 7b's answer does not carry over.**
`profit_psi_stem_TF`, `hydraulic_cost_TF` and `assim_colimited` have **no
out-of-line symbol in `test_golden` at all** — they are fully inlined, where the
supply path was already out-of-line. A virtual λ therefore *adds* a call rather than
disappearing into one. Measure it the way 7b-iii did, and see its ⚠️ about making
the tag runtime-unknowable.

⚠️ **Prentice least-cost does not fit the scheme as stated.** Its cost is a ratio,
`(aE + bV)/A`, not a subtracted term, so the `A − λE` equivalence needs an extra
step. That equivalence "was verified numerically" — but **that code is not in this
repo**, so reproduce it before relying on it.

⚠️ **There are two vulnerability curves** and this item touches the cost function,
which is where the mistake was made before: using the root parameters for the stem
cost carried λ ∝ ψ^3.02 into a manuscript draft where it should have been ψ^0.64.

**A reference implementation exists, and it is the useful find here.** Sabot et
al.'s `TractLSM` implements Medlyn, Tuzet, WUE_LWP, CGain, CMax, LeastCost, CAP, MES,
SOX, ProfitMax and ProfitMax2 in one place — essentially this item's whole family,
against which our λ(state) formulations can be checked rather than merely derived.
That is a stronger check than #40's `plantecophys` comparison, which reaches only the
Medlyn path, and it makes "apples-to-apples by construction" testable.

## #4 — trait gradients: what is left

Stage 1 landed in [#42](https://github.com/traitecoevo/phylloptim/pull/42) and stage
3 in [#46](https://github.com/traitecoevo/phylloptim/pull/46); the measurements are
under 11c–11f. Two things remain.

**Stage 2 — compose the gradient in C++, one R call per operating point.** This
carries the R layer's remaining speed case, and #39 re-sized it: with `leaf_solve()`
fixed, a gradient's ~10 R calls per parameter are what is left of the boundary cost.
`Leaf::operating_point_values()` is the pattern to copy — one call returning a flat
vector.

⚠️ **Stage 2 is also what an HMC user needs, which makes it more than a speed fix.**
Stan cannot evaluate a likelihood that calls compiled C++ inside this package, so the
calibration study reached for a derivative-free sampler instead. A gradient composed
in C++ and reachable without the R boundary is the enabler for that, and nothing
else on this list is.

**`stem_c`'s rebuild.** 11f removed the vulnerability-spline rebuild for `stem_b` by
homogeneity; `stem_c` has no such identity, because it changes the curve's shape
rather than its scale. The honest route is an analytic

```
∂G/∂c = ∫₀^ψ exp(−(s/b)^c) · (−(s/b)^c · ln(s/b)) ds
```

by quadrature, used to **correct** the spline rather than replace it. ⚠️ Read 11f's
rejected closed form first: substituting the exact integral for the spline
differentiates the wrong model, and was wrong by 3.5e-3. Not started, and it needs
an accuracy argument rather than a speed one.

**`root_b` is the same identity in `roots.hpp`**, where the root curve has two
splines, seven read sites and a per-layer cache. Deferred rather than attempted:
nothing measures or fits the root vulnerability curve, since P50 and P88 are stem
measurements.

⚠️ **And the customer is further off than "waits for a fit that frees the curve"
suggests.** The calibration study pins `stem_b` and `stem_c` analytically from
measured P50 and P88 in *every* fit, and derives `psi_crit` from them — so 11f's
`stem_b` machinery has never been exercised by it, and neither has #38's domain
bound. Freeing even the cost *exponent* per species, a much smaller ask than freeing
the curve, gave no convergence. Treat `stem_c`'s rebuild as unmotivated until
something actually frees a vulnerability curve.

## 12. Demonstrate calibration — then consider inversion — #6

**The claim.** ~2.8 µs a solve *and* exact derivatives rather than finite
differences makes this an unusually good calibration target: gradient-based
optimisers and Hamiltonian samplers want many evaluations and gradients that are not
numerical noise. plant hit the finite-difference problem directly — the analytic
`dprofit_droot_collar_psi` exists because the FD gradient was too noisy for TF24f's
acclimation tracking, and plant #576 turned out to be an AD-versus-FD branch
asymmetry.

⚠️ **The "worst case for finite differencing" half of that claim is retracted.**
Before 11a it was true: a central difference in a photosynthetic trait returned
exactly zero below a relative step of 1e-4, and in a hydraulic trait it returned the
**wrong sign**. After 11a a plain central difference gives ~4 correct digits across
steps from 1e-8 to 1e-2, and it stays correct at the 42 pinned rows where the IFT
composite does not. So the argument is **AD-is-exact-and-cheaper against
FD-is-workable**, not AD-works against FD-does-not.

**The deliverable.** A vignette that fits traits to `A`, `gs` and `psi_leaf` across
soil moisture and light by gradient-based optimisation, and then does the same fit by
finite differences **from the same start point with the same optimiser** — reporting
iterations, wall time, and whether FD converges at all. ⚠️ If AD does not win
convincingly, say so; that is a useful result, and better found in a vignette than
in a paper. Achievable precision of a target is **~1e-9**, set by `psi_stem_to_ci`.

**Write it against FD, and it is not throwaway.** The old rule here — that an
FD-only calibration is worthless because the headline result *is* the AD-versus-FD
comparison — assumed FD does not work, which 11a fixed. So this item **specifies**
what is left of #4 rather than waiting on it.

**What the calibration study has established about this item.** It is the customer,
it is private, and only its package-facing findings are recorded here.

- **No model derivative has been used anywhere in it.** Its fits used finite
  differences only, and in the 40-parameter hierarchical fit **~98% of the work was
  computing them**. So the comparison has not been made *on real data* — 12a made it
  on simulated data, and AD lost there — and the note in that project saying the
  package "does not expose trait derivatives" predates #42 and #46. ⚠️ **Re-pointing
  it at `leaf_gradient()` is the first action for this item**, not writing another
  fit. Expect 12a's result to carry: its species have ~81 rows each, so a
  per-observation gradient pays the same #52 construction and the same R boundary.
- ⚠️ **`optim` under-reports the cost of finite differences, and the error is a factor
  of the dimension.** `counts["function"]` excludes the calls made computing
  numerical gradients. The honest total is function calls plus one gradient's worth
  per gradient evaluation — which is exactly the quantity AD removes, so getting it
  right *is* the comparison. This is 11e's "read this before quoting a speedup"
  arrived at independently, from the other side.
- ⚠️ **Two package gaps have to close before a real-data fit means anything.** `R_d`
  is not settable from R (#41) and there is no leaf energy balance (#7). Both act
  directly on A, so a fit absorbs the mismatch into `cost_scale_TF24` — the one
  parameter such a study is usually about. `a`, `curv_fact_elec_trans` and
  `curv_fact_colim` are settable but easy to leave at defaults, with the same effect.
- ⚠️ **The deliverable cannot be a self-contained vignette over real data.** This
  package ships no empirical data — its only tabular file is the golden-file
  regression fixture, which is model output — and CRAN wants vignettes offline, fast
  and dependent only on CRAN packages. So the vignette either synthesises its data,
  or vendors a small extract whose licence permits it, or the real-data study stays
  outside the package and the vignette demonstrates the machinery on simulated data.
  Decide that before writing it; it changes what the vignette can claim.

⚠️ **When the wall-clock comparison is made, there is one specific way for it to come
out wrong.** `leaf_gradient()` constructs a `Leaf` per call, at **204 µs — 73
solves** (11e). A fit that calls it once per observation therefore measures object
construction rather than the gradient method, and would conclude that exact gradients
lose on wall clock even where they win on solve count. Reuse one object across
observations — that is what `set_traits()` exists for — and report what the exact arm
pays per call that the FD arm does not. The solve-count half needs no such care:
central differences cost 2N solves per gradient, so the ratio is bounded below by the
parameter count.

### 12a. The synthetic vignette — DONE, and AD lost

`vignettes/fitting.Rmd` ([#53](https://github.com/traitecoevo/phylloptim/pull/53))
fits `vcmax_25`, `jmax_25` and `cost_scale_TF24` to `A`, `gc` and `psi_stem` over 72
operating points (12 ψ_soil × 3 PPFD × 2 VPD) of noisy simulated data, recovering
them to −3.5% / +7.5% / −1.2%. Same optimiser, same start, same optimum (201.13):

| | objective evals | honest total | per gradient |
|---|---:|---:|---:|
| Nelder-Mead | 106 | 106 | — |
| L-BFGS-B, numerical | 18 fn + 18 gr | **126** | 6.9 ms |
| L-BFGS-B, **exact** | 19 fn + 19 gr | **19** | **33.5 ms** |

`honest total` adds back the `2 × npar` central differences `optim` does not count.
**The exact gradient cuts the solve count ~5× and is 4.9× slower in wall clock**, and
per the instruction above that is what the vignette says.

⚠️ **This supersedes the construction warning above.** That confound is real and is
only a **third** of the gap, so it changes the conclusion rather than explaining it
away: warmed and interleaved, 32.8 ms for the exact arm of which 11.0 ms is the 72
constructions. The other two thirds is the **R call boundary** — the composite trades
one solve for several much cheaper evaluations, each an R call costing more than the
C++ work it wraps. Fixing construction alone would not make the exact arm faster;
only stage 2 (composition in C++) reaches the rest.

⚠️ **And "reuse one object across observations" cannot be done today**, which is why
that advice did not save the measurement: `leaf_gradient()` constructs internally and
takes no model argument. That is **#52**.

Two things the vignette records because both cost time here:

- **Plain `method = "BFGS"` overflows on its first step.** The objective is scaled by
  `1/sigma²` with `sigma_gc = 0.0015`, so the gradient is ~2e4 in log units; from an
  identity Hessian the first trial point is `theta − g` and `exp()` of that is `Inf`,
  so the objective throws before the line search can shrink it. Bound the parameters.
- **A third of an ordinary drought design is pinned** — 48 interior / 24 pinned of 72
  — and the fit still recovers the traits, because the active-set guard differences
  the solve there. The fallback is a normal path, not an exotic one.

⚠️ **Warm up before timing anything in this layer.** An unwarmed `leaf_model()` reads
~3× its steady-state cost. #52 first claimed 89% off exactly that error, and the
wrong number pointed at the wrong fix.

**Two companion vignettes landed with it.** `cpp-interface.Rmd` documents the
linkable half and states that trait gradients are **R-only** — a C++ consumer gets
`set_traits`, `perturb_stem_b`, `dprofit_droot_collar_psi(psi, &feasible)` and
`evaluate_root_collar_psi`, but must compose the theorem and reimplement the
active-set guard itself. `phylloptim.Rmd` gains "Choosing which gradients to
compute": `pars` selects, omitting it computes all sixteen, and all-16 costs ~4.6×
one. ⚠️ The four *outputs* are not selectable — `pars` chooses rows, never columns.

**Then the bigger question: inversion.** `plantecophys::fitaci` inverts A-Ci curves
for Vcmax, Jmax and Rd; this package inverts nothing. Inverting *hydraulic* traits —
P50-type vulnerability parameters, root conductances — from gas exchange plus water
potential would be a genuinely new capability rather than a reimplementation, and no
R package does it. ⚠️ Scope it properly: it is a project, not a task, and it needs
someone to decide what data it is meant to consume.

## What a calibration can and cannot identify

Three structural results from the calibration study. They are properties of **the
model's own parameterisation**, not of any dataset, which is why they belong here:
anything that fits this model will meet them, and two of them suggest what the
package should offer.

- **`kmax` and the soil-to-collar series resistance are not separately identifiable
  from leaf water potential.** Only their total is. On the single-potential path that
  means the split between plant and soil share is whatever the prior says, and it
  should never be reported as an estimate. ⚠️ This bears on #33: if the leaf takes
  resistances rather than root carbon, the interface should make the identified
  quantity the natural thing to supply.
- **`cost_scale_TF24` and `beta2` are confounded by the algebra of the cost itself.**
  The TF24 cost is `cost_scale · (1 − k)^beta2`, so at fixed cost
  `d(log cost_scale)/d(log beta2) = −beta2 · <log(1 − k)>`. The ridge is therefore
  predictable rather than empirical, and **the identified coordinate is the cost at a
  reference loss of conductivity**, not either parameter alone. Worth documenting,
  and arguably worth offering as the parameterisation — these are two of the four
  parameters whose gradients are 100% argmax-mediated (11c), so they are also the two
  a gradient-based fit will struggle with.
- ⚠️ **Reparameterising cannot move a point estimate.** A MAP is invariant to it. It
  buys interpretability and sampling geometry, and nothing else — recorded because it
  was tried as a fix for a misfit and could not have worked.

**And the default `cost_scale_TF24` of 7.5 is not a general value** — it was
calibrated for *Eucalyptus saligna* inside plant. It should be documented as a
species-specific default, since it is the parameter a calibration is most likely to
be estimating.

## 13. Energy balance: the full cut — #7, #28

The Penman-Monteith path (`use_energy_balance_`, default off) is a deliberate
minimal core. What it is missing is not equally valuable, and one piece is
expensive. In priority order:

1. **Leaf-to-air VPD. Do this one.** `saturation_vapour_pressure()` and its slope
   exist and are tested but **are not wired in**, so stomatal conductance still
   divides by the prescribed air VPD however far Tleaf has risen above Tair —
   understating the driving gradient exactly when energy balance matters. Cheap, and
   the pieces are already there.
2. **Temperature-dependent outgoing longwave — this is #28.** Replace the fixed
   `longwave_net_offset = -40 W m-2` with `eps·sigma·Tleaf^4` plus an atmospheric
   term. ⚠️ That makes `Rn` depend on Tleaf, so it needs either plantecophys's
   radiation-conductance linearisation or an iteration. From plant #581 / #567.
3. **Free convection — the expensive one, and probably not worth it.** A Grashof
   term makes boundary-layer conductance depend on `|Tleaf − Tair|`, which makes the
   balance implicit and forces an inner iteration **inside the collar solve**, which
   already runs many inner evaluations. At ~3 µs per solve in a model that calls it
   millions of times that is a poor trade. If leaf-temperature accuracy under low
   wind ever becomes the priority, prefer a one-step correction to a converged inner
   solve.

## 6d. The R interface: what is left — #33, #34

Stages 0 through 3 are done. Both remaining stages are **blocked on plant #591**,
and **2b comes before 4**.

### Stage 2b — take resistances, not root carbon (#33)

**The solve never touches root carbon.** `uptake_impl` and `duptake_dpsi` read
exactly `network_.r_R_H_min` and `network_.r_R_V_sum`, plus `grav_head_z_` and
`max_soil_layer`. `c_r_V`, `c_r_H` and `r_R_V` are carried only because plant
exposes them through RcppR6. So taking carbon makes the leaf own `beta_R_H`,
`beta_R_V`, `dz_` and the 1/3 : 2/3 horizontal/vertical split — none of which is a
leaf-level quantity. `leaf_specific_conductance_max` is the precedent: plant computes
`kmax = K_s·θ/(h·η_c)` and passes it in already reduced.

Two objections, both answerable:

- *It leaves the tested surface.* Keep `root_network_from_carbon` as a documented
  helper, and have the golden grid call helper-then-`set_physiology`, so the tested
  path is still exercised end to end.
- *Cost.* Rebuilding five fresh vectors per call measured **+0.074 µs (0.061 →
  0.135), about +2% of a whole solve**; `std::move` did not help; an in-place
  `RootNetwork& out` overload brought it back to 0.056–0.061. ⚠️ Re-measure
  interleaved, per hazard 5, rather than assuming it still holds.

The R-user argument is the decisive one: a bare-leaf caller has no root carbon
profile, no layer thickness and no opinion about `beta_R_V`, so `set_drivers()`
currently has to invent a default and document it. That is the smell.

### Stage 4 — delete plant's `Leaf` bindings (#34)

Hazard 7 exists only because plant's YAML names `Leaf`'s fields; it **dissolves**
once plant stops doing that. This is the only stage that can break plant, so it is a
coupled PR pair — plant tracks this package's `master` through `Remotes:`.

## 31. Profit below `psi_upstream` is a plausible number on a negative conductance

`profit_psi_stem_TF` evaluated below the upstream potential runs the algebra on a
negative conductance and returns a number that looks fine. Profit is
**discontinuous by 1.58 units** across that boundary. It is unreachable from plant's
solve, which is why it survived — and the first thing an R user does is plot the
profit function. Found writing the stage 2 vignette. 11a's "do not return the raw
bracket bound" finding is the same discontinuity seen from the solver's side.

## 38. `psi_crit` above the curve's domain fails without naming either

The vulnerability curve is built over `[0, b·log(100)^(1/c)]` — 6.05 MPa at the
defaults, against a `psi_crit` of 5.87. Set `psi_crit` beyond that, or shrink
`stem_b` by more than ~3%, and the solve dies with *"Extrapolation disabled and
evaluation point outside of interpolated domain"*, which names neither `psi_crit`
nor the curve. Nothing validates the relationship. This is also what bounds how far
either gradient route can move `stem_b` in one step (11f).

## 40. Cross-validate the Medlyn path against `plantecophys`

The Medlyn path bypasses the hydraulic solve, so it is the one formulation here with
an established external implementation to check against — and #3 needs it
first-class. `plantecophys` also solves the coupled system analytically rather than
by root-finding, which is worth considering on its own merits. Found doing #6.

## 41. `R_d` cannot be set from R

`rd_to_vcmax_ratio_` became a settable C++ member in #15, but it is **absent from
`inst/RcppR6_classes.yml`** and is not a constructor argument, so from R a fixed
0.015 is imposed. `R_d_` is bound read/write, but every `set_physiology()` recomputes
it. Sabot's data imply ratios of **0.0046–0.0302** — a factor of 6.5 — and for
*Q. ilex* the model uses 0.914 against the data's 0.525, a **74% Rd overstatement**
and a systematic ~4% species-varying bias in A. ⚠️ Item 10b's claim that this
parameter "is now settable" is true of C++ and false from R.

## 14. Naming, home, publication

The name is settled — `phylloptim`, in
[#47](https://github.com/traitecoevo/phylloptim/pull/47), closing issue #1 — and
`DESCRIPTION` already points at the right repository. What is left is the paper.

**The software paper writes itself from `COMPARISON.md`**: every R leaf gas-exchange
package assumes a stomatal conductance model, and this one derives stomatal
behaviour from hydraulics.

**The more interesting one is already in progress and is not primarily about this
package.** `Falster-stomatal_analytical_analysis`, in `atelier/2-research/active/`:
*"The marginal cost of water as a common currency for stomatal optimality models:
size dependence, testable contrasts, and a diagnosis"*, targeting New Phytologist or
PC&E. The argument is 7a's — six models share `dA/dE = λ` and differ only in
λ(state), so USO is the generic solution of the family rather than a model, which
relocates the empirical question onto the *shape* of λ.

⚠️ **That manuscript is this package's first customer, not a downstream user** — its
blockers are ours. It has no code artefact beyond TF24-against-plant scripts, and
items 7a, 8 and 9 would give it one. Its highest-priority open science question is
now filed as **#49**: whether `kmax ~ h^-1` is defensible at all, given Koçillari et
al. 2021 find no height trend in leaf-area-specific conductance across 103 plants.

---

# Decisions that still constrain the code

Closed. They are here because something in the tree points at them by number, or
because the measurement is load-bearing and lives nowhere else.

## 1. Validated against plant

The swap was **bit-identical**: plant's full suite 0 fail / 0 error on both builds,
`test-leaf.r` **218 expectations** unchanged, and an SCM regression bit-identical
across **78 of 78** recorded numeric nodes for one- and two-species runs — including
the ODE step sequence, which is the check that mattered, since a perturbation below
tolerance could still have flipped a refinement decision.

⚠️ **The 1 ULP was in the measuring instrument, and the lesson generalises.** 585 of
2592 values differed at 1–2 ULP, decomposing as **345 from R's decimal parser** and
240 from the shutdown NA sentinel. `as.numeric`, `scan` and `read.delim` share a
string-to-double conversion that is **not correctly rounded**, returning a value one
ULP off for about **18% of inputs** — `"26.550866314209998"` parses to
`0x1.a8d0593240001p+4` where the nearest double is `0x1.a8d059324p+4`. The golden
file is written by C++ at `%.17g` and read into R; plant's values are computed
in-process and never touch a string, so the parser perturbed one side only. Reading
through `tests/validate/tsv_to_hex.c` gives bit-identity across all **2352** finite
values, and `compare_primitives.R` confirms it independently on **329**.

⚠️ **Consequence for any new R-side expected value: use C99 hex floats.** Decimal
`%.17g` would fail for ~18% of them against a model that is exactly right.

⚠️ **A flag you set may not reach the compiler.** `R_MAKEVARS_USER` silently did not
apply `-ffp-contract=off` to plant's build during this work. Check the flag appears
on the compile line before concluding anything from it.

**Still open, and small:** a longer, drier SCM run. It buys **state-space
coverage** — reaching shutdown and recovery paths the short run never enters — not
drift protection, which the bit-identity already gives. Belongs with item 9.

## 2. Stale state between solves

All four exits fixed here and in plant (#577/#585). Hazard 8 in the guide has the
four defects and the ⚠️ that a bit-identical golden run is no evidence about them.
What is only here:

**The shutdown semantics, agreed on plant #578.** Zero for `transpiration_`,
`stom_cond_CO2_` and `soil_consumption_`; `ci_` = `gamma_ * umol_per_mol_to_Pa`; and
⚠️ **`-R_d_`, not zero, for `assim_colimited_`** — a shut leaf still respires, and
this is the one an implementer gets wrong.

**How prevalent it was.** Measured over a plant run: **33.78% of cohort-time records
carried at least one stale layer, and 25.73% carried three of five.** Reverting the
one word reports **1.31e-04 and 4.75e-05** kg H₂O m⁻² s⁻¹ of phantom uptake in
layers 2 and 3.

⚠️ Several comments still describe this leak in the present tense
(`tests/validate/compare_with_plant.R`, `tests/cpp/test_golden.cpp`). The 48 × 5
blast radius they cite is right; the tense is not.

## 3. plant consumes this package

`feature/consume-leaf-package` (plant #591) compiles and passes. Three things the
survey got wrong about the risk, none of which was predicted:

⚠️ **Version pinning needs two mechanisms, not one.** plant needed odelia `master`,
not `>= 0.2.0`: `Patch::ode_rates` is non-const while `r_ode_rates` takes a
`const&`. Fixed as odelia#48 → 0.2.1. And this package sat at 0.0.1 across four
merged PRs, so a `LinkingTo (>= 0.0.1)` floor was meaningless. **Pin both ways — a
`LinkingTo` floor and `Remotes: ...@sha` — and bump the upstream version first, or
the floor says nothing.**

⚠️ **`TF24_Environment` sets the `atm_kpa` driver to 100.5, not 101.3.** Deriving
ppm→Pa from it moves offspring production by **+2.4%**, against **+0.10%** for the
entire rest of the swap. **This package's golden grid evaluates at 101.3 and is blind
to it by construction** — a first-class limitation of the safety net.

⚠️ **The `area_leaf` cancellation is right for plant's C++ and a trap for its
tests.** At `area_leaf = 0.05` a test's root system is **20× too weak**; the
critical-demand collar moved **−0.685 → −2.57 MPa** while the zero-uptake collar,
being scale-invariant, stayed put — so **only one of two adjacent guards fired**,
which is how a scaling error hides.

**Method:** build `origin/develop` in a worktree against the *same* installed
dependencies as a control. Without it every failure is ambiguous.

⚠️ **"Checked, nothing upstream" has a shelf life of days.** Re-check with
`git log <base>...origin/develop -- <the files we forked>`.

**Deferred, and user-facing:** renaming plant's `pars.b` / `pars.c` /
`pars.g1_TF24` propagates to hyperpar columns and `inst/scenarios/*.csv`, so plant
still carries the `b`/`c` ambiguity whose cost hazard 1 records.

**One out-of-family consumer exists:**
`Falster-stomatal_analytical_analysis/notes/tf24_closed_form_bench.cpp` includes
`<plant/leaf_model.h>` and reads `l.b`, `l.c`, `l.beta2`. Build it when verifying a
rename.

## 4. The include graph is R-free

odelia's `ode_util.hpp` was the last R touchpoint; fixed upstream as odelia#44, so
`tests/cpp/shim/` is gone and the suite compiles against the real headers. The choice
worth recording: `util.hpp` **throws** rather than taking an `ODELIA_NO_R` guard,
because a guard leaves R as the default, and the wrong default is what produces this
class of bug.

## 5. Where XAD comes from — leave it

Only forward mode (`xad::fwd<double>`) is used, so there is no tape and nothing to
link. Confirmed by `nm`: plant's `leaf_model.o` has **zero** `xad::Tape` symbols.

⚠️ **Reverse mode genuinely does not link header-only, and that is why the decision
holds.** `xad::adj<double>` leaves `xad::Tape<double,1>::active_tape_` undefined; the
definitions live in `odelia/src/Tape.cpp`, the only compiled copy of the tape
runtime. plant reaches it through `PKG_LIBS`; this package's CMake build and any
Python consumer cannot. So reverse mode would cost the no-R build.

**And it would buy nothing here.** `dprofit` is **7.1% of a solve**, so a gradient's
cost is dominated by re-solving and forward mode's O(N) never bites. Only N=1 tapes
are instantiated, so `adj<double,N>` would not link even with the library.

## 6a. "Header-only" — two layers, one-way

`inst/include/` is the model in plain C++ and never includes from `src/`; `src/` and
`R/` are the glue, may use Rcpp, and include downward only. The guide states the rule
and what enforces it — `cpp-tests.yml` building on runners with **no R**.

The rejected alternative is worth keeping, because someone will ask again: a separate
`leafr` package. Two repositories for one model, and the target audience would have
to `install.packages` twice.

## 6c. RcppR6 versus hand-written bindings — both, in layers

Generate the ~**90** low-level accessors (60 field getters/setters plus 30 method
forwards), because hand-writing them is 90 chances to transpose an argument, and that
error class produces plausible numbers with a green C++ suite. Hand-write the surface
above them.

⚠️ **RcppR6 is not on CRAN, and that is defused by not declaring it at all.** The
generated files are committed; RcppR6 appears in neither `Suggests` nor `Remotes`; it
is installed by exactly one CI job, which regenerates and diffs — because a committed
generated file goes stale silently. If it ever has to go, odelia's hand-written
approach is the migration target, and that is cheap because the public API already
lives in the hand-written layer.

## 6d stages 1–3

⚠️ **The C++ suite is the regression baseline and is blind to the R layer**, so a
mistranslated binding argument gives a green suite and plausible R numbers.
`tests/testthat/test-golden.R` ties four golden points back bit-exactly; that is why
it exists.

`leaf_traits()` / `leaf_control()` partition the constructor's **19** arguments, four
of which are tolerances, and a test asserts the partition is exact.

⚠️ **`psi_soil >= 0` is surfaced, not smoothed over.** A friendly wrapper is exactly
where someone would `abs()` it, and that check is all that stands between a pre-#25
script and a plausible wrong number.

⚠️ **The leaf-temperature clamp is deliberately not in `leaf_control()`**, though 6b
and 10b both listed it. It is a finiteness guard on the Arrhenius block, not a
tunable: making it settable lets a caller produce non-finite photosynthetic
parameters and get NaNs back with no indication why.

**Stage 3's footgun was designed out rather than documented.** There is no settable
supply tag at any level; `set_supply_multilayer()` / `set_supply_single()` each
reconfigure completely and clear the solved state; `supply_kind`,
`single_resistance_` and `single_gravity_head_` are bound **read-only**, and a test
asserts the assignments fail. `setup_clean_leaf()` clears **both** paths' soil state,
and `set_drivers()` **refuses** `soil_depth` / `root_carbon_per_leaf_area` on the
single path rather than ignoring them.

## 7b. The supply path — and the multi-layer λ identity

Soil and root transport was **257 of 1117 lines of member-function body (23%)** and
13 state members, entering the solve only as `E_up = f(P_collar)` — which is what
made it separable. Hazard 5 carries the inlining evidence and the dispatch
percentages.

⚠️ **The multi-layer λ identity is the theory this interface was for, and it lives
only here.**

```
λ_multi = λ_TF24 · [1 + kmax·f(ψ_r)/S],    S = dE_up/dψ_r = duptake_dpsi
```

Verified against a central difference: ratio **1.0000**
(`notes/tf24_multilayer_lambda.R`). The bracket is ≥ 1, so **the single-layer λ
always understates the true marginal cost, by a factor of 2 to 12**, and the height
scaling of `g1_eff` flattens from **h^-0.30 to h^-0.15** once it is included. This was
"the one genuine gap, and the blocker for any plant implementation".

⚠️ **Copyability is a hard constraint on any redesign here.**
`make_strategy_ptr(TF24_Strategy s)` takes the strategy **by value** and holds
`Leaf leaf;` as a member (`tf24_strategy.h:375`), so templating `Leaf` on its supply
breaks the `plant::Leaf` alias and ~**12,000 lines** of generated RcppR6.

## 7b-ii. The traps the interface did not account for

1. ⚠️ **`soil_consumption_` is a buffer plant assigns back into.**
   `tf24_strategy.cpp:505-508` writes crown-integrated values into
   `leaf.soil_consumption_[a]`, `leaf.E_up_` and `leaf.profit_`, read back at `:46`.
   So the buffers stay on `Leaf` and are handed over by reference. ⚠️ **The units
   differ across that boundary**: `E_up_` is kg H₂O m⁻² s⁻¹ while
   `soil_consumption_[i]` is mol, converted downstream in plant.
2. **`begin_solve()` is a measured hot-path optimisation**, not structure: it
   collapses ~2 spline evaluations per layer to ~1 by building the soil-side caches
   in one pass.
3. Pointer identity (`&psi_soil == &psi_soil_inverted_`) was the likeliest place to
   lose bit-identity; `roots.hpp` no longer depends on it to be fast.
4. ⚠️ **`duptake_dpsi` returning NaN is a contract, not a failure** — hazard 6 has
   this in full, including the five-method interface a third supply path must
   implement.

## 7b-iii. Staging — and the dispatch measurement

All four stages landed. Stage 1 cost **1.7%** (3.53 vs 3.47 µs, interleaved ×3 at
reps=2000), cause **not attributed**: `uptake_impl`'s argument list went 2 → 6 and it
is called ~10³ times per solve, which is the leading candidate only — the aliasing
hypothesis was **tested, refuted and reverted**.

⚠️ **A zero-carbon layer gets `r_R_H_min = 0`, i.e. infinite conductance where there
are no roots.** Backwards, preserved deliberately, and it looks unreachable from
plant today. Do not "fix" it without checking that.

### Stage 2 — which dispatch mechanism, measured

⚠️ **This block is what hazard 6 points at, and the recipe exists nowhere else.**

Re-measured against the post-stage-1 structure with `SinglePotential` as a genuine
second alternative, all three arms executing the same `MultiLayerRoots` code
(identical bench checksum, so only the *reach* varied), 10 interleaved rounds at
reps=2000:

| option | measured | vs direct | verdict |
|---|---|---|---|
| direct call | 3.507 µs | — | baseline, not an option once there are two paths |
| **enum tag + `switch`** | **3.478 µs** | **−0.8%** | **free.** At or below the direct call in all 10 rounds |
| `std::variant` + `std::visit` | 3.542 µs | +1.0% | real, but a third of the +2.6% the older table claimed |
| template `Leaf<Supply>` | — | — | rejected: breaks the `plant::Leaf` alias and ~12,000 lines of RcppR6 |
| `std::function` | — | — | +0.6%, and a poor fit: this is a stateful four-method interface, not a callback |
| virtual base | — | — | +1.1%, plus clone boilerplate and a heap allocation per `Leaf` copy |

*Why* it is free is a durable property rather than a lucky number: `uptake_impl` is
out-of-line and stays that way (hazard 5), the tag never changes within a solve so
the branch predicts perfectly, and a predictable conditional in front of an
out-of-line call disappears into it. `std::variant` cannot match that because
`std::visit` emits **two out-of-line `__dispatch` thunks** — an indirect jump, which
predicts worse than a direct branch. Visible in the binary:
`nm -C bench | grep __dispatch`.

⚠️ **Reproducing this: the tag must be runtime-unknowable.** Set it from a constant
and the branch constant-folds, every arm measures zero, and it looks like a free
lunch that is not one. The measurement above set it from `argc`. Verify before
believing any arm: `nm -C bench | grep -c SinglePotential` should be 0 for the direct
arm and non-zero for the others — if the dispatch folded, the unused path's symbol
disappears entirely.

The surviving argument for `std::variant` is that it makes an invalid state
unrepresentable. Worth 1.0% only if that risk is real; revisit at three or four
supply paths.

## 8. Report λ and g1_eff as first-class outputs — #50

Done: λ and `g1_eff` are reported, and `g1_eff = sqrt(3·Γ*·Patm/(1.6·λ))` is
comparable to `plantecophys::fitBB` and the values tabulated by Lin et al. (2015).

⚠️ **Open as #50: the unit trap.** This package uses **1.67** for
`H2O_CO2_stom_diff_ratio` where Medlyn (2011) and the g1 literature use **1.6** — a
**2.2% offset** in every `g1_eff` reported for comparison against them.

## 9. The closed-form fast path

Measured against a real `plant::Leaf` (`notes/tf24_closed_form_bench.cpp`):

| route | µs | speedup |
|---|---:|---:|
| exact `optimise_psi_stem_TF` | 2.611 | 1× |
| power law + **one** Newton step | 0.241 | **10.8×** |
| explicit β₂ = 1/c form | 0.056 | **47×** |

**0.051 of that last 0.056 µs is `set_physiology`** — the solve itself has
essentially vanished. Header-only conversion and LTO both measured **zero**.

⚠️ **One Newton step is deliberate. k = 2 is worse in the tail — do not "improve"
it.**

⚠️ **The guard tests an output** (`ci/ca > 0.5`), so the fast path is applied post
hoc with a fallback; it cannot be branched on up front.

⚠️ **Report the realised speedup, not the kernel's.** `1/[φ + (1−φ)/10.8]` — at
φ = 0.2 that is **3.6×**, not 10.8×. Measuring φ is still open.

⚠️ Smoothness: roughness in dA/dh is **0.0015 closed form against 0.0011 exact**. The
old justification for that mattering — golden section versus Brent — is superseded by
11a, but the constraint itself stands.

Being analytically differentiable, the closed form would give an exact demographic
gradient and make TF24f's acclimation apparatus redundant.

## 10a. The naming hazards

`b`/`c` → `stem_b`/`stem_c`, and `g1_TF24` → `cost_scale_TF24`; hazard 1 has the cost
of leaving them unmarked. `mass_root_prop` → `root_carbon_per_layer`, because it is
not a proportion of mass. plant still uses the old names, which is item 3's deferred
rename.

⚠️ **`R` and `n` are gone from the public namespace for a reason.**
`tf24_closed_form_bench.cpp` declares `const double n = l.c*l.beta2 - 1.0` and a
Newton residual `R` inside a scope where `plant::R` and `plant::n` were visible. It
compiles **only because the locals shadow them**.

The strong-types attempt (#23) is closed unmerged; hazard 2 has that history. The one
sentence worth carrying forward: **reach for a type to enforce an invariant you
actually need, not to police one you could remove.**

## 10b. Shrink the input set

Done in #15: five dead entities removed and `set_physiology` taken 14 → 11 arguments,
three further dead constants deleted, and **thirteen** temperature-response
parameters made settable. ⚠️ One of those, `rd_to_vcmax_ratio_`, is settable in C++
but **not bound to R** — that is #41, and this item overstated it.

`area_leaf` is out of the supply contract, not merely out of `set_physiology`. Why
that was safe is a **homogeneity property** rather than a coincidence: every term in
`r_R` is inversely linear in root carbon, so `r_R` is homogeneous of degree −1 in the
carbon vector, and scaling carbon by `1/A` scales `r_R` by exactly `A`. Consequence:
**`SinglePotential::resistance_` is per unit leaf area.** The residual 2 ULP is
reassociation.

**Still open, and here rather than in an issue.**

A **`Control` struct** for `GSS_tol_abs`, `ci_abs_tol`, `ci_niter`,
`vulnerability_curve_ncontrol` and `integration_tol_` — deferred because it is
surgery on a 19-argument constructor, and it overlaps item 6. ⚠️ The
leaf-temperature clamp is **not** on that list; see 6d.

**The constants that are really parameters**, in descending order of case:

- ⚠️ **The six Arrhenius shape parameters** are the strongest case: they are
  `plantecophys`'s `EaV`, `EdVC`, `delsC`, `EaJ`, `EdVJ` and `delsJ`, **whose
  defaults plantecophys changed at v1.4 after a literature review** — and thermal
  acclimation (TF24t) modifies them.
- ⚠️ **`a`, `curv_fact_elec_trans` and `curv_fact_colim` are already settable, and
  that is the problem** — they sit at *E. saligna* defaults, they act directly on A,
  and a calibration that leaves them there absorbs the mismatch into whichever
  parameter it is fitting. Being settable is not the same as being plumbed through
  the surface a fit uses; `leaf_traits()` takes them, so this is a documentation and
  worked-example gap rather than a missing feature. `R_d` is the same family and is
  worse — see #41, where it is not settable from R at all.
- The Penman-Monteith block: `longwave_net_offset = -40` is an explicit placeholder
  (plant #581, and item 13); `sw_abs_per_par = 2.0` folds in an absorptance that
  `tealeaves` exposes as `abs_s`/`abs_l` and `plantecophys` as `LeafAbs`;
  `latent_heat_vap` and `vol_heat_cap_air` are fixed at 25 °C.
- The nine Bernacchi constants are a weaker case, but `bigleaf` exposes them.

## 10c. The two latent inconsistencies — #51

The first is fixed: `umol_per_mol_to_Pa = 0.1013` hard-coded atmospheric pressure
(`0.1013 = 1e-6 × 101300 Pa`) while `atm_kpa_` was live on the conductance side, so
the model was internally inconsistent away from sea level — the same trap
`plantecophys` warns about in bold, but undocumented here. Now derived from
`atm_kpa_`.

The second is **open as #51**: `kg_to_mol_h2o` and `kg_per_mol_h2o` are not
reciprocals — `1/55.4939 = 0.0180200` against `0.018015`, a **0.028%** discrepancy,
deliberate to preserve results.

## 11. `Leaf<T>` — closed, superseded, not going to be built

Both payoffs are spent. The `namespace detail` AD replicas were deleted by 11b using
scalar-generic *member* templates rather than a class template, and "templating is
what turns AD into calibration" was falsified by 11e plus a calibrating #6.

⚠️ **Do not reopen this to get exact trait derivatives; that is done.** If it returns
it should return from the plant side, under a title describing an **outer** AD pass
over plant — which is what plant #537's cut-point rule is about.

**And even that now looks avoidable.** The composite works for a **driver** as well
as a trait: differentiating with respect to `leaf_specific_conductance_max` — which
is how height reaches the leaf — gives ratios of **0.99994–0.99996** against a
resolved slope, with the mixed partial stable to **7 significant figures across five
decades of step**, and the indirect term is 100% of `dA/dk`. ⚠️ Measured at one
interior point on the leaf's share only, so it sizes the idea rather than settling it
— and the active-set hazard transfers, which matters *more* in height than in a
trait, because the SCM transports along height.

## 11a. The collar solve solves its own first-order condition

Golden section is gone: `dprofit == 0` by safeguarded TOMS748. Hazard 3 carries the
smoothness table and the residual improvement. What is only here:

⚠️ **The arbitration, which hazard 3 explicitly delegates to this item.** A wrong
gradient is plausible, so the root-find was checked three ways: a **20001-point
derivative-free scan** of profit agrees with it to **3–7e-10 MPa** and on slope to 5
significant figures; profit never decreases; and `|dprofit|` goes **6.2e-04 →
~1e-14**. Under golden section the collar came back at `dprofit = −6.219675e-04` with
slope −8.9566 MPa⁻¹, i.e. **6.94e-05 MPa** off the true stationary point.

**Why that mattered: the argmax came back smooth, plausible and sign-inverted.** For
traits in the hydraulic path the wandering offset dominated the real response —
`root_b` gave **−2.639e-03 where the truth is +2.587e-04**, with second differences
of 4.7e-08, so it looked *clean*. `beta_R_H` likewise: **+1.443e-05 against
−3.510e-05**. For traits that only enter the objective the difference was exactly
**zero** below a relative step of ~1e-4, and `vcmax_25` came back **8.21e-03 against
a true 1.72e-02** — 52% low. ⚠️ **A staircase announces itself; a smooth wrong sign
does not**, and no step-size tuning or AD-through-the-iterations fixes it, because
the error is in the solved argmax rather than in how it is differentiated. plant #406
predicted the class in the abstract; this is a measured instance.

⚠️ **Do not return the raw bracket bound on a pinned row.** Profit is
**discontinuous** across the feasibility boundary (#31): −4.696 just inside against
−6.136 at `bound_a`, so the obvious safeguarded-solver move moved **22 golden rows'
profit down by 1.44**. Return the stepped-inside point. A wet-pinned answer is
determined to ~1e-6 of bracket width, not to `collar_root_tol`.

⚠️ **The `0.0` shutdown sentinel fires at the wet bracket endpoint**, which is where a
bracketing solver probes first, so it reads as a root (profit −1.897 there against
2.516 at the optimum). The region is narrow — **≤3.46e-07 MPa** into the bracket,
median 1.22e-08, ≤6.2e-07 of bracket width — so the fix is cheap, but it is a
prerequisite rather than defence in depth. `test_leaf.cpp` states the mechanism.

**42 of 240 feasible rows are pinned** (24 wet at `root_zero_E`, 18 dry at
`min(root_crit, supply_psi_crit())`), all at `psi_soil` 3–4, so the collar solve is a
**constrained** one. ∂²profit/∂ψ² ranges **−1.56 to −61.4**. ⚠️ The interior→pinned
transition is a genuine **kink** in trait space — a property of the constrained
problem, present under golden section too, not a regression.

**Two things found while reading XAD and plant, both still live:**

⚠️ **XAD's vector forward mode has a trap.** `xad::fwd<T,N>` exists, but
`ExprTraits<FReal<Scalar,N>>::vector_size` is hard-coded to 1 while `FRealDirect`
reports N (`XAD/Interface.hpp:59`). `xad::fwdd` is the safer carrier.

⚠️ **TF24f's fixed point *is* `dprofit == 0`** (`plant/src/tf24f_strategy.cpp`, plant
#525/#529) — the same equation from the other side. plant #529 claims TF24f's `k→∞`
limit is TF24's optimisation; before 11a that was **false**, the two converging on
points ~**7e-05 MPa** apart with no `k` reconciling them. After 11a they agree.

## 11b. One body for the model and its derivative, and the ci tolerance

The `namespace detail` AD replicas are gone, replaced by scalar-generic member
templates that the model and its derivative both instantiate — bit-identical on the
double path at five entry points × 11 operating points. The guide carries the
amplifier table and the 1e-10 decision. What is only here:

⚠️ **Contraction follows the call graph.** Unifying a replica with the function it
mirrors changes results even when the algebra is identical, because FMA contraction
depends on which expressions end up in one translation unit's inlined body. The
reassociation causes were named: `et/4*((ci−g)/(ci+2g))`, and `pow(s,2)` versus
`s*s`.

**Ordering mattered more than either change.** Tightening `psi_stem_to_ci` *first*
made the replica unification move results by **4.98e-07 → 3.51e-10, a 1400×
reduction** — because the tolerance was the amplifier.

**The tolerance not taken:** 1e-13 costs **+6.2%** (2.73 → 2.90 µs) and is still 17%
faster than golden section's 3.51. Recorded so the next person does not re-derive the
curve.

⚠️ **A false correctness claim, corrected.** The unification was *not* more accurate
— **40 of 82 rows improved, median factor 0.999**, with roots differing by ~1e-15
MPa. It is one body, not a better one.

⚠️ **The kernels are templated on their ARGUMENT, not on traits.** Every trait is a
plain `double` member and there is no templated member storage, so 11b delivered
exact derivatives with respect to **state** only. That is why #4 needed 11c–11f.

## 11c. What measuring first established

⚠️ **The highest-value numbers in this item, and they exist nowhere else — no code
comment leads a reader back here.**

The mixed partial is well conditioned: `M = ∂²profit/∂ψ∂θ` is **stable to 7
significant figures across five decades** of step (relative 1e-7 to 1e-3), because ψ
is an *argument* of `dprofit`, so no argmax is involved. `H` runs −1.56 to −61.4. The
composite agrees with a resolved slope at **ratio 0.9979–1.0000 for all eight
traits** across five soil-moisture levels.

**How much of `dA/dθ` is argmax-mediated** — which is why the second term is not a
correction:

| trait | indirect share |
|---|---:|
| `cost_scale_TF24`, `beta2`, `stem_b`, `stem_c` | **100%** |
| `jmax_25` | 61% |
| `vcmax_25` | 52% |
| `root_b` | 6% |
| `beta_R_H` | −15% |

⚠️⚠️ **The active set is the crux, not the mixed partial.** At a pinned optimum the
constraint binds, `dprofit ≠ 0` there, and `H` is the wrong denominator — so the bare
composite fails silently:

| point | composite / truth |
|---|---:|
| pinned wet | **3.5e+07**, 5.3e+06, −8.5e+03 |
| pinned dry | 1.65, 0.93 |
| interior | 1.0000 |

The true pinned-wet gradient is ~**1e-08** while the composite returns O(1). **A
plain finite difference does not have this failure**, because it differences the
constrained answer. 42 of 240 rows are pinned. This table is what `R/gradient.R`'s
guard and the 42-pinned test both rest on.

⚠️ Achievable precision of any of this: **~1e-9**, set by `psi_stem_to_ci`.

## 11d. Implicit function theorem, forward AD, no tape

**The motivation, measured on the fit:** 10.1 minutes, **10,045 model evaluations,
97.6% of them finite-difference gradients** — 245 gradients × 40 parameters over 16
species (Sabot et al. 2022).

**The costing, which #4's remainder is still sized against:**

| route | µs per gradient | vs today |
|---|---:|---:|
| 2N re-solves (what the fit did) | 288.6 | — |
| IFT: one solve + 2N `dprofit` | 24.1 | **12×** |
| parameter-explicit kernels | 4.0 | **72×** |

⚠️ **Read 11e before quoting either figure.** Both were denominated against the
solve, and no caller spends its time there.

**Reverse mode is not needed**, independent of item 5's linking problem: forward AD
is O(N), vector and scalar were within noise at N=40 (**0.164 vs 0.163 µs**), and
both were slower than a central difference at N ≥ 12.

⚠️ **Three measurements in this item were wrong before they were right, all from the
work not actually running:** a stale `src/RcppR6.o` (R does not track header
dependencies), a stale `bench_solve` after a `sed`, and dead-code elimination
hoisting a benchmark's inner loop out of the timing region — where FD and forward
both timed at exactly 0.000 µs. **Vary the input per iteration, consume the output,
and treat a suspiciously round zero as a broken harness rather than a result.**

**Still unmeasured:** FMA contraction's contribution (1.19e-06 under golden section,
verified by rebuilding both sides with `-ffp-contract=off` and getting bit-identity —
probably much smaller now); hazard 3 on the plant side, blocked on plant #591;
whether the inner root-finds need an IFT of their own.

## 11e. Stage 1, and the cost model that was wrong three times

`leaf_gradient()` and `set_traits()` shipped. The active-set premise is **tested, not
assumed**: the implied Newton step `|dprofit(ψ*)/H|` is a distance in MPa needing no
scale of its own, and over the 288-point grid the two populations are five orders
apart — **worst interior 4.8e-11, mildest pinned 6.3e-06** — so the 1e-08 default
sits in an empty band. Classification: **198 interior / 42 pinned / 48 no-gradient**,
reproducing 11c exactly. A second guard, found by measurement: at a pinned optimum
ψ* sits 1e-06 of a bracket width from its bound, so the step in ψ cannot be centred
without being clamped — which catches all 42 pinned and all 48 shut-down rows on its
own. ⚠️ Not a substitute for the stationarity test: it fires only when the feasible
interval is narrow.

⚠️ **In C++ the composite wins; in R it loses. Say which layer a figure belongs to.**

| per parameter | FD | IFT | |
|---|---:|---:|---:|
| the **11** traits that rebuild no spline | 76.7 µs | **17.5 µs** | **4.39×** |
| `stem_b`, `stem_c`, `root_b`, `root_c` | 168.7 µs | 146.1 µs | 1.15× |
| all 15 | 245.4 µs | 163.6 µs | 1.50× |

**Quote the split, never the 1.50× total.** From R the composite measured **6%
slower** than the finite difference, because a call costs ~1.8 µs against the 0.26 µs
of C++ work it wraps.

**The components, because 11d's table omitted two of them:**

| component | µs | 11d |
|---|---:|---|
| `find_root_collar_psi` | 6.22 | counted |
| `dprofit_droot_collar_psi` | 0.49 | counted at 0.26 — that figure had the caches already seated |
| **`evaluate_root_collar_psi`** | **0.88** | **not counted at all** |
| `set_physiology` | 0.13 | not counted |
| `set_traits`, no rebuild | 0.02 | not counted |
| **`set_traits`, spline rebuilt** | **21.8** | **not counted, and 3.5× a solve** |

The missing one is the composite's *other half*, the direct term `∂A/∂θ|_ψ`:
`dA/dθ = ∂A/∂θ|_ψ + (∂A/∂ψ)(dψ*/dθ)`, and dropping the first term is not an
approximation but a different quantity. **Any future costing must carry both.**

⚠️ **Twelve R calls to report a solve.** Reading the twelve outputs through their own
active bindings cost ~**15 µs** against a 2.8 µs solve — five times the model.
`Leaf::operating_point_values()` returns them in one call. That is the number to
reach for when someone proposes exposing more state to R one field at a time.

⚠️ **And the trap in #44's step rule.** The traits' rule floors the step at 1, which
is right for parameters spanning 0.3 to 9.4e3 and is a **3% perturbation** of
`leaf_specific_conductance_max` at 3.14e-05 — a secant, not a derivative. Those
parameters take a relative step, and the test pins step-*independence*, which is the
assertion that catches this class of error.

⚠️ **Size a cost model against the caller, not against the model.** Three instances
in this item: 11d's 12× was computed against the solve; the R retraction was
denominated in the wrong currency; and stage 3 was promoted for a hydraulic
calibration that turns out to hold `stem_b` and `stem_c` **fixed** — they are data,
derived from each species' measured P50 and P88 — while fitting two parameters that
are not traits at all.

## 11f. The homogeneity result

A gradient in `stem_b` needs no spline rebuild: **35.5 → 1.45 µs**, which puts it in
the same class as `vcmax_25` at 1.64.

**What the rebuild was:** seeding 101 knots costs **11.86 µs of incomplete gammas**,
each interpolator `init` 3.08 µs, and rescaling the knots **0.013 µs**. The cost was
the gamma function, not the spline machinery — and decomposing it is what found the
answer.

**The identity.** `G(ψ; b, c) = b·g(ψ/b; c)`, so `b` enters only as a scale on both
axes:

> **G(ψ; s·b, c) = s · G(ψ/s; b, c)**

The knot grid scales with it too (`ψ_max = b·log(100)^(1/c)`), so this holds for the
**spline** and not merely for the integral it approximates — measured, a rescaled
spline reproduces a rebuilt one to **0–3e-16** at factors from 1.001 to 1.2. The two
derivative accessors carry no factor, since differentiating `s·G(ψ/s)` in ψ cancels
it.

⚠️⚠️ **A closed-form version was built first and rejected. Read this before
rebuilding it.** Reading `G` from `(b/c)·γ_lower(1/c, (ψ/b)^c)` costs 0.117 µs
against 11.9 and covers `stem_c` too. It disagrees with the rebuild route by a
systematic **3.5e-3**, where the rescale agrees to 3e-9–9e-5. **The tell was that
both routes were step-stable across five decades and converged to *different*
values** — 2.62930 against 2.62850. Noise moves with the step; a model mismatch does
not.

**The rule that follows, and it generalises: a gradient must differentiate the model
being evaluated, not the model the evaluated one approximates.** A *more* accurate
perturbed evaluation makes the gradient worse. Here it also mixed routes inside the
composite — `M` from the closed form, `H` from the base state's spline.

**It is quieter as well as cheaper.** A rebuild reseeds 101 gammas per side, so the
two sides of a central difference carry uncorrelated rounding, which dividing by
`h ≈ 4e-6` amplifies: the rebuild route jitters up to **9e-05** between neighbouring
steps — at `step = 1e-06` returning 2.6284252697 where its neighbours agree on
2.6284982855 — while the rescale route stays on the plateau. That is why the
equivalence test's tolerance is loose, and it should not be tightened.

⚠️ **From R it is 1.23×**, not 24.5×, because `leaf_gradient()` pays 204 µs per call
to construct a `Leaf`.

## 15. Housekeeping

CI, the C++ suite under `R CMD check`, and rendered C++ API docs all landed in #20;
the guide's build section is the live description. Two decisions worth keeping:

⚠️ **`tests/cpp.R` compiles the two sources directly rather than calling `make`**,
because `R CMD check` warns about the GNU extensions in `tests/cpp/Makefile`. The
sanctioned fix, `SystemRequirements: GNU make`, was tried and reverted: installing
this package needs no make, so every `LinkingTo: phylloptim` consumer would inherit a
declaration that is **false for them**.

**Doxygen renders the C++ API through `tools/doxygen_filter.awk`**, which converts
plain `//` comments at render time so that ten headers did not have to be rewritten —
and CI asserts every non-comment line survives the filter byte-for-byte, because a
bug there would silently change what Doxygen says the code does. Publishing is opt-in
through the repo variable `PUBLISH_DOCS=true`.
