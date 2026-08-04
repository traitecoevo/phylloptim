# Seed: finish trait gradients, forward mode (#4, PLAN 11d)

Handoff brief for starting the trait-gradient work in a fresh session. Everything
below was established while landing #36 and writing PLAN 11a–11d. **Read this
first so you don't re-derive it — three of the conclusions here reverse an earlier
one, and two of them were reached only after a measurement went wrong.**

`PLAN.md` is the authority for reasoning; items 11a–11d are the relevant ones.
This file is the entry point for this particular next step.

## Where the repo is, as of 2026-08-04

`master` is `748a712`. The substantive code is `d22907a` (PR #36), which landed
three things:

| | what | blast radius |
|---|---|---|
| 11a | collar solve solves `dprofit == 0` (TOMS748), not golden section | 1.5e-03 |
| 11b | `psi_stem_to_ci` tolerance 1e-7 → **1e-10** | 1.82e-07 |
| 11b | `namespace detail` AD replicas **deleted** for scalar-generic member templates | 3.51e-10 |

Plus `dprofit_droot_collar_psi(psi, bool* feasible = nullptr)`.

**Open issues: #1, #3, #4 (this), #6, #7, #28, #31, #33, #34, #38, #39, #40, #41.**
#38–#41 were filed *from* the calibration project and are the requirements it hit.

## The goal

Supply `dA/dθ` (and `dgs/dθ`, `dψ_leaf/dθ`) with respect to **traits**, so
[`leaf-calibration`](../../leaf-calibration) can fit with real gradients. Its
hierarchical fit currently spends **97.6% of 10,045 model evaluations** on
finite-difference gradients — 245 gradients × 40 parameters, 10.1 minutes.

## ⚠️ Decisions already made. Do not reopen without new measurements.

**1. Implicit function theorem, forward mode, no tape.** PLAN 11d. The win is
architectural, not modal: the cost is *re-solving the model 2N times per
gradient*, and IFT removes that.

**2. Reverse mode is NOT wanted, and item 5's no-linking decision was revisited
and stands.** `xad::adj<double>` does not link header-only — the `Tape`
implementation is compiled once in `odelia/src/Tape.cpp`, which says so in a
comment, and consumers link its `.so`. leaf's CMake and Python consumers cannot.
Only N = 1 tapes are instantiated, so `adj<double,N>` would not link even with the
library. **And it is unnecessary**: after IFT, differentiation is ~10% of the cost,
so reverse mode's O(1)-in-N advantage applies to a tenth of the work.

**3. `Leaf<T>` is not on this path.** Its only remaining justification is plant
#537's cut-point rule, which concerns an outer AD pass over *plant*. Decouple them:
#4 can close without it. The issue's title no longer describes what it is for.

**4. Vector forward mode buys nothing here.** Measured on the assim kernel,
`fwd<double,40>` and 40 scalar sweeps are within noise (0.164 vs 0.163 µs), and at
N ≥ 12 both are *slower* than central FD on the same kernel. Forward AD is O(N)
just like a difference quotient. An earlier note claimed vector forward was "what
gets the speedup"; it is not.

## The maths to implement

At the solved operating point, with ψ the root-collar potential and θ a trait:

```
H = ∂²profit/∂ψ²          curvature at the optimum; measured -1.56 to -61.4 over
                          the golden grid, never near zero
M = ∂²profit/∂ψ∂θ         mixed partial, ψ held FIXED while θ moves

dψ*/dθ = -M / H

dA/dθ  = ∂A/∂θ|_ψ  +  (∂A/∂ψ)(dψ*/dθ)
         ^direct       ^through the optimum
```

Why this is cheap: `dprofit_droot_collar_psi` is **exact in ψ and smooth in θ**
(ψ is an argument, so no argmax is involved), and it is closed form — the transport
derivative is analytic via the splines' own `.deriv()`, and the `ci` derivative
comes from IFT on the demand-supply residual rather than from re-running TOMS748.
So `dprofit` costs **0.256 µs against 3.607 µs for a full solve — 7%**.

`M` by a plain difference quotient in θ is stable to **7 significant figures across
five decades** of step size (relative 1e-7 to 1e-3). It does not need AD.

⚠️ **The second term is not a correction.** For `cost_scale_TF24`, `beta2`,
`stem_b` and `stem_c` it is **100%** of `dA/dθ`; for `vcmax_25` 52%. A finite
difference on the solved output dropped it silently before 11a.

## The staged plan

Costs are per trait gradient at N = 40, measured. R overhead is ~15.5 µs per call
regardless, which is why stage 2 exists.

| stage | where | µs | vs baseline | effort |
|---|---|---:|---:|---|
| baseline — FD on the solve | (what the fit does) | 1528 | — | — |
| **1** | R layer, composition in R | **153** | **10×** | small |
| **2** | composition in C++, one R call | **40** | **39×** | medium |
| **3** | parameter-explicit kernels + forward AD | **19** | **78×** | medium |

**Stage 1 first.** It needs no C++ change at all: `evaluate_root_collar_psi` and
`dprofit_droot_collar_psi` are already bound, and R's per-call overhead on
`dprofit` turns out to be only 1.7 µs, so the composition in R keeps most of the
benefit. It validates the active-set branch and the composite against the
references below, and unblocks the calibration project immediately.

Then stage 2 moves the loop into C++ behind **one** R call per operating point —
which is the same lesson as #39: R overhead dominates, so batch it. Stage 3 replaces
the difference quotient in θ with AD on parameter-explicit kernel overloads, and is
the only stage needing new C++ derivative code. Measure before doing it; the
marginal gain is ~2×.

## ⚠️⚠️ The one hazard that will silently ruin this: the active set

**Stationarity is the premise of the entire derivation, and it fails at a
constrained optimum.** When the optimum is pinned to a bracket bound, `dprofit ≠ 0`
at the answer, `dψ*/dθ` is the *bound's* derivative rather than `-M/H`, and `H` is
not the right denominator either.

The formula does not fail loudly. Measured:

| operating point | class | trait | composite / truth |
|---|---|---|---:|
| ψ_soil=4, vpd=2, 5 layers | pinned wet | `stem_b` | **3.5e+07** |
| ψ_soil=3, vpd=4, 3 layers | pinned wet | `stem_b` | 5.3e+06 |
| ψ_soil=4, vpd=2, 5 layers | pinned wet | `vcmax_25` | −8.5e+03 |
| ψ_soil=4, vpd=0.5, 3 layers | pinned dry | `vcmax_25` | 1.65 |
| ψ_soil=2, vpd=2, 1 layer | interior | either | **1.0000** |

At the pinned-wet rows the true gradient is ~1e-08 and the composite returns
**O(1)** — plausible-looking and wrong by seven orders of magnitude. **42 of 240
feasible golden-grid rows are pinned, all at ψ_soil 3–4**, i.e. the dry end a
calibration in drought will visit.

A finite difference on the solved output does *not* have this failure, because it
differences the constrained answer. That is the one thing the crude method does
better, and it is why an explicit active-set test is not optional.

**How to detect it:** the classification `maximise_profit_over_collar` already
makes internally — sign of `dprofit` just inside the wet end and at the dry bound.
Consider exposing it rather than reimplementing it. At a pinned point either return
the bound's own derivative or fall back to differencing the solve.

## Verification references

Measured at ψ_soil = 2.0, PPFD = 900, VPD = 2.0, one layer, default
`leaf_traits()`. `H = −8.9561` there. Ratios are against a least-squares slope over
±2% at n = 41, which agreed to 0.9979–1.0000.

| trait | `dψ*/dθ` | `dA/dθ` | indirect share |
|---|---:|---:|---:|
| `vcmax_25` | 1.7459e-03 | 1.7209e-02 | 52% |
| `jmax_25` | 1.0862e-04 | 9.1320e-04 | 61% |
| `cost_scale_TF24` | −7.6704e-02 | −3.9520e-01 | 100% |
| `beta2` | −4.2610e-02 | −2.1954e-01 | 100% |
| `stem_b` | 5.5247e-01 | 2.8465e+00 | 100% |
| `stem_c` | −1.6390e-01 | −8.4448e-01 | 100% |
| `root_b` | 2.6656e-04 | 2.3895e-02 | 6% |
| `beta_R_H` | 8.1973e-06 | −2.8802e-04 | −15% |

Also verified across ψ_soil ∈ {0.5, 1, 3, 4} for `vcmax_25` and `stem_b`, ratios
0.9999–1.0002. **Test the pinned rows separately** — the table above is all interior.

## Three traps that cost time in the work leading here

All three produced plausible wrong numbers, and all three were the work not
actually running.

1. **R does not track header dependencies.** `R CMD INSTALL` after editing
   `inst/include/` reuses a stale `src/RcppR6.o`, and the R layer goes on running
   the *old* model. Two rounds of diagnostics were taken against the previous solver
   before this was noticed. **`rm -f src/*.o src/*.so` first**, then check one value
   against the C++ suite.
2. **Benchmarks: build both arms in the same tree, and interleave.** A seven-point
   sweep built in one loop in a scratch tree implied +1.5% where the controlled A/B
   said +3.4%. Code layout moves `bench_solve` by ~2%, the same order as the
   effects being measured.
3. **Dead-code elimination will hoist a micro-benchmark's inner loop out of the
   timed region.** FD and scalar forward both timed at exactly `0.000 µs`, which is
   the tell. Vary the input per iteration and consume the output.

## What is downstream

- **#6 / `leaf-calibration`** is the customer and is already running with FD. Its
  `analysis.qmd` has a section headed "Nothing here uses a gradient of the model"
  that this work is meant to retire. It fits **three** responses — `A`, `gs` and
  `ψ_leaf` — so the composite is wanted for each, not just for `A`.
- **#41** (`R_d` unbound) and **#39** (`leaf_solve()` R overhead) are on that
  project's own next list and are independent of this.
- The concept is written up for a general audience in
  [overstorey](https://github.com/traitecoevo/overstorey) —
  `posts/2026-08-04-differentiate-the-equation/`.
