# Trait gradients (#4): what this note seeded, and what came back different

⚠️ **SPENT. Stage 1 landed — `leaf_gradient()` and `set_traits()` are on `master`.
[`PLAN.md` item 11e](../PLAN.md) is the record and the authority; read that first.**

This file is kept rather than deleted because **its central quantitative claim was
wrong, and the shape of the error is worth carrying forward**. It projected a 10×
speedup for stage 1 in the R layer. Measured, stage 1's gradient composite came back
**6% slower** than the finite difference it was built to replace, and the 4× that
does exist came from something this note did not mention at all. Everything below is
marked for what survived.

## ✅ What survived

**The maths.** At the solved operating point, with ψ the root-collar potential:

```
H = ∂²profit/∂ψ²          curvature at the optimum; measured -1.56 to -61.4
M = ∂²profit/∂ψ∂θ         mixed partial, ψ held FIXED while θ moves

dψ*/dθ = -M / H
dA/dθ  = ∂A/∂θ|_ψ  +  (∂A/∂ψ)(dψ*/dθ)
```

Implemented exactly as written. `dprofit_droot_collar_psi` is exact in ψ and smooth
in θ, so `M` needs no AD — a plain difference quotient is stable to 7 significant
figures across five decades of step. The second term is **100%** of `dA/dθ` for
`cost_scale_TF24`, `beta2`, `stem_b` and `stem_c`, and 52% for `vcmax_25`.

**The verification references**, at ψ_soil = 2.0, PPFD = 900, VPD = 2.0, one layer,
default `leaf_traits()`, `H = −8.9561`. These are now a test
(`tests/testthat/test-gradient.R`) and the composite reproduces all eight to 4
digits:

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

**The active-set hazard, which was the right thing to be worried about.** It was
this note's headline warning and it earned the space: 42 of 240 feasible golden rows
are pinned, and the bare composite is wrong there by up to 3.5e+07. What shipped
tests the premise rather than assuming it — the implied Newton step `|dprofit/H|`,
whose two populations are five orders of magnitude apart over the grid — and falls
back to differencing the solve. The suggestion to *reuse* the solver's own
classification was **not** taken: an independent stationarity test needs no new C++
and does not couple the R layer to the solver's internal branch structure.

**All three traps.** Every one of them cost time again in this round:

1. **R does not track header dependencies.** `rm -f src/*.o src/*.so` before
   reinstalling, every time.
2. **Build both benchmark arms in one tree and interleave.** Done; the arms differ
   by 6%, and round-to-round spread is ~0.5%, so this was load-bearing.
3. **Dead-code elimination, and in R its cousin.** The first timing harness written
   for this work used `f <- function() expr` and measured a **lazily-evaluated
   promise that fires once** — reporting 0.07 µs for everything, including a
   204 µs object construction. Same tell as the C++ case: a suspiciously uniform
   number. Vary the input per iteration and consume the output.

## ❌ What was wrong

**1. "Stage 1 needs no C++ change at all."** It needed one, for a reason this note
did not consider: **the traits are not settable from R**, so perturbing θ meant
constructing a new `Leaf`, at **204 µs** — 73 solves. The note's 153 µs estimate
assumed a θ-perturbation costs one `dprofit` call (1.7 µs). `set_traits()` was
added to fix that, and it could not be fifteen settable fields, because
`set_physiology` derives `vcmax_`/`jmax_`/`R_d_` behind a cache keyed on
`(leaf_temp, atm_o2_kpa)` alone — so "change the trait, then set the drivers again"
silently does not recompute them.

**2. "Stage 1: 153 µs, 10× — R's per-call overhead on `dprofit` is only 1.7 µs."**
The premise is right and the conclusion does not follow. 1.7 µs of R overhead per
call is **seven times** the 0.26 µs of C++ work `dprofit` does, so a route that
trades one solve for one gradient evaluation *plus an extra R call per side* comes
out behind. Measured per parameter at one operating point: FD with a fresh `Leaf`
343 µs, FD with reuse **84 µs**, the IFT composite **89 µs**.

**3. The whole staged table's currency.** Stages 1–3 were costed in µs of C++ work.
Of the 84 µs per parameter that stage 1 actually costs, **~6 µs is C++** and the
rest is the R boundary. The table is not so much wrong as denominated in the wrong
thing, and the correction is #39's: batch the calls.

**4. "Unblocks the calibration project immediately."** Overstated. `leaf_predict()`
builds one `Leaf` per species and loops ~32 observations, so it amortises the
construction this work removed. Measured on that shape, `set_traits` is worth
**26%**, not 4× — and three quarters of what remains is still the R boundary. #6's
next win is #39, not a better gradient.

## Still open, unchanged by this round

- **#6 / `leaf-calibration`** is the customer. It fits **three** responses — `A`,
  `gs` and `ψ_leaf` — and all three are differentiated, plus `ψ*` itself.
- **#41** (`R_d` unbound) and **#39** (`leaf_solve()` R overhead) are independent of
  the gradient work — and #39 is now the higher-value one of the two.
- **`Leaf<T>` is still not on this path.** 11d's reasoning stands and stage 1 closed
  without it.
- **Reverse mode is still not wanted.** `xad::adj<double>` does not link
  header-only; and after IFT, differentiation is a few percent of the cost.
- The concept is written up for a general audience in
  [overstorey](https://github.com/traitecoevo/overstorey) —
  `posts/2026-08-04-differentiate-the-equation/`. ⚠️ **If that post quotes the 12×
  or the 10×, it needs the correction in 11e.**
