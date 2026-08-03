# leaf

A header-only C++ leaf gas-exchange model in which **stomatal behaviour emerges
from hydraulics** instead of an empirical conductance function.

Farquhar-von Caemmerer-Berry photosynthesis is coupled to an explicit
soil → root → stem → leaf water transport path — Weibull vulnerability curves
for both xylem and roots, multi-layer soil, per-layer root resistance,
gravitational head — and the operating point is chosen by Sperry-style gain-risk
**profit maximisation** over the root-collar water potential. Forward-mode
automatic differentiation (XAD) supplies exact derivatives of profit with respect
to that potential, for models that need to track acclimation.

A full solve costs about **4 µs**, which is what makes it usable inside a
demographic model that calls it millions of times.

This code was developed as the TF24 strategy inside
[traitecoevo/plant](https://github.com/traitecoevo/plant) and is extracted here
so it can be tested, profiled, extended and embedded on its own.

## Why a separate package

**A home for several stomatal models, not just ours.** The package carries our
hydraulic gain-risk formulation, and it already contains two alternatives — the
Sperry et al. (2017) cost formulation and the Medlyn et al. (2011) optimal
stomatal model — inherited from plant. Today those are second-class:
`optimise_psi_stem_Sperry` is hardwired to a single soil layer, and the Medlyn path
bypasses the hydraulic solve altogether, so neither can be swapped in as a
like-for-like alternative. The goal is to make each one a **first-class member**,
selectable and runnable against identical drivers, alongside Prentice et al. (2014)
least-cost and Cowan-Farquhar.

The right way to do that turns out to be one level deeper than swapping cost
functions. Those models all maximise a profit, so they all satisfy the same
first-order condition `dA/dE = λ` and differ **only** in the function λ(state) —
the marginal cost of water. So what should be pluggable is λ. Six models become
six small functions sharing one tested numerical core, which makes a comparison
apples-to-apples by construction rather than by careful bookkeeping. This is the
central result of a companion manuscript (see [PLAN.md](PLAN.md) item 14), and it
is something none of the existing R packages can support, because each commits to
a single hydraulically explicit scheme or to none. See [PLAN.md](PLAN.md) item 7a.

The same refactor applies from the other side, to the **water supply path**. The
gas-exchange core is already entirely soil-agnostic — the multi-layer soil and
root system enter the solve only as a single supply function `E_up = f(P_collar)`,
so pulling them behind an interface would let the multi-layer root system be
swapped for a single soil water potential. That both lowers the barrier for a
bare-leaf user (no root-mass profile to construct) and is what makes comparison
*fair*: the alternative formulations worth comparing against are all written for
one ψ_soil, so you have to be able to hold the supply side fixed. [PLAN.md](PLAN.md)
item 7b.

**Fast and differentiable, so it is built for calibration.** A full hydraulic
solve costs ~4 µs, and forward-mode AD (XAD) gives exact derivatives rather than
finite differences. That combination is what calibration wants: gradient-based
optimisers and Hamiltonian samplers need many evaluations *and* clean gradients,
and finite-differencing a nested root-find through a golden-section search is
exactly the case where numerical gradients are noisiest. Two honest caveats: the
AD currently differentiates with respect to the collar potential only, not with
respect to traits, and getting trait gradients needs the templated `Leaf<T>` of
[PLAN.md](PLAN.md) item 11; and there is no worked example yet — a calibration
vignette is item 12 and is the demonstration this claim needs.

## Status

**v0.0.1 — early.** The model itself is mature and in production use inside
plant. The *packaging* is new, and two things are not done yet:

- The values this build produces have **not yet been cross-checked against
  plant's compiled build**. That is job one; see [PLAN.md](PLAN.md).
- There is **no R-level API**. This is a headers-only package, consumed from C++.

It also carries one inherited defect worth knowing about before you use it: on
the hydraulic-shutdown path, transpiration, assimilation and soil water uptake
are left holding the previous solve's values (plant #578). `main` reproduces
plant's behaviour exactly, including this; the fix is on the `feature/api-cleanup`
branch alongside the other changes that alter results or the API.

## Use from C++

There is nothing to link against — one include is the whole library.

```cpp
#include <leaf.hpp>

leaf::Leaf l;                   // default Eucalyptus saligna traits
l.setup_transpiration(100);     // build the xylem vulnerability splines
l.setup_root_vulnerability(100);

std::vector<double> psi_soil{2.0};             // positive suction, MPa
std::vector<double> soil_depth{1.0};           // m
// kg C per m2 LEAF, not absolute carbon: the leaf is purely intensive.
std::vector<double> root_carbon_per_leaf_area{20.0};

l.set_physiology(root_carbon_per_leaf_area, /*PPFD*/ 900,
                 psi_soil, soil_depth,
                 /*leaf_specific_conductance_max*/ 3.14e-5,
                 /*atm_vpd*/ 2.0, /*ca*/ 40.0,
                 /*leaf_temp*/ 25.0, /*atm_o2_kpa*/ 21.0, /*atm_kpa*/ 101.3);

l.find_root_collar_psi();       // solve

l.opt_psi_stem_;                // leaf water potential at the optimum, MPa
l.opt_root_psi_;                // root-collar potential, MPa (positive magnitude)
l.assim_colimited_;             // A, umol m-2 s-1
l.transpiration_;               // E, kg H2O m-2 s-1
l.stom_cond_CO2_;               // gc, mol CO2 m-2 s-1
l.profit_;                      // A - hydraulic cost
l.soil_consumption_;            // per-layer water uptake
```

Compile with C++20 and three include paths — this package, `odelia`, and Boost:

```sh
c++ -std=c++20 -O2 \
  -I /path/to/leaf/inst/include \
  -isystem /path/to/odelia/inst/include \
  -isystem /path/to/boost \
  my_program.cpp -o my_program
```

## Use from R

As a headers-only dependency, the way `BH` is used. Add to your DESCRIPTION:

```
LinkingTo: BH, odelia (>= 0.2.0), leaf
```

`LinkingTo` is **not** transitive in R, so you must name `BH` and `odelia`
yourself even though it is `leaf` that includes them — including the odelia
version, for the same reason.

There is no R-level interface to `leaf::Leaf` yet — see [PLAN.md](PLAN.md).

## Dependencies

Deliberately few, and all header-only:

| | why | how |
|---|---|---|
| **odelia** (>= 0.2.0) | cubic-spline interpolator for the pre-integrated vulnerability curves, and the vendored **XAD** automatic-differentiation library | `LinkingTo` |
| **BH** (Boost) | TOMS748 root finder, incomplete gamma for the closed-form vulnerability integral | `LinkingTo` |

Nothing else, and **neither dependency needs R**. The leaf model itself does not
touch **Rcpp** or the R C API: `leaf/util.hpp` replaced plant's `util::stop`
with a plain `std::runtime_error` and `NA_REAL` with a quiet NaN. odelia's
solver core was the last R touchpoint in the include graph, via `ode_util.hpp`;
that was removed upstream in traitecoevo/odelia#44, so the test suite now builds
against the real headers with nothing standing in for R at all. **odelia 0.2.0 is
the first release with that fix**, hence the version requirement: an older odelia
would otherwise fail deep in the build with `RcppCommon.h: No such file or
directory`, which does not point at the cause.

Both dependencies are already required by plant, so plant pays nothing new for
depending on this package.

## Tests

The test suite is plain C++ and needs neither R nor a test framework:

```sh
make -C tests/cpp
```

It discovers BH and odelia through `Rscript` if R is installed, and otherwise
falls back to a sibling `odelia/` checkout and Homebrew Boost. Override with
`make BH_INC=... ODELIA_INC=...`.

`R CMD check` runs the same suite, compiled with R's own configured compiler
against the installed headers — so a package that `LinkingTo`s this one finds out
from its own check when a header stops compiling.

## API documentation

```sh
doxygen        # docs/html/index.html
```

Doxygen rather than roxygen because roxygen documents R objects and this package
has none. The headers' comments are the substantive documentation here, and
`tools/doxygen_filter.awk` presents them to Doxygen without modifying a single
source file.

## How this compares to other leaf models

See [COMPARISON.md](COMPARISON.md) for a feature-by-feature comparison against
`plantecophys`, `bigleaf` and `tealeaves`. The short version: those packages are
stronger on empirical stomatal models, leaf energy balance and fitting to
measured data; this one is the only one with an explicit hydraulic architecture
and a profit-maximisation solve, and the only one written to be embedded in a
larger model.

## Licence

AGPL (>= 3), inherited from plant.
