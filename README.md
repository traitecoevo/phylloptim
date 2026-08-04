# leaf

A C++ leaf gas-exchange model, callable from R, in which **stomatal behaviour
emerges from hydraulics** instead of an empirical conductance function.

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

**v0.1.0 — early, but validated.** The model itself is mature and in production
use inside plant. The *packaging* is what is new.

- **Cross-checked against plant's compiled build**, and the swap was bit-identical
  at the point it was made: plant's full suite 0 fail / 0 error on both builds,
  and the SCM regression identical across 78/78 nodes. The 1-ULP disagreement
  that held this up turned out to be R's decimal parser rather than either model.
- **The shutdown defect is fixed.** On the hydraulic-shutdown path, transpiration,
  assimilation and uptake used to be left holding the previous solve's values
  (plant #578). That, and three further stale-state exits ported from plant #585,
  are all fixed here.
- **Results now differ from plant's own leaf, deliberately** — see
  [NEWS.md](NEWS.md). The most consequential single change is deriving the
  ppm→Pa conversion from the actual atmospheric pressure instead of a hard-coded
  101.3 kPa, which moves TF24 offspring production by 2.4% in plant, because
  plant's driver default is 100.5.

## Two ways in, and the model does not need R

The model is a set of self-contained C++ headers under `inst/include`. They use
no R and no Rcpp, depend only on Boost and the header-only parts of odelia, and
compile and run with no R installed — so the same model is available to a C++
program, to a Python extension, and to R, with none of those paying for the
others. The R layer (`src/`, `R/`) sits on top of those headers and is never
included by them; the dependency runs one way only. See [PLAN.md](PLAN.md) item
6a for the decision, and `.github/workflows/cpp-tests.yml` for what enforces it —
it builds the whole C++ suite on a runner with no R on it.

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

Or use the CMake package, which handles those three paths for you. `odelia` is
distributed as an R package but the headers used here are plain C++, so a git
checkout of it is enough — nothing needs installing or building:

```sh
git clone https://github.com/traitecoevo/odelia
cmake -B build -DLEAF_ODELIA_INCLUDE_DIR=$PWD/odelia/inst/include
cmake --build build
ctest --test-dir build          # runs the C++ suite, including the golden file
cmake --install build --prefix /usr/local
```

```cmake
find_package(leaf REQUIRED)
target_link_libraries(my_program PRIVATE leaf::leaf)
```

`leaf::leaf` is an INTERFACE target — headers, an include path and `cxx_std_20`,
with nothing to link. `add_subdirectory(leaf)` works the same way if you would
rather vendor it.

⚠️ **Build with optimisation on if you intend to compare against the golden
file.** Measured on macOS/arm64: `-O1`, `-O2` and `-O3` all reproduce
`tests/cpp/golden/operating_points.tsv` bit-for-bit and `-O0` does not, missing
by about 13 ULP because it declines to contract `a*b + c` into an FMA. A debug
build that fails `test_golden` by ~1e-15 has found nothing. The CMake build
therefore defaults to `Release` rather than to CMake's flagless default.

## Use from Python

The same headers, through [pybind11](https://pybind11.readthedocs.io/) — no R
anywhere in the picture. A minimal extension module:

```cpp
// pyleaf.cpp
#include <leaf.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

PYBIND11_MODULE(pyleaf, m) {
  py::class_<leaf::Leaf>(m, "Leaf")
      .def(py::init<>())
      .def("set_physiology", &leaf::Leaf::set_physiology)
      .def("find_root_collar_psi", &leaf::Leaf::find_root_collar_psi)
      .def_readonly("profit", &leaf::Leaf::profit_)
      .def_readonly("opt_psi_stem", &leaf::Leaf::opt_psi_stem_)
      .def_property_readonly("g1_eff", &leaf::Leaf::g1_eff);
}
```

```cmake
find_package(leaf REQUIRED)
find_package(pybind11 REQUIRED)
pybind11_add_module(pyleaf pyleaf.cpp)
target_link_libraries(pyleaf PRIVATE leaf::leaf)
```

```python
>>> import pyleaf
>>> l = pyleaf.Leaf()
>>> l.set_physiology([20.0], 900, [2.0], [1.0], 3.14e-5, 2.0, 40.0, 25.0, 21.0, 101.3)
>>> l.find_root_collar_psi()
>>> l.profit
2.5158434693939866
```

That value is the golden file's `profit` at this operating point, to the last
bit — which is the point of the example. `util::stop` throws
`std::runtime_error`, so pybind11 turns the model's input validation into a
`RuntimeError` with no extra work.

`std::vector<double>` needs `pybind11/stl.h`, as above; swap it for
`pybind11/numpy.h` and an `Eigen`-style binding if you want the soil profile to
arrive as an array without a copy.

## Use from R

Drivers in, operating point out. `leaf_solve()` is vectorised, so a response
curve is one call:

```r
library(leaf)

leaf_solve(psi_soil = 2.0, PPFD = 900)
#>   psi_soil layers PPFD atm_vpd ca leaf_temp atm_kpa psi_stem  collar    ci
#> 1        2      1  900       2 40        25   101.3 3.595247 2.92039 10.49
#>          A         E         gc   profit ... lambda    g1_eff
#> 1 5.599511 1.142e-05 0.01921993 2.515843 ... 159884.6 0.5025448

# a drought response
leaf_solve(psi_soil = seq(0.5, 5, length.out = 20), PPFD = 900)
```

`gc` is not from a fitted conductance model — it is what falls out of maximising
profit over the hydraulic path. `lambda` is the marginal cost of water, dA/dE, at
the operating point, and `g1_eff` re-expresses the solved conductance as a Medlyn
`g1`, which is a convenient common scale for comparison.

Traits and numerical settings are separate, so a calibration loop varying traits
never has to know which of the C++ constructor's nineteen arguments are
tolerances:

```r
leaf_solve(psi_soil = 3.0, PPFD = 900,
           traits  = leaf_traits(vcmax_25 = 120, stem_b = 2.5),
           control = leaf_control(GSS_tol_abs = 1e-5))
```

For the stateful interface — which is what plant uses, and what you want if you
care about intermediate state:

```r
l <- leaf_model()                          # or leaf_model(traits, control)
set_drivers(l, psi_soil = 2.0, PPFD = 900)
l$find_root_collar_psi()

operating_point(l)   # the same one-row data.frame
l$profit_            # or reach into the object directly
l$lambda             # marginal cost of water, dA/dE
```

`Leaf()` is also exported: it is the raw C++ constructor, nineteen positional
arguments and no defaults. `leaf_model()` is that with the arguments named,
defaulted and split into traits versus tolerances, and is what you should use.

**All water potentials are positive magnitudes in MPa.** One representation
throughout, and it is asserted rather than documented — a negative `psi_soil` is
an error, not a sign convention the model quietly accepts.

⚠️ **`root_carbon_per_leaf_area` is the one argument with no good default.** A
bare-leaf user does not have a root carbon profile, and the value the R layer
supplies is a stand-in rather than a recommendation. The single-soil-potential
supply path that removes the need for it exists in the C++ and is being brought
to the R side.

See `vignette("leaf")` for the whole tour.

### As a dependency of another R package

Name it in `LinkingTo` to compile against the headers, the way `BH` is used:

```
LinkingTo: BH, odelia (>= 0.2.0), leaf (>= 0.1.0)
```

`LinkingTo` is **not** transitive in R, so you must name `BH` and `odelia`
yourself even though it is `leaf` that includes them — including the odelia
version, for the same reason. A `LinkingTo` consumer gets `<leaf.hpp>`, which is
R-free; `<leaf.h>` is the R binding layer's own umbrella and is not for you.

## Dependencies

Deliberately few. The two the *model* needs are header-only:

| | why | how |
|---|---|---|
| **odelia** (>= 0.2.0) | cubic-spline interpolator for the pre-integrated vulnerability curves, and the vendored **XAD** automatic-differentiation library | `LinkingTo` |
| **BH** (Boost) | TOMS748 root finder, incomplete gamma for the closed-form vulnerability integral | `LinkingTo` |

**Rcpp** and **R6** are needed by the R layer only. They are not in the model's
include graph and a C++ or Python consumer never sees them.

Nothing else, and **neither model dependency needs R**. The leaf model itself does not
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

Two suites, and the C++ one is the regression baseline.

```sh
make -C tests/cpp        # plain C++: no R, no test framework
```

It discovers BH and odelia through `Rscript` if R is installed, and otherwise
falls back to a sibling `odelia/` checkout and Homebrew Boost. Override with
`make BH_INC=... ODELIA_INC=...`. `ctest --test-dir build` runs the same two
programs through CMake.

At its centre is `tests/cpp/golden/operating_points.tsv`: 288 operating points
recorded at full precision and compared **bit-exactly**, which is what makes a
large refactor of this code checkable rather than hopeful. It is bit-exact on the
platform that generated it (macOS/arm64) and compared with per-field tolerances
elsewhere, because libm's `exp`/`pow` are not bit-reproducible across platforms.

```sh
R CMD check .            # the C++ suite, plus the R layer's own tests
```

`R CMD check` runs the C++ suite compiled with R's own configured compiler
against the installed headers — so a package that `LinkingTo`s this one finds out
from its own check when a header stops compiling. It also runs
`tests/testthat/`, which ties the R layer back to the same golden points. That
tie-back matters more than it looks: the C++ suite never loads the R layer, so a
mistranslation in the bindings would otherwise produce a green suite and
plausible R numbers. Those expected values are written as **C99 hex floats**, on
purpose — R's decimal parser is not correctly rounded and returns a value one ULP
off for roughly 18% of full-precision inputs, so decimals there would fail
against a model that is exactly right.

## API documentation

```sh
doxygen        # docs/html/index.html
```

Doxygen for the C++ API and roxygen for the R one. The headers' comments are the
substantive documentation here, and `tools/doxygen_filter.awk` presents them to
Doxygen without modifying a single source file.

## How this compares to other leaf models

See [COMPARISON.md](COMPARISON.md) for a feature-by-feature comparison against
`plantecophys`, `bigleaf` and `tealeaves`. The short version: those packages are
stronger on empirical stomatal models, leaf energy balance and fitting to
measured data; this one is the only one with an explicit hydraulic architecture
and a profit-maximisation solve, and the only one written to be embedded in a
larger model.

## Licence

AGPL (>= 3), inherited from plant.
