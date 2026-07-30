# Next steps

Ordered by what blocks what. Items 1-3 are the ones that have to happen before
this package can be trusted or used; everything after is improvement.

---

## 1. Validate against plant — bit-identical or explain why not

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

## 2. Fix the shutdown-state leak

**A real defect, inherited from plant unchanged, and it affects plant today.**

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

Branch `feature/consume-leaf-package` in plant does this: delete
`inst/include/plant/leaf_model.h` and `src/leaf_model.cpp`, add
`LinkingTo: leaf`, and provide a compatibility shim so that `plant::Leaf` and the
handful of leaf constants plant reads (`kg_per_mol_h2o`) keep resolving. Nothing
in plant's own sources should need to change beyond that shim.

Remaining work on that branch: build it, run the test suite, and decide whether
to keep the `plant::Leaf` alias permanently or migrate plant's ~12,000 lines of
generated RcppR6/RcppExports code to `leaf::Leaf`. The alias is much cheaper and
costs nothing at runtime.

## 4. Drop the last R coupling

The leaf model is R-free. odelia is not: `odelia/ode_util.hpp` includes
`RcppCommon.h` for `Rcpp::stop`/`Rcpp::warning` and for `Rcpp::as`/`wrap`
declarations on `odelia::util::index`, and `odelia/interpolator.hpp` includes it
for `util::stop`. That is the only reason this package cannot claim to be plain
C++ outright.

`tests/cpp/shim/RcppCommon.h` is a 15-line stand-in that satisfies all of it,
which both keeps the test suite R-free and specifies exactly how small the
problem is. The upstream fix is for odelia to guard its Rcpp bits — e.g. an
`ODELIA_NO_R` switch, or moving `stop`/`warning` to a throw the way
`leaf/util.hpp` does. File against odelia.

## 5. Settle where XAD comes from

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
which is the obvious thing to want. Make each one selectable and give all of them
the same contract: multi-layer soil, the same `set_physiology` inputs, the same
outputs, driven through one entry point. Then add **Prentice et al. (2014)**
least-cost alongside them.

Design notes:

- Prefer compile-time dispatch (a policy template parameter, or a
  `Leaf<Formulation>`) over a runtime `enum` and a `switch`. This is a
  header-only library whose whole point is that the hot loop inlines; a virtual
  call or a branch per profit evaluation would land inside a golden-section search
  running ~10³ inner evaluations. Composes naturally with item 8.
- Medlyn is the awkward one, because it is not a profit model — it prescribes
  `gs` and has no hydraulic cost. It needs either a shim that maps its `gs` onto a
  `psi_stem` through the supply curve, or an honest admission that it sits at a
  different interface. Decide which before writing the dispatch, not after.
- Keep `optimise_psi_stem_*`'s single-layer forms as the unit-test entry points
  even after the multi-layer versions exist; they are much easier to reason about.

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
conductance. Template `Leaf` on it so the calls still inline; a `std::function`
or a virtual base would put an indirect call inside a loop running ~10³ inner
evaluations, which is the thing to avoid.

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

Do 7b before 7a: the supply interface is the thing the alternative cost functions
need to plug into.

## 8. Template `Leaf` on its scalar type

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

There is a second, larger payoff, which is item 9: templating on the scalar type
is what turns "we have AD" into "we can calibrate". See below.

## 9. Demonstrate calibration — and then consider inversion

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

Do this after item 8. Attempting it before will produce a vignette that
finite-differences trait gradients, which is precisely the thing being argued
against.

**Then the bigger question: inversion.** The clearest capability gap against
`plantecophys` is that `fitaci` inverts A-Ci curves for Vcmax/Jmax/Rd and this
package cannot invert anything. Calibrating *hydraulic* traits — P50-type
vulnerability parameters, root conductances — from gas-exchange and water-potential
data would be a genuinely new capability rather than a reimplementation of
`fitaci`, and the AD story is what would make it tractable. Scope it properly
before committing; it is a project, not a task, and it needs someone to decide
what data it is meant to consume.

## 10. Energy balance — the full cut

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

## 11. Naming, home, publication

- **Package name.** `leaf` is clear inside this family and too generic outside
  it. Decide before anything is published: `leafhydro`, `hydroleaf` and
  `leafoptim` all say more.
- **Repository home.** DESCRIPTION points at `traitecoevo/leaf`; create it, or
  change the URL.
- **Paper.** The natural framing is the one gap in the existing R landscape:
  every leaf gas-exchange package assumes a stomatal conductance model, and this
  one derives stomatal behaviour from hydraulics instead. See
  [COMPARISON.md](COMPARISON.md).

## 12. Housekeeping

- **CI.** GitHub Actions matrix building `tests/cpp` on gcc and clang, Linux and
  macOS. The suite needs no R, so this is fast and catches the portability
  problems (`long double`, FMA contraction) that have bitten plant before —
  see the arm64 note in `hydraulic_cost_Sperry`.
- **Run the C++ suite under `R CMD check`** so `LinkingTo` consumers get told
  when a header breaks.
- **Doxygen or roxygen for the C++ API.** The header comments are unusually
  good — they explain why, not what — and are worth rendering.
