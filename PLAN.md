# Plan: tightening phylloptim after 0.27.0

Written for: agents picking up work on this package. Each item names what to change, where, what evidence says it is worth doing, and what "done" is. The claims were verified against commit e898a0f (branch `collar-dry-pin-continuity`, 0.27.0) on 2026-11-02; where an item says "measure first", the number quoted came from an audit agent's prototype and has to be reproduced before it is acted on.

## Rules that apply to every item

- **A change to any model goes through `inst/model-registry.yml`**, and `tests/testthat/test-registry.R` has to agree with the code afterwards. The registry is a record, not a check: nothing fails when the two drift, so the agent is the check.
- **No slowdown unless it is essential, and every claim of cost is measured interleaved.** Run `tests/cpp/bench_solve` on master and on the branch alternately, three times each, on an idle machine, and quote both. For plant-facing changes also time a plant stand (the seasonal stand in plant's `test-scm.R` is the near-end-heavy one). Hazard 20 in [`.claude/hazards.md`](.claude/hazards.md) has what a wrong bench reads.
- **Every trait has gradient infrastructure in R.** A new differentiable parameter needs its `par_*` slot before `par_PPFD`, a row in `.solve_price_fields` and `.gradient_owned_pars` if it is model-owned, a write in `.gradient_setter`, and a row in `.level_defaults` if it has a default. The `fitting` vignette's "Fitting a different collection of parameters" says why.
- **Bit-identical first.** A refactor that claims to change no number is asserted on all three goldens (`make -C tests/cpp build && ./test_golden && ./test_primitives`) before anything else is believed. Reductions (`TF24_floor` at `lambda_o = 0` is `TF24`) are asserted bit for bit, not to a tolerance.
- **Comments carry contracts and hazards, not history.** State the reason, never the issue number. The nineteen hazards in `.claude/hazards.md` are the list; a new one is appended as 20.
- Open an issue when you pick an item up, and close it from the PR. Measurements, alternatives rejected and suite counts go in the PR's first comment, not its body ([plant-meta's commit-message rule](https://github.com/traitecoevo/plant-meta/blob/main/governance/commit-messages.md)).

## Order, and why

1. Two latent bugs (section A). Cheap, and each is a wrong number waiting for a caller.
2. The CI gaps (section B). Cheap, and they stop the rot the audit found in sections D and E from recurring.
3. The measured speed win (section C1), then the two smaller ones.
4. The documentation sweep (section D). Large but mechanical; do it after B so the new checks hold it.
5. The infrastructure generalisation (section E). The largest piece. It needs its own design pass before code; do not start it from this document alone.

## A. Latent bugs — required

**A1. `evaluate_psi_stem` leaves `lambda_emergent_` stale.** `inst/include/phylloptim/leaf_model.hpp`, `Leaf::evaluate_psi_stem` (the stem-route evaluation at a prescribed stem potential, ~line 6358) writes every output of the point except `lambda_emergent_`, so a caller reading the shadow price after it gets the previous solve's. Every other exit of a solve writes it (`lambda_collar_for<K>()` or `lambda_for<K>()`); hazard 8 says why every path must. Done when the function writes it, `test_leaf.cpp`'s stale-state sweep covers the stem route's prescribed evaluation, and the goldens are bit-identical.

**A2. `leaf_batch()` never checks the prices.** `R/gradient-batch.R`, `leaf_batch()` (line 85) takes the model-owned prices as arguments and seats them without `.check_prices()`, which `leaf_solve()` (`R/leaf-model.R`, line 1083) calls on the same arguments. A batch built with a price its model does not own, or a missing one, fails later or not at all. Done when the batch calls the same check and `test-gradient-batch.R` asserts the same refusal `leaf_solve()` gives.

**A3. The C++ fixtures' `stem_P50`.** `tests/cpp/test_transpose.cpp:398`, `test_leaf.cpp:1727` and `:4863` hard-code `stem_P50 = 3.4287` where `leaf_traits()` defaults to 3.4 and plant's trait-built default is 2.889. Nothing is wrong numerically, but a fixture that is "the default" in its comment and is not in its value misleads the next reader. Done when each fixture either reads the default or says in one line why it does not.

## B. Contracts the CI does not hold — required

**B1. `tests/testthat/test-readme.R` never runs on CI.** It skips with `skip_if_not(file.exists(readme))`, and under `R CMD check` the README is not reachable from the test directory. The `cpp-tests.yml` workflow has the source tree; run the file from there with `testthat::test_file()` on the checkout, or copy the README into `tests/testthat/` at check time. Done when a deliberately broken README chunk fails CI.

**B2. The registry assertions that read the source tree skip on CI** for the same reason: `test-registry.R` lines 134 and 148 (`skip_if_not(file.exists(path))`, `skip_if_not(... "index.xml")`). Same fix as B1. Done when a registry entry pointing at a symbol that does not exist fails CI.

**B3. No cost gate.** `cpp-tests.yml` compiles `bench_solve` and `bench_gradient` but asserts nothing; `tools/cost-baseline.tsv` carries the history but stops at an old API (`tail -1` it). Decide one of: (a) a gate that runs `bench_solve` on CI and fails past a generous bound (shared runners are noisy; `test-cost.R` chose `skip_on_ci()` for that reason), or (b) a required line in the PR template that quotes the interleaved bench for any change under `inst/include`, with the baseline file appended per release. (b) is the recommendation: it costs nothing on CI and the rule above already demands the measurement. Done when the PR template carries the line and the baseline is current.

**B4. The PR template does not mention the registry.** `.github/PULL_REQUEST_TEMPLATE.md` says nothing about `inst/model-registry.yml`. Add one line: a change to any model's equations, parameters, defaults or operating-point classification updates the registry in the same PR. Done when it does.

**B5. The derivative surface's FD-versus-IFT check covers a minority of parameters.** `test_leaf.cpp`'s referees compare the implicit-function rows against a difference of the solve for the parameters each fixture happens to perturb. Extend one referee to every `par_*` at every operating-point kind on a small grid, with the per-kind tolerances the existing tests measured. Done when the referee lists every parameter by name and a planted wrong row in any one of them fails.

## C. Speed — measure first

**C1. The collar objective walks the root uptake three times per candidate collar.** `Leaf::dprofit_at_collar_psi` and the functions it calls (`E_from_Soil_to_Root_Collar`, the uptake derivative) each rebuild the per-layer uptake for the same collar. The audit's prototype shared one walk and measured 15% faster on `bench_solve`, bit-exact on all three goldens. Reproduce that: share the walk through the existing per-solve cache, assert bit-identity, bench interleaved. Done when the bench moves by what the PR claims and the goldens do not.

**C2. The vulnerability-curve rows are rebuilt per parameter row in the gradient.** `gradient.hpp` differentiates the pre-integrated tables' rows for each of the 36 parameters, and the table's own trait rows (`cumulative_vulnerability_complement_derivatives_at`) do not depend on which parameter is being driven. Memoise them per operating point. The audit did not measure this one; measure `bench_gradient` first, and drop the item if it is under 5%.

**C3. The ci solve's Newton slope.** `psi_stem_to_ci` lifts a slope it already has one call later; the audit estimated ~1%. Do this only alongside C1, on the same bench run; it is not worth a PR of its own.

**C4. `leaf_batch()` crosses to C++ once per supply kind.** `R/gradient-batch.R` branches on `supply_kind` (lines 152, 336) in R and reaches C++ separately for the single and layered paths. `test-cost.R` bounds the per-row cost, so this is only worth doing if that bound is what stops a fit; measure a 64-row batch first.

## D. Documentation — required, mechanical

**D1. Issue numbers in code comments.** `leaf_model.hpp` carries 110 `#NNN` references, and `R/`, `gradient.hpp` and the tests more. Each is replaced by the reason it stood for, or deleted where the surrounding comment already gives it. Do it file by file, one PR per file, with no code change in the same PR, so the diff is reviewable. The rule is in `.claude/CLAUDE.md` ("Writing in this guide").

**D2. History narration.** Comments that say what the code used to do ("the probe that sat here is gone", "was X before #NNN") are rewritten as what the code does and what would go wrong otherwise, or cut. Same PRs as D1.

**D3. `COMPARISON.md` is stale on its headline claims.** It calls the package `leaf`, cites "PLAN 7b", says Prentice14 is "planned" (LeastCost ships), lists two runnable optimality schemes (nine ship), says the energy balance is "minimal", and quotes a cost per solve. Rewrite the table against `vignettes/the-models.Rmd` and the current bench, or delete the file and point at the vignette. Deleting is acceptable if nothing links to it but `.claude/CLAUDE.md`.

**D4. Stale counts.** The reference grid is 672 points, not 576: `README.md:601`, `.claude/measurements.md:7`, `vignettes/the-models.Rmd:27` and `:202` say 576 (`.claude/CLAUDE.md:54` is right). `vignettes/cpp-interface.Rmd:169` quotes `static_assert(g::n_theta == 19)`; `n_theta` is 35 (`gradient.hpp:103`), and the same vignette's line 151 says so. `README.md`'s Status section is at 0.10.0; `NEWS.md` is at 0.27.0. `README.md:557`, `:570` and `:583` put odelia's floor at 0.7.0; `DESCRIPTION` says `>= 0.10.0`. Fix each from the code, not from another document, and add a line to `test-readme.R` (after B1) that greps the two counts out of the goldens.

**D5. Thermal damage has no worked example outside `the-models.Rmd`.** The `fitting` and `phylloptim` vignettes do not mention it. One short section in `the-models.Rmd` showing the gate on and off at one driver is enough; it is optional.

## E. One roof for models and gradients — a design pass, then code

These are the audit's findings on where a model or a gradient is bespoke rather than generic. They are listed so the design pass starts from facts; none is to be done piecemeal, because each touches the dispatch every other item depends on. The design has to say what a model IS (the enum arm, the cost/profit/lambda trio, the registry entry, the R price field) and make everything else derive from that, so that adding a tenth curve touches the places `TF24_floor` touched and no others (hazards.md, "One solver: the long form").

- **E1.** `find_root_collar_psi()` solves TF24 whatever `set_model()` seated (`leaf_model.hpp:2318` documents it as the TF24 shorthand). plant calls `find_root_collar_psi_for<K>` for that reason. Decide whether the shorthand survives; if it does, it dispatches on the seated model.
- **E2.** `ProfitMax` is name-branched at 13 sites (`K == CostCurve::ProfitMax`). Each is a property (needs normalisers, reference potential, no closed form) that belongs in one `if constexpr` table beside the others.
- **E3.** The derivative surface (`supply_draw_at`, `operating_collar`, `billed_outputs_at`, `dry_end_collar_at`, ...) is instantiated for `TF24` and `TF24_floor` only, and `cost_deriv<K>` is `double`-only (`leaf_model.hpp:2584`), so no other curve can be differentiated through a run. Decide which curves plant will ever differentiate and either template the cost derivative on the scalar or state the restriction in the registry.
- **E4.** The closed form (`closed_form.hpp`) is unreachable from `leaf_solve()`; it exists for the bench and plant's seeding. Either bind it or say in the registry that it is internal.
- **E5.** Medlyn is present in C++ (62 mentions) but is not a `CostCurve`; R reaches it only through `g1_eff`. `COMPARISON.md` says "present, not dispatched". Decide: a curve, or deleted.
- **E6.** The damage enum overlaps the operating-point kinds in meaning (shutdown versus damage-limited). One enumeration, or two with disjoint names.
- **E7.** The model list is copied into `.claude/CLAUDE.md`, `the-models.Rmd`, the README and `COMPARISON.md`, with counts that disagree (seven, eight, nine). After the design pass the registry is the one list and the documents quote it (`test-registry.R` already reads the registry; extend it to grep the documents' counts).
- **E8.** Fifteen exits of the collar solve write `operating_point_kind_` by hand; eight `profit_psi_stem_*` bodies share their shape; the temperature finite-difference block appears twice; the closed-form Newton is duplicated. These are the places a tenth curve would have to touch today. Collapse them only after E1–E3 have decided the dispatch; collapsing first bakes today's dispatch in.

## Done

- 0.27.0: the dry-end decision is made on the end solved to collar precision, on the collar's own distance; `boundary-dry` and `interior` meet continuously; plant's sweep agrees with a tangent of the run to 7e-8. Measurements on [phylloptim#205](https://github.com/traitecoevo/phylloptim/pull/205).
