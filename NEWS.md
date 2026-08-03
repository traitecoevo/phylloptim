# leaf 0.1.0

**The first version bump since the package was created, and it exists so downstream can
pin.** `0.0.1` has been the version through four merged PRs that changed the API and the
results, so a `leaf (>= 0.0.1)` requirement was satisfied by the original extraction —
which `plant` cannot compile against at all. That is the same trap odelia 0.2.1 was cut
for: a version that names several different header sets is not something a consumer can
depend on.

A minor bump rather than a patch, because the C++ API broke in ways a consumer sees at
compile time, and the numbers moved.

## Breaking API changes

- **`set_physiology` takes 10 arguments, not 14** (#15). `area_leaf`, `rho`, `a_bio` and
  `sapwood_volume_per_leaf_area` were dead stores — assigned, never read — and are gone,
  along with the four members of the same names. `mass_root_prop` became
  `root_carbon_per_leaf_area`, and it is the old value **divided by `area_leaf`**: the
  model is purely intensive, and uptake is exactly homogeneous in that ratio. Passing the
  old absolute carbon compiles and runs, and gives a root system too weak by a factor of
  `1/area_leaf`.
- **`root_collar_psi_` is now `opt_root_psi_`, and it is a positive magnitude** (#25).
  Renamed rather than reused deliberately: keeping the name with a flipped sign is the one
  outcome where an old analysis reads the wrong value in silence, and a rename gives a
  binding error instead.
- **Every ψ in the package is a positive magnitude in MPa** (#25).
  `psi_soil_inverted_`, `psi_soil_inverted_vec_` and `supply_psi_soil_inverted()` are
  deleted. `E_from_Soil_to_Root_Collar`, `find_root_psi`, `find_psi_stem_from_psi_root`,
  `dE_from_soil_dpsi_collar` and `transpiration_to_psi_stem` keep their signatures but
  take magnitudes, and `find_root_psi`'s bracket ends swap (wettest layer first). They
  validate the soil vector and stop on a negative entry, so a pre-#25 caller fails loudly
  rather than returning a wrong number.
- **`dE_from_soil_dpsi_collar` returns a positive conductance** (#25), so callers that
  negated it to recover one must stop.
- **Renames** (#15): `b`/`c` → `stem_b`/`stem_c` (there are two Weibull curves and the
  unmarked pair was the source of a real error), `g1_TF24` → `cost_scale_TF24`.
- **`umol_per_mol_to_Pa` is no longer a namespace-scope constant** (#15). It was
  `0.1013` = 101.3 kPa in disguise; it is now the member `Leaf::umol_per_mol_to_Pa_`,
  derived per call from `atm_kpa_`, so the model is self-consistent away from sea level.

## Behaviour changes

Each landed with its own measured blast radius against the golden file; see the PRs.

- **Four exits no longer leave stale state** (#15, #26). `Leaf` is reused for every
  individual in a patch, so any output a branch declined to write became the previous
  plant's value. Fixed: the three `set_shutdown_state` call sites, the `assim_max_ < 0`
  exit, `soil_consumption_` cleared with `.assign` not `.resize`, and
  `dprofit_droot_collar_psi` (which segfaulted without a prior solve and returned a
  spurious gradient in the reversed-gradient state).
- **The collar bracket is clamped to `root_psi_crit`** (#24, plant #584). The clamp
  compared a magnitude against a signed potential and could never bind. The window is
  empty at this package's defaults and 1.2 MPa wide at plant's.

## Downstream

`plant` requires `leaf (>= 0.1.0)`; anything below it will not compile against plant's
`feature/consume-leaf-package`.
