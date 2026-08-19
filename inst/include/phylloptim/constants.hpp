// -*-c++-*-
#ifndef PHYLLOPTIM_CONSTANTS_HPP_
#define PHYLLOPTIM_CONSTANTS_HPP_

// Physical and physiological constants for the leaf model.
//
// Split out of leaf_model.h and changed from `static const double` (internal
// linkage, so one copy per translation unit) to `inline constexpr double`,
// which is the right spelling for a header-only library: one entity, usable in
// constant expressions, no per-TU duplication.
//
// Two names changed on the way out of plant, both because they were namespace
// -scope one-letter identifiers in a public header:
//   * `R` (the universal gas constant) -> `gas_constant`. `R` at plant namespace
//     scope was a live collision hazard in a project where `R` also names the
//     language, `R_` prefixes the C API, and several strategies declare local
//     `R`s of their own.
//   * `n` ("number of integration steps") was dead -- nothing read it -- and is
//     dropped. It shadowed local `n`s in every consumer that had one.

#include <cmath>

namespace phylloptim {

// double check best namespace for constants (private vs global)
// kJ mol ^-1
inline constexpr double vcmax_ha = 60000;
// kJ mol ^-1
inline constexpr double vcmax_H_d = 200000;
// kJ mol ^-1
inline constexpr double vcmax_d_S = 650;

// kJ mol ^-1
inline constexpr double jmax_ha = 30000;
// kJ mol ^-1
inline constexpr double jmax_H_d = 200000;
// kJ mol ^-1
inline constexpr double jmax_d_S = 650;

// umol ^ -1 mol ^ 1
inline constexpr double gamma_25 = 42.75;
// kJ mol ^-1
inline constexpr double gamma_ha = 37.83e3;

// umol mol ^-1
inline constexpr double kc_25 = 404.9 ;
// kJ mol ^-1
inline constexpr double kc_ha = 79.43e3;

// umol mol ^-1
inline constexpr double ko_25 = 278400 ;
// kJ mol ^-1
inline constexpr double ko_ha = 36.38e3;

// umol mol^-1 -> Pa is deliberately NOT a constant here: it is 1e-6 * P. See
// Leaf::umol_per_mol_to_Pa_, derived from atm_kpa_. The old `0.1013` was 101.3
// kPa hard-coded; don't reinstate it.

// The molar mass of water, and the ONE place it is written down (#51). Everything
// that crosses between a kg basis and a mol basis is derived from this, so the two
// directions are reciprocal by construction and a round trip is the identity to
// within one rounding.
//
// ⚠️ THIS USED TO BE TWO INDEPENDENT NUMBERS THAT DISAGREED BY 0.0277%, and the
// header said so on purpose: `kg_to_mol_h2o` was 55.4939 while `kg_per_mol_h2o` was
// 0.018015, with a comment recording that the first was "intentionally distinct
// from 1/kg_to_mol_h2o ... kept at the historical 0.018015 to preserve results".
// Two constants naming the same physical quantity in opposite directions, used in
// opposite halves of the model: the demand side converted transpiration kg -> mol
// with the first, the supply side converted uptake mol -> kg with the second.
//
// WHICH ONE WAS WRONG IS NOT A MATTER OF CONVENTION, which is what let this be
// settled rather than argued. 55.4939 is 1/0.018020, i.e. it encodes a molar mass of
// 18.0200 g/mol. The molar mass of water is 18.015 g/mol -- from the standard atomic
// weights, 2(1.008) + 15.999 -- so 0.018015 is the physical value and the forward
// constant was the odd one. Unifying therefore moves `kg_to_mol_h2o` from 55.4939 to
// 55.509298..., which is the LARGER of the two possible moves, and it is still the
// right one: the alternative preserves more digits by adopting a molar mass water
// does not have.
//
// It also happens to leave plant's `using ::phylloptim::kg_per_mol_h2o` --
// tf24_strategy.cpp's water consumption rate, the only live use of either name over
// there -- multiplying by an unchanged constant. That is a convenience, not the
// reason; plant's results still move, because the leaf's operating point does.
//
// Both old names are kept, and deliberately: plant `using`-declares both, so
// deleting either would break its build for no benefit.
inline constexpr double molar_mass_h2o = 0.018015; // kg mol^-1
// kg mol^-1, for converting a molar water flux back to kg.
inline constexpr double kg_per_mol_h2o = molar_mass_h2o;
// mol H2O kg^-1 -- the reciprocal, now genuinely so.
inline constexpr double kg_to_mol_h2o = 1.0 / molar_mass_h2o;
// mol mol ^-1 / (umol mol ^-1)
inline constexpr double umol_to_mol = 1e-6;
// Pa kPa^-1
inline constexpr double kPa_to_Pa = 1000.0;

// universal gas constant J mol^-1 K^-1
inline constexpr double gas_constant = 8.314;

//convert deg C to deg K
inline constexpr double C_to_K = 273.15;

//H20:CO2 stomatal diffusion ratio
inline constexpr double H2O_CO2_stom_diff_ratio = 1.67;

inline constexpr double gravity_head = 9.8e-3; // MPa / m

// --- Penman-Monteith leaf energy balance (minimal core; #523) -----------------
// See notes/penman-monteith/. These back Leaf::leaf_temp_from_E and the es/Delta
// helpers. Only used on the (default-off) use_energy_balance_ path.
// latent heat of vaporisation of water, J kg^-1 (fixed at 25 deg C)
inline constexpr double latent_heat_vap = 2.45e6;
// volumetric heat capacity of air, J m^-3 K^-1
inline constexpr double vol_heat_cap_air = 1200.0;
// PAR energy conversion: ~4.57 umol photons per J of PAR (shortwave)
inline constexpr double umol_par_per_joule = 4.57;
// shortwave absorbed ~= 2 * absorbed PAR (PAR ~= 50% of shortwave; doc 3.3)
inline constexpr double sw_abs_per_par = 2.0;
// fixed net longwave (cooling) offset, W m^-2 (clear-sky approximation; doc 3.3).
// A faithful treatment makes outgoing longwave depend on Tleaf (would make Rn
// implicit in Tleaf); staged for a sensitivity gate and tracked in #581.
inline constexpr double longwave_net_offset = -40.0;
// fixed aerodynamic resistance fallback, s m^-1 (doc 6/7.4; used when the wind
// model is unavailable, e.g. a bare Leaf with no wind/d set)
inline constexpr double aerodynamic_resistance_fixed = 50.0;
// leaf boundary-layer coefficient C_ra, s^0.5 m^-1, in ra = C_ra*sqrt(d/U) (doc 4.1)
inline constexpr double aerodynamic_resistance_coef = 200.0;
// Physical clamp on the energy-balance leaf temperature (deg C). The linear
// balance can return absurd temperatures at non-equilibrium operating points an
// optimiser may probe (a large transpiration driving Tleaf below absolute zero,
// making the Arrhenius block non-finite). Real leaves operate far inside this
// range; clamping keeps A(Tleaf) finite so such points get a finite (poor)
// profit and are simply rejected.
inline constexpr double leaf_temp_min = -40.0;
inline constexpr double leaf_temp_max = 70.0;


} // namespace phylloptim

#endif
