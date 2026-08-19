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

// The molar mass of water -- 18.015 g/mol, from the standard atomic weights -- and
// the ONE place it is written down. Both conversion directions derive from it, so
// they are reciprocal by construction; do not reintroduce a second literal (#51).
//
// Both names below are kept because plant `using`-declares both.
inline constexpr double molar_mass_h2o = 0.018015; // kg mol^-1
// kg mol^-1, for converting a molar water flux back to kg.
inline constexpr double kg_per_mol_h2o = molar_mass_h2o;
// mol H2O kg^-1 -- the reciprocal, now genuinely so.
inline constexpr double kg_to_mol_h2o = 1.0 / molar_mass_h2o;
// mol mol ^-1 / (umol mol ^-1)
inline constexpr double umol_to_mol = 1e-6;
// Pa kPa^-1
inline constexpr double kPa_to_Pa = 1000.0;

// universal gas constant J mol^-1 K^-1.
//
// ⚠️ 8.314 until the #51 audit, which is a truncation rather than a convention:
// since the 2019 SI redefinition R = N_A * k_B is **exact**, at
// 8.314462618153240 J mol^-1 K^-1. So there is no version of this where 8.314 is
// the right number, and the only question was whether the 5.56e-05 relative
// correction was worth a results change.
//
// It is, because the amplification is what matters and it is not 5.56e-05. This
// constant appears ONLY in arrh_curve and peak_arrh_curve, always as `Ea/(R*T)`
// with Ea/RT of order 24 at the defaults, so a relative change eps in R moves the
// exponent by ~24*eps and the rate by ~1.3e-03 -- an order above this package's
// 1e-04 "real difference" threshold. Measured on the golden grid: worst 2.9e-03.
//
// ⚠️ AND IT MOVES NOTHING AT 25 C, exactly. `arrh_curve` has (leaf_temp - 25) in
// its numerator, so it returns ref_value with exp(0) whatever R is; and
// peak_arrh_curve's arg2 and arg3 are the same expression at the same temperature,
// so their ratio is exactly 1. Every reference value in this model is DEFINED at
// 25 C, which is why the golden grid carries a second temperature -- and this
// change is the cleanest demonstration of why: the 25 C block comes out
// byte-identical and the whole diff is at 40 C.
inline constexpr double gas_constant = 8.314462618153240;

//convert deg C to deg K
inline constexpr double C_to_K = 273.15;

//H20:CO2 stomatal diffusion ratio
inline constexpr double H2O_CO2_stom_diff_ratio = 1.67;

// rho_w * g, as a head in MPa per metre of depth.
//
// ⚠️ NOT CHANGED BY THE #51 AUDIT, AND FLAGGED RATHER THAN FIXED. 9.8e-3 is 6.78e-04
// below `1000 * 9.80665 / 1e6 = 9.80665e-03`, which is above this package's 1e-04
// "real difference" threshold — so it is not a rounding, and it does move results
// through the per-layer head `gravity_head * z_soil_mid_`.
//
// It is left alone because unlike the gas constant there is no forced answer: the
// value depends on a water-density convention nobody has stated. rho_w = 1000
// (nominal) gives 9.80665e-03; at 4 C, 999.97 gives 9.80636e-03; at 25 C, 997.05
// gives 9.77772e-03 — a 0.3% spread, four times the discrepancy being corrected.
// Picking one is a modelling decision, and 9.8e-3 may well be a deliberate round
// number rather than a slip.
//
// So: whoever needs this decides and says which density, in an issue of its own.
// Recorded here so the next audit does not have to re-derive the spread.
inline constexpr double gravity_head = 9.8e-3; // MPa / m

// --- Penman-Monteith leaf energy balance (minimal core; #523) -----------------
// See notes/penman-monteith/. These back Leaf::leaf_temp_from_E and the es/Delta
// helpers. Only used on the (default-off) use_energy_balance_ path.
// latent heat of vaporisation of water, J kg^-1.
//
// ⚠️ THE VALUE AND ITS OWN COMMENT DISAGREE, and the #51 audit left this for #28
// rather than picking one. lambda(T) ~= 2.501e6 - 2370*T (J/kg, T in C), so at the
// 25 C this claimed to be "fixed at" the value is 2.442e6 — 2.45e6 is 3.3e-03 above
// that, and corresponds to about 21.5 C. The comment said 25 C; the number does not.
//
// Not corrected here for two reasons. It reaches only the Penman-Monteith path,
// which is default-off, so the golden grid cannot see it either way — and on that
// path lambda is arguably the wrong SHAPE rather than the wrong value, because it is
// a function of the leaf temperature the balance is solving for. Fixing the constant
// and then making it temperature-dependent would move PM results twice. #28 is where
// that belongs; the comment is corrected so it no longer states a temperature the
// value does not have.
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
