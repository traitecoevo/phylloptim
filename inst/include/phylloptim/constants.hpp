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

// universal gas constant J mol^-1 K^-1. Exact since the 2019 SI redefinition
// (R = N_A k_B), so it is not a value to round.
//
// ⚠️ Read only by arrh_curve and peak_arrh_curve, always as `Ea/(R*T)` with Ea/RT
// of order 24 -- so a relative change here reaches the rates ~24x amplified, and a
// change too small to look like it matters does.
inline constexpr double gas_constant = 8.314462618153240;

//convert deg C to deg K
inline constexpr double C_to_K = 273.15;

//H20:CO2 stomatal diffusion ratio
inline constexpr double H2O_CO2_stom_diff_ratio = 1.67;

// rho_w * g, as a head in MPa per metre of depth.
//
// ⚠️ 9.8e-3 is 6.8e-04 below `1000 * 9.80665 / 1e6`, above this package's 1e-04
// "real difference" threshold, and it reaches results through the per-layer head.
// Deliberately not "corrected": the right value needs a water-density convention
// stated first (nominal 1000, 4 C and 25 C span 0.3%, four times the discrepancy),
// so it is a modelling decision rather than arithmetic.
inline constexpr double gravity_head = 9.8e-3; // MPa / m

// --- Penman-Monteith leaf energy balance (minimal core; #523) -----------------
// See notes/penman-monteith/. These back Leaf::leaf_temp_from_E and the es/Delta
// helpers. Only used on the (default-off) use_energy_balance_ path.
// latent heat of vaporisation of water, J kg^-1, at about 21.5 C.
//
// ⚠️ NOT a 25 C value, which an earlier comment here claimed: lambda(T) ~= 2.501e6
// - 2370*T gives 2.442e6 at 25 C. Left as-is because on the (default-off)
// Penman-Monteith path lambda should become a function of the leaf temperature the
// balance solves for, not a better constant -- see #28.
inline constexpr double latent_heat_vap = 2.45e6;
// volumetric heat capacity of air, J m^-3 K^-1
inline constexpr double vol_heat_cap_air = 1200.0;
// PAR energy conversion: ~4.57 umol photons per J of PAR (shortwave)
inline constexpr double umol_par_per_joule = 4.57;
// shortwave absorbed ~= 2 * absorbed PAR (PAR ~= 50% of shortwave; doc 3.3)
inline constexpr double sw_abs_per_par = 2.0;
// --- Longwave radiation (#97 / #28) ------------------------------------------
//
// ⚠️ A `longwave_net_offset = -40.0` USED TO STAND HERE and the net longwave was
// that constant. It is not a constant: over Tair 20-60 C the physical term runs
// about -80 to -135 W m^-2, so a fixed -40 is out by 40-95 W m^-2 at the hot end,
// which is where the energy balance matters most.
//
// The treatment now: an ISOTHERMAL net radiation (as if the leaf sat at air
// temperature) plus a RADIATIVE CONDUCTANCE for the leaf's own departure from it.
//
//   Rn_iso = absorbed shortwave - (1 - ema)*sigma*Tair_K^4
//   g_rad  = 4*eps_leaf*sigma*Tair_K^3
//
// ⚠️ THIS IS WHAT KEEPS THE BALANCE EXPLICIT, which is the objection #97 and PLAN
// 13.2 both raise against a faithful longwave. A leaf's own emission is
// sigma*eps*Tleaf_K^4, which would make Rn a function of the temperature being
// solved for; linearising it about Tair turns that dependence into one extra term
// in the DENOMINATOR of the balance. So Tleaf is still algebraic in E, dTleaf/dE
// is still a constant, and no inner iteration is introduced. The linearisation is
// the standard radiation-conductance one (Monteith & Unsworth; the same step
// plantecophys takes) and is accurate to the extent that (Tleaf - Tair) is small
// against Tair_K -- 20 K against 300 K is a 4th-order term of order 1%.
inline constexpr double stefan_boltzmann = 5.670374419e-8;  // W m^-2 K^-4
// Longwave emissivity of a leaf. Thermal-infrared leaf emissivities are measured
// at 0.94-0.99; 0.97 is the usual working value.
inline constexpr double leaf_emissivity = 0.97;
// Clear-sky atmospheric emissivity, Brutsaert (1975) as Campbell & Norman write
// it: ema = 0.642*(ea/Tair_K)^(1/7).
//
// ⚠️ ea IS IN PASCALS HERE, and getting that wrong is silent. The same relation is
// more often quoted as `1.24*(ea/T)^(1/7)` with ea in hPa; the two agree because
// 0.642*100^(1/7) = 0.642*1.9307 = 1.2395. With ea in kPa the coefficient would
// have to be 1.723, and using 0.642 with kPa returns ~0.31 where the answer is
// ~0.82 -- a plausible-looking number that triples the longwave loss.
inline constexpr double atmos_emissivity_coef = 0.642;
inline constexpr double atmos_emissivity_exponent = 1.0 / 7.0;
// Floor on the actual vapour pressure used above, Pa. `ea = esat(Tair) - D_air`
// goes non-positive for a deficit larger than saturation, which is not a physical
// atmosphere but is reachable from a driver sweep; the emissivity would then be
// NaN through pow() of a negative base.
inline constexpr double vapour_pressure_min = 1.0;
inline constexpr double zero_celsius_in_kelvin = 273.15;
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
// Floor on the leaf-to-air vapour pressure deficit (kPa) used by Fick's law
// (PLAN 13.1). The same argument as the temperature clamp above, one step
// downstream: a leaf cooled below the dew point has a NEGATIVE deficit, and
// dividing a positive transpiration by it reports a negative stomatal
// conductance. Flooring keeps such points finite and unattractive rather than
// sign-inverted. 0.01 kPa is far below any operating deficit -- the driver
// defaults are 1.5-2 kPa -- so it binds only where the model is already outside
// what a one-way diffusion equation describes.
inline constexpr double vpd_leaf_min = 0.01;


} // namespace phylloptim

#endif
