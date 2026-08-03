// -*-c++-*-
#ifndef LEAF_HPP_
#define LEAF_HPP_

// leaf: a header-only leaf gas-exchange and hydraulics model.
//
// One include gets you leaf::Leaf, which couples Farquhar-von Caemmerer-Berry
// photosynthesis to an explicit soil-root-stem-leaf hydraulic path and picks the
// operating point by gain-risk profit maximisation. There is nothing to link
// against.
//
// Typical use:
//
//   #include <leaf.hpp>
//
//   leaf::Leaf l;                       // default Eucalyptus saligna traits
//   l.setup_transpiration(100);         // build the vulnerability splines
//   l.setup_root_vulnerability(100);
//   l.set_physiology(area_leaf, mass_root_prop, rho, a_bio, PPFD,
//                    psi_soil, soil_depth, leaf_specific_conductance_max,
//                    atm_vpd, ca, sapwood_volume_per_leaf_area,
//                    leaf_temp, atm_o2_kpa, atm_kpa);
//   l.find_root_collar_psi();           // solve; results in l.opt_psi_stem_ etc.

#include <leaf/closed_form.hpp>
#include <leaf/constants.hpp>
#include <leaf/leaf_model.hpp>
#include <leaf/optimize.hpp>
#include <leaf/quadrature.hpp>
#include <leaf/roots.hpp>
#include <leaf/uniroot.hpp>
#include <leaf/util.hpp>
#include <leaf/vulnerability.hpp>

#endif
