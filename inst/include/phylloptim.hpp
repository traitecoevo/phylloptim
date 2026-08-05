// -*-c++-*-
#ifndef PHYLLOPTIM_HPP_
#define PHYLLOPTIM_HPP_

// leaf: a header-only leaf gas-exchange and hydraulics model.
//
// One include gets you phylloptim::Leaf, which couples Farquhar-von Caemmerer-Berry
// photosynthesis to an explicit soil-root-stem-leaf hydraulic path and picks the
// operating point by gain-risk profit maximisation. There is nothing to link
// against.
//
// Typical use:
//
//   #include <phylloptim.hpp>
//
//   phylloptim::Leaf l;                       // default Eucalyptus saligna traits
//   l.setup_transpiration(100);         // build the vulnerability splines
//   l.setup_root_vulnerability(100);
//   l.set_physiology(root_carbon_per_leaf_area, PPFD, psi_soil, soil_depth,
//                    leaf_specific_conductance_max, atm_vpd, ca,
//                    leaf_temp, atm_o2_kpa, atm_kpa);
//   l.find_root_collar_psi();           // solve; results in l.opt_psi_stem_ etc.

#include <phylloptim/closed_form.hpp>
#include <phylloptim/constants.hpp>
#include <phylloptim/leaf_model.hpp>
#include <phylloptim/optimize.hpp>
#include <phylloptim/quadrature.hpp>
#include <phylloptim/roots.hpp>
#include <phylloptim/single_potential.hpp>
#include <phylloptim/uniroot.hpp>
#include <phylloptim/util.hpp>
#include <phylloptim/vulnerability.hpp>

#endif
