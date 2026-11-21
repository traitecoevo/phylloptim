// -*-c++-*-
#ifndef PHYLLOPTIM_CLAMP_SITES_HPP_
#define PHYLLOPTIM_CLAMP_SITES_HPP_

#include <cstddef>
#include <memory>
#include <vector>

namespace phylloptim {

// Where this model replaces a value by a bound, counted. Two bounds: the stomatal
// corner E_s >= 0, and the leaf-to-air vapour deficit held at vpd_leaf_min; every
// other bound is a refusal at the entry point or a smooth form. The corner is
// counted, not refused, because it is part of the model: the value function is
// C1 there and the argmax has one-sided derivatives. The deficit's floor is a
// kept exception (below it the leaf condenses), and the marginal profit has a
// kink where it binds, so a consumer differentiating across it needs to know.
// The floor has two sites: every read that floors, which is the search's
// incidence, and the seat a solve leaves the leaf at, which is the state's.
//
// A consumer keeping its own tally folds these in rather than sharing an enum: a
// header-only model cannot depend on the consumer that includes it.
// ⚠️ Append only: a consumer mirrors this list by position.
enum clamp_site {
  // The stomatal flux read as zero where the hydraulic flux falls short of the
  // residual flux: at the operating point only by rounding at the wet bound, but
  // routinely in the derivative probes just wet of it.
  CLAMP_STOMATAL_FLUX = 0,
  // The vapour deficit held at vpd_leaf_min, where a leaf cooled by the energy
  // balance reaches the dew point: every read that floors it, on every candidate
  // a search tries.
  CLAMP_VPD_LEAF_FLOOR,
  // The deficit a solve leaves the leaf with on that floor: once as each solve,
  // replay or prescribed point returns, whatever its kind, and once as the
  // derivative surface places the held coordinates at that point, at any
  // scalar. A restore after a probe counts none. The surface places once per
  // operating_collar or outputs call that is handed no HeldCoords, so a caller
  // sharing one placement between them counts one, a caller sharing none two.
  CLAMP_VPD_LEAF_SEATED,
  CLAMP_SITE_COUNT
};

inline const char* clamp_site_name(int site) {
  switch (site) {
  case CLAMP_STOMATAL_FLUX:          return "stomatal_flux";
  case CLAMP_VPD_LEAF_FLOOR:         return "vpd_leaf_floor";
  case CLAMP_VPD_LEAF_SEATED:        return "vpd_leaf_seated";
  }
  return "unknown";
}

// The tally, behind a shared pointer so a copy still reports: plant copies its
// strategy per cohort and discards it, which would take a plain member's counts
// with it.
struct clamp_counter {
  std::shared_ptr<std::vector<std::size_t>> counts =
    std::make_shared<std::vector<std::size_t>>(CLAMP_SITE_COUNT, 0);

  void note(int site) const { ++(*counts)[site]; }
  std::size_t at(int site) const { return (*counts)[site]; }
  void clear() const { counts->assign(CLAMP_SITE_COUNT, 0); }
};

}  // namespace phylloptim

#endif
