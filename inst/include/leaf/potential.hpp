// -*-c++-*-
#ifndef LEAF_POTENTIAL_HPP_
#define LEAF_POTENTIAL_HPP_

namespace leaf {

// ---------------------------------------------------------------------------
// WATER POTENTIAL: TWO CONVENTIONS, NOW IN THE TYPE  (issue #8)
// ---------------------------------------------------------------------------
// This model carries water potential in two conventions, each natural to its
// domain, and until now the only thing keeping them apart was a comment block
// plus the `_inverted` suffix:
//
//   Potential -- SIGNED, <= 0, MPa. The soil -> root-collar transport works in
//                real signed gradients (psi_soil - P_collar - gravity*z), so
//                this is the convention the supply physics wants.
//   Suction   -- a POSITIVE MAGNITUDE, -MPa. The stem/leaf demand side, the
//                vulnerability splines and every published trait value
//                (psi_crit, root_psi_crit) are quoted this way.
//
// The two differ by a minus sign, which is exactly why mixing them is invisible:
// every wrong answer is the right magnitude. Three incidents so far, all real:
//
//   * plant #584 (and it is LIVE in this package, see prepare_collar_solve):
//     `std::max(-root_crit, -root_psi_crit)` compares a magnitude against a
//     signed potential, so the intended clamp to the root critical potential can
//     never bind. Reproduced here at psi_soil = 5.90 with psi_crit = 5.91988: the
//     bracket runs to 5.906974, past root_psi_crit = 5.870283.
//   * the multi-layer lambda came out NEGATIVE while PLAN 8 was being written,
//     because dE_from_soil_dpsi_collar differentiates with respect to the signed
//     potential and the identity wanted a positive conductance.
//   * psi_soil_ / psi_soil_inverted_ have to be kept straight by suffix alone.
//
// WHAT THESE TYPES DO AND DO NOT COVER
//
// They cover the *supply* side: soil -> root collar, where both conventions meet
// and every flip lives. They deliberately stop at the demand side (transpiration,
// psi_stem_to_ci, hydraulic_cost_TF, proportion_of_conductivity, profit_*), which
// has exactly one convention -- positive magnitudes -- and no recorded incident.
// Typing it too would put `.value` on the whole photosynthesis core to guard
// against a confusion that cannot arise there.
//
// They also do NOT cover the lambda incident above, and it is worth being clear
// about that: the wrong-signed quantity there was a *derivative*, dE_up/dpsi,
// which is neither a Potential nor a Suction. Catching that needs a conductance
// type on the supply contract's derivative method; it is not in this pair.
//
// THE R BOUNDARY STAYS `double`. plant's RcppR6 bindings list psi_soil_,
// root_collar_psi_, opt_psi_stem_ and psi_stem as `access: field` with type
// double, and generate Rcpp::as<double> / Rcpp::wrap against them by name --
// so retyping those members would break plant's generated glue, and teaching
// Rcpp about these types would reintroduce the R dependency hazard 9 exists to
// keep out. Every plant-bound field and method therefore keeps its `double`
// signature, and the typed layer sits behind thin adapters. That is also the
// answer issue #5 needs: the R interface is doubles, the types are internal.
//
// ZERO COST BY CONSTRUCTION. Each is one scalar in a trivially-copyable struct,
// so it is passed in a register exactly as the bare scalar was, and every
// conversion is a negation the compiler already emitted. Asserted below, and
// measured: see the PR for the interleaved bench.
//
// TEMPLATED ON THE SCALAR so issue #4 (template Leaf on its scalar type, to
// delete the hand-maintained AD replicas) does not have to undo this. The
// aliases at the bottom are what the model actually spells.
// ---------------------------------------------------------------------------

template <class T> struct SuctionT;

// Signed water potential, <= 0, MPa.
template <class T>
struct PotentialT {
  T value{};

  PotentialT() = default;
  constexpr explicit PotentialT(T v) : value(v) {}

  // The ONE sanctioned way to reach the other convention. Every minus sign that
  // used to be an unmarked leading `-` is now one of these.
  constexpr SuctionT<T> magnitude() const;
};

// Water potential as a positive magnitude, -MPa.
template <class T>
struct SuctionT {
  T value{};

  SuctionT() = default;
  constexpr explicit SuctionT(T v) : value(v) {}

  constexpr PotentialT<T> signed_potential() const {
    return PotentialT<T>(-value);
  }
};

template <class T>
constexpr SuctionT<T> PotentialT<T>::magnitude() const {
  return SuctionT<T>(-value);
}

// --- comparison -------------------------------------------------------------
// Same-convention only, which is the whole point: `Potential < Suction` does not
// compile, and that is what makes the prepare_collar_solve clamp visible. Only
// operator< is needed by std::min/std::max; the rest are here so call sites read
// as they did before.
#define LEAF_POTENTIAL_COMPARE(op)                                             \
  template <class T>                                                           \
  constexpr bool operator op(PotentialT<T> a, PotentialT<T> b) {                \
    return a.value op b.value;                                                  \
  }                                                                             \
  template <class T>                                                           \
  constexpr bool operator op(SuctionT<T> a, SuctionT<T> b) {                     \
    return a.value op b.value;                                                  \
  }
LEAF_POTENTIAL_COMPARE(<)
LEAF_POTENTIAL_COMPARE(>)
LEAF_POTENTIAL_COMPARE(<=)
LEAF_POTENTIAL_COMPARE(>=)
LEAF_POTENTIAL_COMPARE(==)
LEAF_POTENTIAL_COMPARE(!=)
#undef LEAF_POTENTIAL_COMPARE

// --- differences ------------------------------------------------------------
// A difference of two potentials is a plain scalar in MPa, and deliberately so:
// it is a gradient or an interval width, and unlike a potential it carries no
// sign convention to get wrong. Returning a bare scalar is what lets the physics
// in uptake_impl stay the plain arithmetic it already was.
template <class T>
constexpr T operator-(PotentialT<T> a, PotentialT<T> b) {
  return a.value - b.value;
}
template <class T>
constexpr T operator-(SuctionT<T> a, SuctionT<T> b) {
  return a.value - b.value;
}

// Midpoint of two suctions. Spelled out as a named function rather than allowing
// `a + b`, because the sum of two potentials is not a potential and nothing
// should be able to write one by accident.
template <class T>
constexpr SuctionT<T> midpoint(SuctionT<T> a, SuctionT<T> b) {
  return SuctionT<T>(0.5 * (a.value + b.value));
}

// --- what the model spells --------------------------------------------------
using Potential = PotentialT<double>;
using Suction = SuctionT<double>;

static_assert(sizeof(Potential) == sizeof(double),
              "Potential must be exactly a double: it is passed in registers on "
              "the ~10^3-evaluations-per-solve supply path.");
static_assert(sizeof(Suction) == sizeof(double), "likewise Suction");

}  // namespace leaf

#endif
