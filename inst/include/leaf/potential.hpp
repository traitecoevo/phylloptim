// -*-c++-*-
#ifndef LEAF_POTENTIAL_HPP_
#define LEAF_POTENTIAL_HPP_

namespace leaf {

// Water potential in two representations, one type each, so the compiler keeps
// them apart instead of a comment block (issue #8).
//
// WHY THERE ARE TWO OF THESE AT ALL. There is ONE physical quantity, psi, negative
// in a transpiring plant. But this model evaluates two functions of it whose
// domains are incompatible:
//
//   TRANSPORT is linear in psi and needs the sign.
//       E_i = (psi_soil_i - P_collar - rho*g*z_i) / r_i
//     The sign is physical: E_i < 0 means layer i GAINS water (hydraulic
//     redistribution). Rewritten in magnitudes you would have to remember to
//     reverse the subtraction while the gravity term does not reverse.
//
//   VULNERABILITY is a power law and CANNOT take the sign.
//       f(psi) = exp(-(psi/b)^c),   b = 3.898245, c = 2.680147
//     c is not an integer, and x^c for real non-integer c is undefined for x < 0.
//     At psi = -2.0:  exp(-(2.0/b)^c) = 0.846046, exp(-(-0.513)^c) = nan.
//
// The solve computes the COMPOSITION of the two -- transport supplies the potential
// at which vulnerability is evaluated -- so a conversion is mathematically forced.
// No naming and no refactor removes it. These types only make each crossing
// explicit, and stop the two being confused.
//
// HENCE THE NAMES, and why they are deliberately asymmetric. These are not two
// conventions of equal standing: `Psi` is the quantity, `AbsPsi` is |psi|. The
// second is derived from the first by nothing more interesting than an absolute
// value, and the name should say so. An earlier draft called them `Potential` and
// `Suction` -- both real terms, but the pairing let one member claim the physics,
// and it misled its own author into writing "psi_leaf is not a potential", which is
// false. If you rename these, keep the property that the second name reads as a
// function of the first.
//
// NEVER STORE AN AbsPsi, AND NEVER PUT ONE IN A SIGNATURE. It is a transient at the
// point of a `pow` or a spline lookup, nothing more. Everything the model holds or
// passes is a `Psi`. E_column used to take its psi_leaf as a magnitude, because
// that is the form transpiration() wants; that was accurate and wrong-headed, since
// it put both representations in one argument list and invited the reading that
// psi_leaf is a different kind of quantity from psi_soil.
//
// LONGER TERM, THE SECOND TYPE SHOULD NOT EXIST. The reason it does is not psi vs
// |psi| -- it is that this model stores psi negative and the Weibull scale b
// POSITIVE, a mixed convention inside a single formula. Store b in the same
// convention and the problem evaporates: (-2.0)/(-3.898) = +0.513, positive by
// construction and exactly equal bitwise to the ratio of the magnitudes. Then there
// is one type, no conversion, and plant #584 is unwriteable not because the mixed
// comparison fails to compile but because there is no second representation to mix
// in. What remains after that is the spline knot grids (built increasing over
// magnitudes) and the R boundary. See PLAN 10a for the sequencing.
//
// The two differ by a minus sign, which is why mixing them is invisible: every
// wrong answer has the right magnitude. Three incidents so far, all real:
//
//   * plant #584 (and it is LIVE in this package, see prepare_collar_solve):
//     `std::max(-root_crit, -root_psi_crit)` compares a magnitude against a
//     signed potential, so the intended clamp to the root critical potential can
//     never bind. Reproduced here at psi_soil = 5.90 with psi_crit = 5.91988: the
//     bracket runs to 5.906974, past root_psi_crit = 5.870283.
//   * the multi-layer lambda came out NEGATIVE while PLAN 8 was being written,
//     because dE_from_soil_dpsi_collar differentiates with respect to the signed
//     form and the identity wanted a positive conductance.
//   * psi_soil_ / psi_soil_inverted_ had to be kept straight by suffix alone.
//
// WHAT THESE TYPES DO AND DO NOT COVER
//
// They cover the *supply* side: soil -> root collar, where both representations
// meet and every flip lives. They deliberately stop at the demand side
// (transpiration, psi_stem_to_ci, hydraulic_cost_TF, proportion_of_conductivity,
// profit_*), which uses one representation throughout and has no recorded
// incident. Typing it too would cost about twelve more `double` forwarders -- that
// many are plant-bound -- plus ten `.value` unwraps where a psi meets pow, exp or
// a spline, and it would need the templated form through the AD replicas and
// closed_form.hpp, so it collides with issue #4 rather than helping it.
//
// The rule for a signature that spans both, since supply-versus-demand does not
// decide it: type a parameter whenever both representations are in play at that
// point. E_column holds one of each, so both are typed.
//
// They also do NOT cover the lambda incident above, and it is worth being clear
// about that: the wrong-signed quantity there was a *derivative*, dE_up/dpsi,
// which is neither of these. Catching that needs a conductance type on the supply
// contract's derivative method; it is not in this pair. See issue #3.
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
// measured: the bench_solve binaries built with and without this header are
// byte-identical.
//
// TEMPLATED ON THE SCALAR so issue #4 (template Leaf on its scalar type, to
// delete the hand-maintained AD replicas) does not have to undo this. Tested over
// XAD's forward-mode active type, derivatives included. The aliases at the bottom
// are what the model actually spells.

template <class T> struct AbsPsiT;

// A water potential in its SIGNED form: <= 0, MPa.
template <class T>
struct PsiT {
  T value{};

  PsiT() = default;
  constexpr explicit PsiT(T v) : value(v) {}

  // The ONE sanctioned way to reach the other convention. Every minus sign that
  // used to be an unmarked leading `-` is now one of these.
  constexpr AbsPsiT<T> abs() const;
};

// The same quantity as a MAGNITUDE: >= 0, -MPa.
template <class T>
struct AbsPsiT {
  T value{};

  AbsPsiT() = default;
  constexpr explicit AbsPsiT(T v) : value(v) {}

  constexpr PsiT<T> as_psi() const {
    return PsiT<T>(-value);
  }
};

template <class T>
constexpr AbsPsiT<T> PsiT<T>::abs() const {
  return AbsPsiT<T>(-value);
}

// --- comparison -------------------------------------------------------------
// Same representation only, which is the whole point: `Psi < AbsPsi`
// does not compile, and that is what makes the prepare_collar_solve clamp
// visible. Only
// `operator<` is needed by `std::min` and `std::max`; the rest are here so call
// sites read as they did before.
#define LEAF_POTENTIAL_COMPARE(op)                                             \
  template <class T>                                                           \
  constexpr bool operator op(PsiT<T> a, PsiT<T> b) {                \
    return a.value op b.value;                                                  \
  }                                                                             \
  template <class T>                                                           \
  constexpr bool operator op(AbsPsiT<T> a, AbsPsiT<T> b) {                     \
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
// A difference of two psi values in the SAME representation is a plain scalar in
// MPa, and deliberately so: it is a gradient or an interval width, and unlike a
// potential it carries no sign convention to get wrong. Returning a bare scalar is what lets the physics
// in uptake_impl stay the plain arithmetic it already was.
template <class T>
constexpr T operator-(PsiT<T> a, PsiT<T> b) {
  return a.value - b.value;
}
template <class T>
constexpr T operator-(AbsPsiT<T> a, AbsPsiT<T> b) {
  return a.value - b.value;
}

// Midpoint of two magnitudes. A named function rather than allowing `a + b`,
// because the sum of two potentials is not a potential and nothing should be able
// to write one by accident.
template <class T>
constexpr AbsPsiT<T> midpoint(AbsPsiT<T> a, AbsPsiT<T> b) {
  return AbsPsiT<T>(0.5 * (a.value + b.value));
}

// --- what the model spells --------------------------------------------------
using Psi = PsiT<double>;
using AbsPsi = AbsPsiT<double>;

static_assert(sizeof(Psi) == sizeof(double),
              "Psi must be exactly a double: it is passed in registers on "
              "the ~10^3-evaluations-per-solve supply path.");
static_assert(sizeof(AbsPsi) == sizeof(double), "likewise AbsPsi");

}  // namespace leaf

#endif
