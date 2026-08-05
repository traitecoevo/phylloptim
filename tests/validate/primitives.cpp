// ⚠️⚠️ THIS FILE DOES NOT COMPILE, AND HAS NOT SINCE #15. Read this before
// trying to use it.
//
// The set_physiology call below passes the PRE-#15 thirteen-argument signature
// (area_leaf, root carbon, rho, a_bio, ..., sapwood volume per leaf area). #15
// deleted those four dead arguments; #33 then replaced the root carbon profile
// with a RootNetwork. So this harness is two interface changes behind, and
// nothing builds it -- neither `make -C tests/cpp`, nor CMake, nor CI.
//
// It is left in place rather than deleted because the METHOD it encodes is still
// the right one (see the reasoning below: compare the primitives, not the solve,
// because the nested solvers destroy the localisation). Reviving it means
// updating the call and deciding whether the primitive list is still current.
// Filed as #64.
//
// Dump the leaf model's PRIMITIVES at fixed inputs, for issue #13.
//
// `compare_with_plant.R` compares end-to-end solve outputs against plant and
// finds 585 of 2592 values differing at 1-2 ULP. That comparison cannot localise
// the cause, and the issue says why: these nested solvers amplify perturbations
// up to about GSS_tol_abs (1e-3), so by the time a difference reaches the output
// it has been through a golden-section search and a ci iteration and no longer
// points anywhere. The fix is to stop comparing the solve and start comparing
// the functions it calls.
//
// So this evaluates each primitive DIRECTLY, at inputs chosen here rather than
// found by a solver, and prints every value at full precision. Its counterpart
// `compare_primitives.R` does the same through plant's R bindings and diffs the
// two. Because the functions are emitted in call-tree order -- pure arithmetic
// first, then things that call it -- the FIRST row class that disagrees is the
// answer.
//
// Build (from the repository root):
//
//   c++ -std=c++20 -O2 -I inst/include \
//       -isystem $(Rscript -e 'cat(system.file("include", package="odelia"))') \
//       -isystem $(Rscript -e 'cat(system.file("include", package="BH"))') \
//       -o /tmp/leaf_primitives tests/validate/primitives.cpp
//
// It writes TSV on stdout: function, arg1, arg2, value.

#include <phylloptim.hpp>

#include <cstdio>
#include <vector>

namespace {

// Identical to tests/cpp/test_golden.cpp and to new_leaf() in
// compare_with_plant.R. Any drift between the three invalidates the comparison,
// which is why all three spell every value out rather than relying on defaults.
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0;
const double kAreaLeaf = 0.05, kRho = 608.0, kABio = 0.0245;
const double kCa = 40.0, kO2 = 21.0, kTleaf = 25.0, kPatm = 101.3;

// Every number is emitted TWICE: once as decimal, for a human reading the file,
// and once as a C99 hex float, which is what the comparison actually uses.
//
// This is not belt and braces, it is the whole point. R's decimal parser is not
// correctly rounded -- `as.numeric`, `scan` and `read.delim` all share it, and
// on a random sample it returns a double one ULP off the correctly rounded value
// for about 18% of inputs. Reading a C++-generated decimal file into R therefore
// perturbs the C++ side by 1 ULP while leaving R-computed values untouched, which
// manufactures exactly the disagreement issue #13 set out to explain. Hex floats
// are exact in both directions: verified on 4000 values, 4000 exact, against
// 3265 of 4000 for `%.17g`.
void row(const char *fn, double a1, double a2, double v) {
  std::printf("%s\t%a\t%a\t%.17g\t%a\n", fn, a1, a2, v, v);
}

} // namespace

int main() {
  phylloptim::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);

  // One representative configuration: a single 1 m layer at psi_soil = 2 MPa,
  // mid-range light and VPD. Single-layer keeps the root network out of it --
  // the multi-layer path is a separate question and would confound this one.
  std::vector<double> ps{2.0}, depth{1.0}, root{1.0};
  l.set_physiology(kAreaLeaf, root, kRho, kABio, /*PPFD*/ 900.0, ps, depth,
                   kKs * kTheta / kH, /*vpd*/ 2.0, kCa, kTheta * kH,
                   kTleaf, kO2, kPatm);

  std::printf("function\targ1\targ2\tvalue_dec\tvalue\n");

  // --- Tier 1: pure arithmetic, no state, no splines, no iteration ----------
  // If anything here disagrees, the cause is libm or FMA contraction and the
  // search is over.
  // Swept over temperature, not evaluated at kTleaf: both curves are normalised
  // to 25 C, so at 25 C they return the reference value unchanged and would
  // compare equal no matter how broken the exp() underneath was.
  const double eas[] = {58550.0, 62500.0, 36380.0};
  const double refs[] = {96.0, 157.44, 40.0};
  for (double t = 5.0; t <= 45.0 + 1e-12; t += 2.5) {
    for (double ea : eas) {
      row("arrh_curve", ea, t, l.arrh_curve(ea, refs[0], t));
    }
    for (double ref : refs) {
      row("peak_arrh_curve", ref, t,
          l.peak_arrh_curve(43900.0, ref, t, 200000.0, 640.0));
    }
  }

  // --- Tier 2: the vulnerability curve. Calls exp/pow; no spline, no solve ---
  for (double psi = 0.0; psi <= 6.0 + 1e-12; psi += 0.25) {
    row("proportion_of_conductivity", psi, 0.0,
        l.proportion_of_conductivity(psi));
  }

  // --- Tier 3: assimilation. Rational arithmetic, then a quadratic root ------
  for (double ci = 1.0; ci <= 39.0 + 1e-12; ci += 1.0) {
    row("assim_rubisco_limited", ci, 0.0, l.assim_rubisco_limited(ci));
    row("assim_electron_limited", ci, 0.0, l.assim_electron_limited(ci));
    row("assim_colimited", ci, 0.0, l.assim_colimited(ci));
  }

  // --- Tier 4: transpiration. The first function that touches a spline -------
  // This is where the interesting money is: the spline is built by odelia, and
  // odelia's interpolator is the one dependency the two builds share by header
  // rather than by copy.
  for (double psi_stem = 0.25; psi_stem <= 5.5 + 1e-12; psi_stem += 0.25) {
    row("transpiration", psi_stem, 0.0, l.transpiration(psi_stem, 0.0));
    row("stom_cond_CO2", psi_stem, 0.0, l.stom_cond_CO2(psi_stem, 0.0));
  }

  // --- Tier 5: functions containing their own iteration ----------------------
  // psi_stem_to_ci runs the ci solve; transpiration_to_psi_stem inverts the
  // spline. A disagreement that FIRST appears here, with tiers 1-4 clean, would
  // mean the primitives agree and the iteration count differs.
  for (double psi_stem = 0.25; psi_stem <= 5.5 + 1e-12; psi_stem += 0.25) {
    row("psi_stem_to_ci", psi_stem, 0.0, l.psi_stem_to_ci(psi_stem, 0.0));
  }
  // The inverse spline's domain ends at the transpiration reached at psi_crit,
  // and asking beyond it throws rather than extrapolating. Derive the ceiling
  // from the model instead of hard-coding one, so this cannot silently go out of
  // range if a trait changes.
  const double e_max = l.transpiration(l.supply_psi_crit(), 0.0);
  for (int i = 1; i < 20; ++i) {
    const double e = e_max * i / 20.0;
    row("transpiration_to_psi_stem", e, 0.0,
        l.transpiration_to_psi_stem(e, 0.0));
  }

  return 0;
}
