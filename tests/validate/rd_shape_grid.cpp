// The #41 respiration-shape A/B, on this package's own grid.
//
// Three arms over the 1152 operating points of tests/cpp/golden/, differing only
// in which model computes them. It exists because the golden file can only say
// "something moved": this says WHICH of the four #41 commits moved it, by holding
// the shape fixed and letting everything else vary.
//
//   arm            model                                     expectation
//   master         phylloptim at origin/master               the reference
//   branch         this branch, default (declining Q10)      moves off 25 C
//   branch_old     this branch, rd_tracks_vcmax_ = true      bit-identical to master
//
// The third arm is the point. If it is bit-identical then the whole of #41 --
// the cache key, the settable R_d_25, the shut-down exit, the R API -- is
// behaviour-neutral, and every difference in the second arm is the respiration
// shape and nothing else. Measured: 0 of 10368 cells differ. See PLAN item 41.
//
// ⚠️ ONE TREE, THREE BINARIES, as hazard 5 requires of any A/B here. The master
// arm is built from master's HEADERS against this file, not from master's own
// test_golden.cpp -- master's grid is 288 points at 25 C and would answer a
// different question.
//
//   W=$(git rev-parse --show-toplevel)
//   mkdir -p /tmp/ab && cd /tmp/ab
//   cp $W/tests/cpp/root_network.hpp .
//   git -C $W archive origin/master inst/include | tar -x --strip-components=2 \
//     -C master_inc                      # mkdir it first
//   ODELIA=$(Rscript -e 'cat(system.file("include", package="odelia"))')
//   BH=$(Rscript -e 'cat(system.file("include", package="BH"))')
//   for arm in master branch branch_old; do
//     case $arm in
//       master)     INC=master_inc;         DEF="" ;;
//       branch)     INC=$W/inst/include;    DEF="" ;;
//       branch_old) INC=$W/inst/include;    DEF="-DOLD_SHAPE" ;;
//     esac
//     c++ -std=c++20 -O2 -I. -I$INC -isystem $ODELIA -isystem $BH $DEF \
//       $W/tests/validate/rd_shape_grid.cpp -o grid_$arm
//     ./grid_$arm > out_$arm.tsv
//   done
//
// Then diff the TSVs. They are written with %.17g in the golden file's own column
// order, so `cmp` is a valid first pass and a per-field relative comparison is the
// second. ⚠️ Do NOT read them back with R's `read.delim` and compare as doubles
// without going through tests/validate/tsv_to_hex.c: R's decimal parser is not
// correctly rounded and perturbs about 18% of full-precision values, which is the
// mistake decision 1 in PLAN.md records.
//
// Not built by `make` or by CMake, and deliberately so: it needs two versions of
// the package at once, which no build in this repo can express.
#include <phylloptim.hpp>

#include "root_network.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main() {
  // The golden grid's fixed drivers and axes, from tests/cpp/test_golden.cpp.
  const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;
  const double kCa = 40.0, kO2 = 21.0, kPatm = 101.3;
  const double temps[] = {15.0, 25.0, 35.0, 40.0};
  const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
  const double ppfds[] = {100.0, 500.0, 900.0, 1500.0};
  const double vpds[] = {0.5, 1.0, 2.0, 4.0};
  const int layer_counts[] = {1, 3, 5};

  printf("psi_soil\tppfd\tvpd\tleaf_temp\tlayers\tpsi_stem\topt_root_psi\tci\t"
         "assim\ttranspiration\tgc\tprofit\te_up\tuptake\n");
  for (double t : temps) {
    for (double p : psi_soils) {
      for (double q : ppfds) {
        for (double v : vpds) {
          for (int n : layer_counts) {
            phylloptim::Leaf l;
            l.setup_transpiration(100);
            l.setup_root_vulnerability(100);
#ifdef OLD_SHAPE
            // ⚠️ The whole reason this arm exists. Absent on master, where the
            // shape is unconditional -- so `-DOLD_SHAPE` must never be combined
            // with master's headers, and the compiler will say so.
            l.rd_tracks_vcmax_ = true;
#endif
            std::vector<double> ps(n), depth(n), root(n);
            for (int i = 0; i < n; ++i) {
              ps[i] = p + 0.25 * i;
              depth[i] = 1.0 * (i + 1);
              root[i] = 1.0 / n / kAreaLeaf;
            }
            l.set_physiology(fixture::root_network(root, depth), q, ps, depth,
                             kKs * kTheta / kH, v, kCa, t, kO2, kPatm);
            l.find_root_collar_psi();
            double uptake = 0.0;
            for (double s : l.soil_consumption_) {
              if (std::isfinite(s)) {
                uptake += s;
              }
            }
            printf("%.17g\t%.17g\t%.17g\t%.17g\t%d", p, q, v, t, n);
            for (double x : {l.opt_psi_stem_, l.opt_root_psi_, l.ci_,
                             l.assim_colimited_, l.transpiration_,
                             l.stom_cond_CO2_, l.profit_, l.E_up_, uptake}) {
              printf("\t%.17g", x);
            }
            printf("\n");
          }
        }
      }
    }
  }
  return 0;
}
