// Golden-file regression test over a grid of operating points.
//
//   make -C tests/cpp golden        # regenerate golden/operating_points.tsv
//   make -C tests/cpp && ./test_golden
//
// Why this exists. The refactors in PLAN.md items 7-11 are meant to be
// behaviour-preserving, and the only way to know is to pin the behaviour first.
// PLAN.md item 1 -- cross-checking against plant's compiled build -- is the real
// validation and is still outstanding; this file is the weaker but immediately
// available version: it freezes what THIS implementation produces so that a
// refactor which changes any of it fails loudly.
//
// Comparison is bit-exact by default. Values are written with %.17g, which
// round-trips an IEEE double exactly, so a passing run means the refactor did not
// perturb a single floating-point operation. If a change is *meant* to alter
// results, regenerate deliberately and say so in the commit.
//
// Bit-exactness holds on the platform that generated the file (macOS/arm64) and
// cannot hold on any other -- `--rtol` is for those. See the comment above main().
//
// Note: a fresh Leaf is constructed for every grid point. That is not for tidiness
// -- it is required, because the shutdown-state leak (PLAN.md item 2) makes a
// reused Leaf order-dependent, which would make this file ill-defined.

#include <leaf.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace {

const char *kGoldenPath = "golden/operating_points.tsv";

struct Row {
  // inputs
  double psi_soil, ppfd, vpd;
  int layers;
  // outputs
  double psi_stem, collar, ci, assim, transpiration, gc, profit, e_up, uptake;
};

// Trait values and fixed drivers from plant's tests/testthat/test-leaf.r.
const double kTheta = 0.000157, kKs = 1.0, kH = 5.0;
const double kAreaLeaf = 0.05, kRho = 608.0, kABio = 0.0245;
const double kCa = 40.0, kO2 = 21.0, kTleaf = 25.0, kPatm = 101.3;

Row solve(double psi_soil, double ppfd, double vpd, int layers) {
  leaf::Leaf l;
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);

  // Spread the soil profile over `layers` equal 1 m layers, drying with depth so
  // that multi-layer runs are not just a repeated single layer, and split root
  // carbon evenly.
  std::vector<double> ps(layers), depth(layers), root(layers);
  for (int i = 0; i < layers; ++i) {
    ps[i] = psi_soil + 0.25 * i;
    depth[i] = 1.0 * (i + 1);
    root[i] = 1.0 / layers;
  }

  l.set_physiology(kAreaLeaf, root, kRho, kABio, ppfd, ps, depth,
                   kKs * kTheta / kH, vpd, kCa, kTheta * kH, kTleaf, kO2, kPatm);
  l.find_root_collar_psi();

  double uptake = 0.0;
  for (double s : l.soil_consumption_) {
    if (std::isfinite(s)) {
      uptake += s;
    }
  }
  return Row{psi_soil,        ppfd,   vpd,        layers,     l.opt_psi_stem_,
             l.root_collar_psi_, l.ci_, l.assim_colimited_, l.transpiration_,
             l.stom_cond_CO2_,   l.profit_, l.E_up_,        uptake};
}

std::vector<Row> run_grid() {
  const double psi_soils[] = {0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
  const double ppfds[] = {100.0, 500.0, 900.0, 1500.0};
  const double vpds[] = {0.5, 1.0, 2.0, 4.0};
  const int layer_counts[] = {1, 3, 5};

  std::vector<Row> rows;
  for (double p : psi_soils) {
    for (double q : ppfds) {
      for (double d : vpds) {
        for (int n : layer_counts) {
          rows.push_back(solve(p, q, d, n));
        }
      }
    }
  }
  return rows;
}

const char *kHeader =
    "psi_soil\tppfd\tvpd\tlayers\tpsi_stem\tcollar\tci\tassim\ttranspiration\t"
    "gc\tprofit\te_up\tuptake\n";

void write_row(FILE *f, const Row &r) {
  fprintf(f, "%.17g\t%.17g\t%.17g\t%d", r.psi_soil, r.ppfd, r.vpd, r.layers);
  for (double v : {r.psi_stem, r.collar, r.ci, r.assim, r.transpiration, r.gc,
                   r.profit, r.e_up, r.uptake}) {
    fprintf(f, "\t%.17g", v);
  }
  fprintf(f, "\n");
}

int generate() {
  const std::vector<Row> rows = run_grid();
  FILE *f = fopen(kGoldenPath, "w");
  if (f == nullptr) {
    fprintf(stderr, "cannot write %s (run from tests/cpp/)\n", kGoldenPath);
    return 1;
  }
  fputs(kHeader, f);
  for (const Row &r : rows) {
    write_row(f, r);
  }
  fclose(f);
  printf("wrote %zu operating points to %s\n", rows.size(), kGoldenPath);
  return 0;
}

// Exact equality, with NaN treated as equal to NaN -- some grid points shut down
// and legitimately produce the NA sentinel (see PLAN.md item 2).
bool same(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) {
    return true;
  }
  return a == b;
}

// Relative difference, used both to decide the tolerant comparison and to report
// how far apart the two builds actually are. Returns 0 for the NaN/NaN and
// exactly-equal cases, and infinity when one side is NaN and the other is not.
double rel_diff(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) {
    return 0.0;
  }
  if (std::isnan(a) != std::isnan(b)) {
    return std::numeric_limits<double>::infinity();
  }
  if (a == b) {
    return 0.0;
  }
  const double scale = std::max(std::abs(a), std::abs(b));
  if (scale == 0.0) {
    return 0.0;
  }
  return std::abs(a - b) / scale;
}

// rtol < 0 means bit-exact. See the comment above main() for why the tolerant
// mode exists and where it is legitimate to use it.
int compare(double rtol) {
  FILE *f = fopen(kGoldenPath, "r");
  if (f == nullptr) {
    fprintf(stderr,
            "FAIL: %s is missing. Generate it with `make golden` and commit it.\n",
            kGoldenPath);
    return 1;
  }
  char line[4096];
  if (fgets(line, sizeof line, f) == nullptr) {
    fprintf(stderr, "FAIL: %s is empty\n", kGoldenPath);
    fclose(f);
    return 1;
  }

  const std::vector<Row> rows = run_grid();
  int failures = 0;
  int inexact = 0;
  double worst_rel = 0.0;
  double worst_got = 0.0, worst_want = 0.0;
  const char *worst_desc = "";
  Row worst_row{};
  size_t i = 0;
  for (; i < rows.size(); ++i) {
    if (fgets(line, sizeof line, f) == nullptr) {
      fprintf(stderr, "FAIL: golden file has only %zu rows, grid has %zu\n", i,
              rows.size());
      ++failures;
      break;
    }
    Row g{};
    // 4 inputs then 9 outputs
    const int n = sscanf(
        line, "%lg\t%lg\t%lg\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg",
        &g.psi_soil, &g.ppfd, &g.vpd, &g.layers, &g.psi_stem, &g.collar, &g.ci,
        &g.assim, &g.transpiration, &g.gc, &g.profit, &g.e_up, &g.uptake);
    if (n != 13) {
      fprintf(stderr, "FAIL: row %zu is malformed (%d fields)\n", i, n);
      ++failures;
      continue;
    }
    const Row &r = rows[i];
    struct Field {
      const char *name;
      double got, want;
    };
    const Field fields[] = {{"psi_stem", r.psi_stem, g.psi_stem},
                            {"collar", r.collar, g.collar},
                            {"ci", r.ci, g.ci},
                            {"assim", r.assim, g.assim},
                            {"transpiration", r.transpiration, g.transpiration},
                            {"gc", r.gc, g.gc},
                            {"profit", r.profit, g.profit},
                            {"e_up", r.e_up, g.e_up},
                            {"uptake", r.uptake, g.uptake}};
    for (const Field &fd : fields) {
      const double rd = rel_diff(fd.got, fd.want);
      // Track the worst deviation over the whole grid regardless of pass/fail,
      // so the summary reports a magnitude rather than leaving the reader to
      // infer one from a truncated list.
      if (rd > worst_rel) {
        worst_rel = rd;
        worst_desc = fd.name;
        worst_row = r;
        worst_got = fd.got;
        worst_want = fd.want;
      }
      if (!same(fd.got, fd.want)) {
        ++inexact;
      }
      const bool ok = rtol < 0.0 ? same(fd.got, fd.want) : (rd <= rtol);
      if (!ok) {
        if (failures < 20) {
          fprintf(stderr,
                  "FAIL psi_soil=%g ppfd=%g vpd=%g layers=%d %s: got %.17g, "
                  "want %.17g  (rel %.3g)\n",
                  r.psi_soil, r.ppfd, r.vpd, r.layers, fd.name, fd.got, fd.want,
                  rd);
        }
        ++failures;
      }
    }
  }
  if (fgets(line, sizeof line, f) != nullptr) {
    fprintf(stderr, "FAIL: golden file has more rows than the grid (%zu)\n", i);
    ++failures;
  }
  fclose(f);

  // Always report the worst deviation. On a bit-exact pass it is 0 and says so;
  // on a tolerant pass it is the number that matters, and watching it drift is
  // the point of running the tolerant mode at all.
  char worst[256] = "none (bit-identical)";
  if (worst_rel > 0.0) {
    snprintf(worst, sizeof worst,
             "%.3g  at psi_soil=%g ppfd=%g vpd=%g layers=%d %s "
             "(got %.17g, want %.17g)",
             worst_rel, worst_row.psi_soil, worst_row.ppfd, worst_row.vpd,
             worst_row.layers, worst_desc, worst_got, worst_want);
  }

  if (failures == 0) {
    if (rtol < 0.0) {
      printf("golden: %zu operating points, all bit-identical\n", rows.size());
    } else {
      printf("golden: %zu operating points within rtol=%.1g\n"
             "  %d of %zu values differ in the last bits; worst relative "
             "difference %s\n",
             rows.size(), rtol, inexact, rows.size() * 9, worst);
    }
    return 0;
  }
  fprintf(stderr,
          "\ngolden: %d mismatches over %zu operating points%s.\n"
          "  worst relative difference %s\n",
          failures, rows.size(),
          rtol < 0.0 ? "" : " beyond tolerance", worst);
  if (rtol < 0.0) {
    fprintf(stderr,
            "If this change was intended, regenerate with `make golden` and say "
            "so in the commit message.\n"
            "If you are on a different platform from the one that generated the "
            "file, do NOT regenerate -- see --rtol.\n");
  }
  return 1;
}

} // namespace

// Usage:
//   test_golden                 bit-exact compare (the default, and the one that
//                               matters on the platform the file was generated on)
//   test_golden --rtol 1e-9     tolerant compare, for a DIFFERENT platform
//   test_golden --generate      overwrite the golden file (see the warning above)
//
// Why --rtol exists. The golden file is a *drift* detector: it pins what this
// implementation produces so a refactor that perturbs it fails loudly. That job
// requires bit-exactness, and gets it, on the platform that generated the file.
//
// It cannot be bit-exact across platforms, and it was never going to be. libm's
// exp/pow are not bit-reproducible between glibc on x86-64 and Apple's libm on
// arm64, x87/FMA contraction differs, and these nested solvers amplify the
// resulting last-bit differences. Measured across that pair: ~1e-15 relative,
// 13 ULP at worst. Regenerating the file on the second platform is NOT the fix --
// it just moves the failure to the first one.
//
// So on other platforms, compare with a tolerance wide enough to absorb
// reassociation and far narrower than a real difference. The calibration is in
// PLAN.md item 1: ~1e-16 is reassociation, ~1e-4 is a real bug, and four orders
// of magnitude separate them. rtol=1e-9 sits in that gap with room on both sides.
// The tolerant mode still prints the worst relative difference it saw, so the
// number is in the log and a drift in it is visible.
int main(int argc, char **argv) {
  if (argc > 1 && std::strcmp(argv[1], "--generate") == 0) {
    return generate();
  }
  double rtol = -1.0; // negative = bit-exact
  if (argc > 2 && std::strcmp(argv[1], "--rtol") == 0) {
    rtol = std::atof(argv[2]);
    if (!(rtol > 0.0)) {
      fprintf(stderr, "--rtol needs a positive number, got '%s'\n", argv[2]);
      return 2;
    }
  } else if (argc > 1 && std::strcmp(argv[1], "--rtol") == 0) {
    fprintf(stderr, "--rtol needs a value, e.g. --rtol 1e-9\n");
    return 2;
  } else if (argc > 1) {
    fprintf(stderr, "unknown argument '%s'\n", argv[1]);
    return 2;
  }
  return compare(rtol);
}
