// Timing harness for a TRAIT GRADIENT, with no R in the way.
//
//   make -C tests/cpp bench_gradient && ./bench_gradient [reps]
//   ./bench_gradient --count <arm> <reps>
//
// The gradient is one forward pass of the derivative surface carrying a lane
// per slot of theta (gradient.hpp), so its cost is flat in the number of
// parameters asked for. The arms put it beside the alternative a caller has,
// a central difference of the whole solve, at one interior point: TF24 on the
// collar route over one layer, the drivers of plant's tests/testthat/test-leaf.r.
// Every arm seats theta and solves first, so "less none" is the gradient alone.
//
// `--count <arm> <reps>` runs one arm alone, for the instructions retired that
// `/usr/bin/time -l` (macOS) or `perf stat` reports: a count repeats where
// wall-clock drifts.
//
// This harness does NOT check the gradient: test-gradient-referee.R and the
// surface's referees in test_leaf do. The printed checksum folds in every value
// the arms compute, so a change that moved the arithmetic shows here too.
//
// tools/bench_history.sh reads the `vcmax_25` line: the difference of the
// solve, then the gradient, each in microseconds less the seat and solve.

#include <phylloptim.hpp>

#include "root_network.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

namespace pl = phylloptim;
namespace gr = phylloptim::gradient;
using clock_type = std::chrono::steady_clock;

const double kTheta = 0.000157, kKs = 1.0, kH = 5.0, kAreaLeaf = 0.05;

gr::Drivers drivers() {
  std::vector<double> root{1.0 / kAreaLeaf}, psi_soil{2.0}, depth{1.0};
  gr::Drivers d;
  d.root_network = fixture::root_network(root, depth);
  d.PPFD = 900.0;
  d.psi_soil = psi_soil;
  d.soil_depth = depth;
  d.atm_vpd = 2.0;
  d.ca = 40.0;
  d.leaf_temp = 25.0;
  d.atm_o2_kpa = 21.0;
  d.atm_kpa = 101.3;
  return d;
}

// theta at the leaf's defaults, varied per rep in TF24_beta2 by an
// inconsequential amount so nothing can be hoisted out of a timed loop.
void fill_theta(const pl::Leaf& l, double* theta, long r) {
  const pl::leaf_pars<double> p = l.passive_pars();
  for (int j = 0; j < gr::n_theta; ++j) theta[j] = p[std::size_t(j)];
  theta[pl::par_kmax] = kKs * kTheta / kH;
  theta[pl::par_g_min] = 0.0;
  theta[pl::par_TF24_beta2] *= 1.0 + double(r) * 1e-12;
}

enum class Kind { None, Gradient, Difference };

struct Arm {
  const char* label;
  Kind kind;
  int npars;  // the first npars slots, every one TF24 reads up to 17
};

const Arm arms[] = {
    {"none (seat + solve)", Kind::None, 0},
    {"gradient, 1 par", Kind::Gradient, 1},
    {"gradient, 4 pars", Kind::Gradient, 4},
    {"gradient, 17 pars", Kind::Gradient, 17},
    {"difference of the solve, 1 par", Kind::Difference, 1},
    {"difference of the solve, 4 pars", Kind::Difference, 4},
    {"difference of the solve, 17 pars", Kind::Difference, 17},
};
constexpr int n_arms = int(sizeof(arms) / sizeof(arms[0]));

gr::Settings settings() {
  gr::Settings s;
  s.model = {pl::Leaf::CostCurve::TF24, gr::Route::Collar};
  return s;
}

double one(pl::Leaf& l, const Arm& arm, long r, const gr::Drivers& d) {
  double theta[gr::n_theta];
  fill_theta(l, theta, r);
  int pars[gr::n_theta];
  for (int i = 0; i < arm.npars; ++i) pars[i] = i;
  const gr::Settings s = settings();
  if (arm.kind == Kind::Gradient) {
    gr::Result res;
    gr::at(l, theta, d, false, pars, std::size_t(arm.npars), s, res);
    double sum = 0.0;
    for (double v : res.grad) sum += v;
    return sum;
  }
  gr::apply(l, theta, d, false);
  gr::route_seat(l, s.model);
  l.find_root_collar_psi();
  double sum = l.profit_;
  if (arm.kind == Kind::Difference) {
    double th[gr::n_theta], up[gr::n_outputs], dn[gr::n_outputs];
    for (int k = 0; k < arm.npars; ++k) {
      // A slot at zero (g_min, CMax_b) may not go negative: a forward step.
      const bool forward = theta[k] == 0.0;
      const double h = std::max(std::abs(theta[k]), 1.0) * 1e-6;
      for (int side = 0; side < 2; ++side) {
        std::copy(theta, theta + gr::n_theta, th);
        th[k] += side == 0 ? h : (forward ? 0.0 : -h);
        gr::apply(l, th, d, false);
        l.find_root_collar_psi();
        gr::outputs(l, side == 0 ? up : dn);
      }
      for (int j = 0; j < gr::n_outputs; ++j) {
        sum += (up[j] - dn[j]) / (forward ? h : 2.0 * h);
      }
    }
  }
  return sum;
}

double us_per(clock_type::time_point a, clock_type::time_point b, long n) {
  return std::chrono::duration<double, std::micro>(b - a).count() / double(n);
}

const Arm* find(const char* name) {
  for (const Arm& a : arms) {
    if (std::string(a.label).rfind(name, 0) == 0) return &a;
  }
  return nullptr;
}

}  // namespace

int main(int argc, char** argv) {
  pl::Leaf l = fixture::multilayer_leaf();
  const gr::Drivers d = drivers();
  double sink = 0.0;

  if (argc == 4 && std::strcmp(argv[1], "--count") == 0) {
    const Arm* arm = find(argv[2]);
    if (arm == nullptr) {
      fprintf(stderr, "no arm starts with '%s'\n", argv[2]);
      return 2;
    }
    const long n = std::atol(argv[3]);
    for (long r = 0; r < n; ++r) sink += one(l, *arm, r, d);
    printf("%s x %ld: checksum %.6f\n", arm->label, n, sink);
    return 0;
  }

  const long reps = argc > 1 ? std::atol(argv[1]) : 3000;
  const long per_round = std::max(1L, reps / 3);
  // Min of three interleaved rounds: the least-contaminated estimate, and
  // sequential A-then-B is not good enough at these effect sizes.
  double best[n_arms];
  for (int a = 0; a < n_arms; ++a) best[a] = 1e300;
  for (int round = 0; round < 3; ++round) {
    for (int a = 0; a < n_arms; ++a) {
      const auto t0 = clock_type::now();
      for (long r = 0; r < per_round; ++r) sink += one(l, arms[a], r, d);
      best[a] = std::min(best[a], us_per(t0, clock_type::now(), per_round));
    }
  }
  printf("%-36s %10s %12s\n", "arm (TF24, collar, 1 layer)", "us", "less none");
  for (int a = 0; a < n_arms; ++a) {
    printf("%-36s %10.2f %12.2f\n", arms[a].label, best[a], best[a] - best[0]);
  }
  // The line bench_history.sh reads: difference of the solve, then the gradient.
  printf("\nvcmax_25 %10.3f %10.3f\n", best[4] - best[0], best[1] - best[0]);

  // Per observation through the batch, which is what a calibration pays: the
  // seat, the solve and the pass, for 20 observations at once.
  {
    const std::size_t n_obs = 20;
    const std::vector<gr::Drivers> obs(n_obs, d);
    const gr::Settings s = settings();
    printf("\n%-24s %10s\n", "batch, 20 observations", "us/obs");
    for (int npars : {1, 4, 17}) {
      int pars[gr::n_theta];
      for (int i = 0; i < npars; ++i) pars[i] = i;
      double best_obs = 1e300;
      const long batch_reps = std::max(1L, per_round / 20);
      for (int round = 0; round < 3; ++round) {
        const auto t0 = clock_type::now();
        for (long r = 0; r < batch_reps; ++r) {
          double theta[gr::n_theta];
          fill_theta(l, theta, r);
          const std::vector<gr::Result> g =
              gr::batch(l, theta, 1, obs, false, pars, std::size_t(npars), s);
          sink += g[0].grad[0];
        }
        best_obs = std::min(best_obs,
                            us_per(t0, clock_type::now(), long(n_obs) * batch_reps));
      }
      printf("%-24s %10.2f\n", (std::to_string(npars) + " pars").c_str(), best_obs);
    }
  }
  printf("\nchecksum %.6f\n", sink);
  return 0;
}
