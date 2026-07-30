// -*-c++-*-
#ifndef LEAF_UTIL_HPP_
#define LEAF_UTIL_HPP_

// Minimal, R-free utilities for the leaf model.
//
// This replaces plant/util.h, which reached into Rcpp for `stop()` and into
// R_ext/Arith.h for the NA_REAL sentinel. Neither is needed: `stop()` throws a
// std::runtime_error, which Rcpp converts into an ordinary R error at the
// package boundary, and the sentinel is a quiet NaN. Keeping this header free
// of R is what lets the leaf model compile and run as plain C++ (see
// tests/cpp/), which in turn is what makes it testable and profilable without
// an R session.

#include <string>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace leaf {
namespace util {

// Throws; never returns. The attribute lets the compiler see that callers whose
// catch blocks end in util::stop() do not fall through, which avoids spurious
// -Wreturn-type warnings.
[[noreturn]] inline void stop(const std::string &msg) {
  throw std::runtime_error(msg);
}

template <typename T> std::string to_string(T x) { return std::to_string(x); }

inline bool is_finite(double x) { return std::isfinite(x); }

// Use this to be explicit when a deliberate floating-point equality test is
// being made (e.g. checking a value is exactly the one we set).
inline bool identical(double a, double b) { return a == b; }

template <typename T> T clamp(T x, T min_val, T max_val) {
  return std::max(std::min(x, max_val), min_val);
}

// Missing-value sentinel for uninitialised leaf state. Stands in for R's
// NA_REAL, which is itself a NaN with a payload; nothing in the leaf model
// inspects the payload, it only ever tests !isfinite(), so a quiet NaN is an
// exact substitute.
inline constexpr double na_value = std::numeric_limits<double>::quiet_NaN();

// Integer counterpart, for max_soil_layer before any layer is seen. Matches
// R's NA_INTEGER so the value round-trips unchanged through plant's bindings.
inline constexpr int na_value_int = -2147483648;

} // namespace util
} // namespace leaf

#endif
