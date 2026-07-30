// -*-c++-*-
#ifndef LEAF_TEST_SHIM_RCPPCOMMON_H
#define LEAF_TEST_SHIM_RCPPCOMMON_H

// Stand-in for Rcpp's RcppCommon.h, so the leaf test suite builds and runs with
// no R installation.
//
// This file is also a specification. The leaf model itself is R-free -- see
// leaf/util.hpp -- and the ONLY thing that still reaches for Rcpp is odelia's
// odelia/ode_util.hpp, which interpolator.hpp includes for its `util::stop` and
// its Rcpp::as/wrap declarations for odelia::util::index. What you see below is
// the entirety of that remaining coupling: a SEXP typedef, two function
// templates that are only ever declared, and stop/warning.
//
// If odelia makes those bits optional (see PLAN.md, "Drop the last R coupling"),
// this shim can be deleted and `leaf` becomes usable as plain C++ with nothing
// standing in for R at all.

#include <string>
#include <vector>
#include <stdexcept>

typedef void *SEXP;

namespace Rcpp {
template <typename T> SEXP wrap(const T &);
template <typename T> T as(SEXP);
[[noreturn]] inline void stop(const std::string &msg) {
  throw std::runtime_error(msg);
}
inline void warning(const std::string &) {}
} // namespace Rcpp

#endif
