// -*-c++-*-
#ifndef _LEAF_H_
#define _LEAF_H_

// ============================================================================
// THIS IS THE R-FACING UMBRELLA HEADER, AND IT PULLS IN Rcpp.
//
// If you are embedding the model in C++ -- which is what `LinkingTo: leaf`
// means, and what plant does -- you want <leaf.hpp>, not this file. That one is
// plain C++ with no R anywhere in its include graph, and keeping it that way is
// issue #11's guarantee.
//
// Only src/RcppR6.cpp includes this header. Nothing in inst/include/leaf/
// includes it, and nothing there may: the dependency runs one way, downward
// from the R layer into the model, and the moment it runs the other way plant
// starts inheriting Rcpp through a header it compiles into its own translation
// units.
//
// The .h / .hpp split carries that whole distinction, which is more weight than
// a file extension should have to bear. It is not a style choice -- RcppR6
// hardwires both `#include <leaf.h>` in the code it generates and
// `<leaf/RcppR6_pre.hpp>` in the umbrella it expects to find, so the name is
// fixed and the model's umbrella had to be the one that moved. Hence this
// banner. See PLAN item 6a for the decision and for what actually enforces it:
// .github/workflows/cpp-tests.yml builds the C++ suite with no R installed, so
// a header that reaches for Rcpp turns three jobs red before review.
// ============================================================================

// The model itself, R-free.
#include <leaf.hpp>

// Forward declarations for the wrap/as specialisations. Must come before
// anything includes Rcpp.h.
#include <leaf/RcppR6_pre.hpp>

// Anything below here may include Rcpp.h.

// The definitions. Safe to be the last line, and is.
#include <leaf/RcppR6_post.hpp>

#endif
