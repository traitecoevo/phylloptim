# Run the C++ suite under `R CMD check` (issue #12, PLAN item 15).
#
# This package ships headers and nothing else -- no R code, no compiled code --
# so without this file `R CMD check` validates essentially nothing. The point is
# not to duplicate the GitHub Actions workflow: it is that a `LinkingTo: leaf`
# consumer's own `R CMD check`, on their machine and their toolchain, tells them
# when a header stops compiling. That is the failure mode that matters for a
# header-only package, and it is invisible to every other check R runs here.
#
# Deliberately uses R's *configured* compiler (`R CMD config CXX20`) rather than
# whatever `c++` happens to be on PATH, because that is the compiler a consumer
# will build against.

cpp_dir <- if (dir.exists("cpp")) "cpp" else file.path("tests", "cpp")

if (!dir.exists(cpp_dir)) {
  stop("cannot find the C++ suite; looked for 'cpp' and 'tests/cpp' from ",
       getwd())
}

# --- Is there a toolchain to test with? -------------------------------------
#
# A binary-only R installation with no Rtools / no command-line tools has no
# compiler, and there is nothing this file can meaningfully do there. Skip
# loudly rather than fail: a missing toolchain is not a broken header.

bail <- function(...) {
  message("SKIP: ", ..., "\n  The C++ suite was not run.")
  quit(save = "no", status = 0)
}

make <- Sys.which("make")
if (!nzchar(make)) {
  bail("no `make` on PATH.")
}

cxx <- system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX20"),
               stdout = TRUE, stderr = FALSE)
cxxstd <- system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX20STD"),
                  stdout = TRUE, stderr = FALSE)
cxx <- paste(cxx, collapse = " ")
cxxstd <- paste(cxxstd, collapse = " ")

if (!nzchar(trimws(cxx))) {
  bail("R has no C++20 compiler configured (`R CMD config CXX20` is empty).")
}

# --- Where the two header dependencies live ---------------------------------
#
# Both are LinkingTo, so they are installed whenever this package is checked.

deps <- vapply(c("BH", "odelia"),
               function(p) system.file("include", package = p),
               character(1))
missing <- names(deps)[!nzchar(deps)]
if (length(missing)) {
  bail("LinkingTo dependency not installed: ", paste(missing, collapse = ", "),
       ". Both BH and odelia are needed to compile the headers.")
}
bh <- deps[["BH"]]
odelia <- deps[["odelia"]]

# Compile against the INSTALLED headers, not ../../inst/include. Under
# `R CMD check` the two are different directories, and the installed one is what
# a `LinkingTo: leaf` consumer actually gets -- including whether the install
# shipped every header the others include.
leaf_inc <- system.file("include", package = "leaf")
if (!nzchar(leaf_inc)) {
  bail("this package's own headers are not installed.")
}

# --- Bit-exact where that is achievable, tolerant everywhere else -----------
#
# The golden file was generated on macOS/arm64 and is bit-exact only there;
# libm's exp/pow are not bit-reproducible across platforms and never were. This
# mirrors the same decision in .github/workflows/cpp-tests.yml -- see the long
# comment above main() in tests/cpp/test_golden.cpp for why it is per-field and
# why it is not a loosening.

sysinfo <- Sys.info()
on_generating_platform <-
  identical(unname(sysinfo[["sysname"]]), "Darwin") &&
  identical(unname(sysinfo[["machine"]]), "arm64")
golden_args <- if (on_generating_platform) "" else "--cross-platform"

message("Building the C++ suite with:")
message("  CXX        ", cxx, " ", cxxstd)
message("  leaf       ", leaf_inc)
message("  BH         ", bh)
message("  odelia     ", odelia)
message("  golden     ", if (nzchar(golden_args)) golden_args else "bit-exact")

status <- system2(
  make,
  c("-C", shQuote(cpp_dir),
    paste0("CXX=", shQuote(cxx)),
    paste0("CXXSTD=", shQuote(cxxstd)),
    paste0("LEAF_INC=", shQuote(leaf_inc)),
    paste0("BH_INC=", shQuote(bh)),
    paste0("ODELIA_INC=", shQuote(odelia)),
    paste0("GOLDEN_ARGS=", shQuote(golden_args)))
)

if (status != 0) {
  stop("the C++ suite failed (make exited ", status, ").\n",
       "  If only the golden comparison failed, read the worst-relative-\n",
       "  difference line it printed first: ~1e-15 is cross-platform libm\n",
       "  and not a regression, ~1e-4 is a real change in behaviour.")
}
