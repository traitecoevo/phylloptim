##' A fingerprint of the model's numerical behaviour
##'
##' A short string that changes when this package's numbers change, and does not
##' change when anything else about it does. Intended as a cache key: a consumer
##' that stores computed results can depend on this one value instead of
##' reimplementing "has phylloptim changed?" for itself.
##'
##' @section Why the package version is not enough:
##' It should be, and from 0.3.0 the discipline is that a results-moving PR moves
##' the minor version. But that is a promise about future PRs, and the history it
##' has to sit on top of is `0.1.0` covering a shut-down fix, an argmax relocation
##' worth 1.5e-03 and four gradient PRs, and `0.2.1` covering the dark-respiration
##' reallocation in #41 and five more merges. A consumer needs something that was
##' derived from the numbers rather than remembered alongside them.
##'
##' @section What it is derived from:
##' The two recorded baselines, together:
##'
##' - `tests/cpp/golden/operating_points.tsv` — 576 solved operating points
##' - `tests/testthat/gradient_golden.tsv` — 20 recorded trait gradients
##'
##' Those files already *are* this package's definition of "the numbers", they are
##' committed, and they are regenerated deliberately rather than incidentally. So
##' the fingerprint inherits the existing discipline instead of needing new
##' enforcement — and `test-fingerprint.R` is what holds it to that, by recomputing
##' the digest from the files on disk and failing if this constant disagrees. A PR
##' that regenerates a golden file and forgets the fingerprint does not pass.
##'
##' Both files are hashed by *content lines* with the trailing carriage return
##' stripped, not by raw bytes, so a CRLF checkout cannot change the answer. Their
##' *names* are deliberately not in the digest, so a rename cannot either.
##'
##' @section ⚠️ What it does NOT tell you:
##' - **Not a source hash.** Renames, comments, documentation, new *non-numerical*
##'   API and performance work all leave it alone. That is the intent — the
##'   calibration study's own workaround hashed `inst/include`, `R/` and `src/`,
##'   and #47's rename would have moved that hash without moving a single number.
##' - **Not a promise that equal fingerprints mean equal output on your machine.**
##'   The golden files are bit-exact only on macOS/arm64; elsewhere the same code
##'   agrees to ~1.4e-04 at the argmax. The fingerprint identifies the *code's*
##'   behaviour, not your platform's realisation of it.
##' - **Not complete coverage.** It sees what the golden grids reach. Notably they
##'   construct a fresh `Leaf` per point, so a stale-state bug of the kind hazard 8
##'   describes can be fixed or introduced without moving either file (all three of
##'   the fixes in #15's class were golden-bit-identical). A change in a code path
##'   no grid point visits will not move this.
##'
##' @return A single string: the digest of the two baselines' digests, truncated to
##'   12 hex characters. Stable across platforms and R versions, since it is a
##'   committed constant rather than a computation.
##' @seealso `NEWS.md`, which is the human-readable form of the same information
##' @examples
##' leaf_behaviour_fingerprint()
##' @export
leaf_behaviour_fingerprint <- function() {
  # ⚠️ GENERATED. Regenerate with tools/fingerprint.R after regenerating a golden
  # file, and never by hand -- test-fingerprint.R recomputes it and fails on a
  # mismatch, which is the whole mechanism.
  "3e36d9b94fd0"
}

# The digest itself, in one place so the exported constant above, the regeneration
# tool and the test that ties them together cannot disagree about how it is
# computed.
#
# Hashes each file's CONTENT LINES rather than its bytes, with the trailing carriage
# return stripped, so a CRLF checkout on Windows gives the same answer as an LF one
# on macOS. ⚠️ `readLines` does NOT do that stripping for you on Unix -- it keeps the
# `\r` at end of line, and the first version of this function relied on it not to.
# `test-fingerprint.R` has the case.
#
# ⚠️ THE FILE NAMES ARE NOT IN THE DIGEST, deliberately. Only the contents are, in
# name order. #58's whole complaint about hashing the source tree is that "#47's
# rename would have silently changed the hash without changing behaviour" -- so
# putting `basename()` in here would reproduce the defect the fingerprint exists to
# avoid. Name order is used to make the result independent of the order the caller
# passed them, which needs the names but does not need to record them.
#
# `tools::md5sum` is file-only, hence the tempfile. It is base R, which is the point
# -- a fingerprint a consumer cannot compute without installing `digest` first is
# one more dependency than this is worth.
golden_fingerprint_of <- function(paths) {
  stopifnot(length(paths) > 0L, all(file.exists(paths)))
  digest_text <- function(txt) {
    f <- tempfile()
    on.exit(unlink(f), add = TRUE)
    writeLines(txt, f)
    unname(tools::md5sum(f))
  }
  content_digest <- function(p) {
    digest_text(sub("\r$", "", readLines(p, warn = FALSE)))
  }
  parts <- vapply(paths[order(basename(paths))], content_digest, character(1))
  substr(digest_text(unname(parts)), 1L, 12L)
}

# The two files, resolved from a directory that is somewhere inside the package
# tree. Not usable from an INSTALLED package -- `tests/` is shipped in the tarball
# but `R CMD INSTALL` does not put it in the library -- which is exactly why the
# exported function above returns a constant instead of computing this at runtime.
golden_fingerprint_paths <- function(root = ".") {
  file.path(root, c("tests/cpp/golden/operating_points.tsv",
                    "tests/testthat/gradient_golden.tsv"))
}
