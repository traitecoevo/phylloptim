#!/usr/bin/env Rscript
# Rewrite the constant `leaf_behaviour_fingerprint()` returns, from the two golden
# files as they stand on disk.
#
#   Rscript tools/fingerprint.R
#
# Run it after regenerating a golden file, in the same commit. `test-fingerprint.R`
# recomputes the digest and fails if the constant disagrees, so forgetting this is
# caught rather than shipped -- which is the only reason the constant is trustworthy.
#
# ⚠️ Unlike tools/gradient_golden.R this needs no installed build and touches no
# model code: it reads two committed text files and edits one string literal. So it
# is safe to run from a working tree, and it does not care which phylloptim is on
# the library path.
#
# ⚠️ It edits R/fingerprint.R in place. The replacement is anchored on the quoted
# string inside the function body, so reformatting that function will break this --
# if you do, fix the pattern here rather than pasting a digest by hand.

root <- normalizePath(".")
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("run this from the package root", call. = FALSE)
}

source(file.path(root, "R/fingerprint.R"))  # for the two internal helpers

paths <- golden_fingerprint_paths(root)
missing <- paths[!file.exists(paths)]
if (length(missing)) {
  stop("golden file(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
}

new <- golden_fingerprint_of(paths)
old <- leaf_behaviour_fingerprint()

target <- file.path(root, "R/fingerprint.R")
src <- readLines(target, warn = FALSE)
# The one string literal in the function body: a line that is nothing but an
# indented quoted token. Anchored tightly on purpose -- a looser pattern would also
# match a quoted string in the roxygen or in the helpers below.
hit <- grep('^  "[A-Za-z0-9]+"$', src)
if (length(hit) != 1L) {
  stop("expected exactly one fingerprint literal in R/fingerprint.R, found ",
       length(hit), call. = FALSE)
}
src[hit] <- sprintf('  "%s"', new)
writeLines(src, target)

if (identical(old, new)) {
  message("fingerprint unchanged: ", new)
} else {
  message("fingerprint ", old, " -> ", new)
  message("  from: ", paste(basename(paths), collapse = ", "))
  message("⚠️  the numbers moved. Say what moved, and by how much, in the commit.")
}
