# The fingerprint against the files it claims to be derived from (#58).
#
# This is the enforcement, not a formality. `leaf_behaviour_fingerprint()` returns a
# committed CONSTANT -- it has to, because `R CMD INSTALL` does not put `tests/` in
# the library, so the installed package cannot see the golden files at all. A
# constant that nothing recomputes is a constant that goes stale, and a stale
# behaviour fingerprint is worse than none: a consumer's cache would look valid
# while the numbers underneath it had moved.
#
# So: recompute from the files on disk, and fail if they disagree. A PR that
# regenerates a golden file and forgets `Rscript tools/fingerprint.R` stops here.

# Both testthat entry points put the working directory at tests/testthat/, and the
# tarball ships tests/cpp/golden/ (it is not in .Rbuildignore -- tests/cpp.R needs
# it). Resolved relative to the test file rather than to getwd() so it does not
# depend on which of the two ran.
package_root <- function() normalizePath(file.path(testthat::test_path(), "..", ".."))

test_that("the fingerprint matches the golden files it is derived from", {
  paths <- phylloptim:::golden_fingerprint_paths(package_root())
  # If this skips, the test is measuring nothing -- so say which file is missing
  # rather than skipping silently.
  expect_true(all(file.exists(paths)),
              label = paste("golden files present:",
                            paste(basename(paths[!file.exists(paths)]),
                                  collapse = ", ")))

  expect_identical(
    leaf_behaviour_fingerprint(),
    phylloptim:::golden_fingerprint_of(paths),
    label = paste0("leaf_behaviour_fingerprint() is stale -- run ",
                   "`Rscript tools/fingerprint.R` and commit the result")
  )
})

test_that("the fingerprint is a stable, usable cache key", {
  fp <- leaf_behaviour_fingerprint()
  expect_type(fp, "character")
  expect_length(fp, 1L)
  # A cache key goes in file names and data frames, so it must be plain and short.
  expect_match(fp, "^[0-9a-f]{12}$")
  # Deterministic within a session, and it is a constant so there is nothing to
  # make it otherwise -- asserted because the whole value of the thing is that two
  # calls agree.
  expect_identical(fp, leaf_behaviour_fingerprint())
})

test_that("the digest depends on content, not on line endings or argument order", {
  # The CRLF case, which is why the digest hashes readLines() output rather than
  # bytes: a Windows checkout must not report different behaviour from a macOS one.
  lf <- tempfile(); crlf <- tempfile()
  on.exit(unlink(c(lf, crlf)), add = TRUE)
  writeLines(c("a\t1", "b\t2"), lf)
  con <- file(crlf, open = "wb"); writeChar("a\t1\r\nb\t2\r\n", con, eos = NULL)
  close(con)
  expect_identical(phylloptim:::golden_fingerprint_of(lf),
                   phylloptim:::golden_fingerprint_of(crlf))

  # Order-independent, since the caller should not have to remember which file
  # comes first.
  paths <- phylloptim:::golden_fingerprint_paths(package_root())
  skip_if_not(all(file.exists(paths)))
  expect_identical(phylloptim:::golden_fingerprint_of(paths),
                   phylloptim:::golden_fingerprint_of(rev(paths)))

  # ⚠️ A RENAME MUST NOT MOVE IT. #58's complaint about hashing the source tree is
  # precisely that "#47's rename would have silently changed the hash without
  # changing behaviour", so the file names are not in the digest. Same bytes under a
  # different name, same fingerprint.
  same_bytes <- tempfile(fileext = ".tsv")
  on.exit(unlink(same_bytes), add = TRUE)
  file.copy(paths[[1]], same_bytes, overwrite = TRUE)
  expect_identical(phylloptim:::golden_fingerprint_of(paths[[1]]),
                   phylloptim:::golden_fingerprint_of(same_bytes))

  # And it does move when a number does -- otherwise it is decoration. One digit of
  # one value in a copy of the real file.
  tweaked <- tempfile()
  on.exit(unlink(tweaked), add = TRUE)
  txt <- readLines(paths[[1]], warn = FALSE)
  txt[[2]] <- sub("1", "2", txt[[2]], fixed = TRUE)
  writeLines(txt, tweaked)
  expect_false(identical(phylloptim:::golden_fingerprint_of(paths[[1]]),
                         phylloptim:::golden_fingerprint_of(tweaked)))

  # Both files together must move if EITHER moves -- a fingerprint that only
  # tracked the solved outputs would miss a change confined to the gradient
  # composite, which is a route a calibration depends on.
  expect_false(identical(phylloptim:::golden_fingerprint_of(paths),
                         phylloptim:::golden_fingerprint_of(c(tweaked, paths[[2]]))))
})
