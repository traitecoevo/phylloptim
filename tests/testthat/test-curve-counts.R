# The documents' count of the cost curves, held to the enumeration. Prose does
# not fail when a curve is added, so "eight" stood in four documents for a
# release after the ninth arrived.
#
# Matched on phrases about curves or models, since "eight" also counts the
# operating-point kinds. The documents count a subset of the curves in small
# numbers ("the four curves bounded there", "two curves that take one"), so a
# phrase with seven or more counts all of them and must be the enumeration's.
test_that("every document's count of the cost curves is the enumeration's", {
  root <- testthat::test_path("..", "..")
  skip_if_not(dir.exists(file.path(root, "vignettes")),
              "the source tree is not reachable from here")

  files <- c(file.path(root, c("README.md", "COMPARISON.md",
                               "inst/RcppR6_classes.yml",
                               ".claude/CLAUDE.md", ".claude/hazards.md")),
             list.files(file.path(root, "vignettes"), "\\.Rmd$", full.names = TRUE),
             list.files(file.path(root, "R"), "\\.R$", full.names = TRUE))
  files <- files[file.exists(files)]
  words <- c("seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen",
             "fourteen", "fifteen")
  pattern <- paste0("\\b(", paste(words, collapse = "|"),
                    ")( cost| optimality| stomatal)? (curves|models)\\b")
  found <- do.call(rbind, lapply(files, function(f) {
    text <- tolower(readLines(f, warn = FALSE))
    hit <- unlist(regmatches(text, gregexpr(pattern, text, perl = TRUE)))
    if (!length(hit)) return(NULL)
    data.frame(file = basename(f), phrase = hit,
               word = sub(" .*", "", hit), stringsAsFactors = FALSE)
  }))
  # The check is worthless if it matched nothing.
  expect_gt(nrow(found), 5L)
  n <- length(cost_curve_names())
  numbers <- c(seven = 7, eight = 8, nine = 9, ten = 10, eleven = 11,
               twelve = 12, thirteen = 13, fourteen = 14, fifteen = 15)
  wrong <- found[numbers[found$word] != n, ]
  expect_identical(nrow(wrong), 0L,
                   info = paste0(wrong$file, ": \"", wrong$phrase, "\"",
                                 collapse = "; "))
})
