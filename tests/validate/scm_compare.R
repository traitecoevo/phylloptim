# Compare two scm_regression.R outputs, exactly.
#
#   Rscript tests/validate/scm_compare.R /tmp/scm_ref.rds /tmp/scm_consume.rds
#
# Walks both structures in parallel and reports, per leaf node, the number of
# differing values and the worst relative difference. Exact comparison: the
# interesting question is whether the 1-ULP leaf difference survives the SCM's
# adaptive stepper and discrete node schedule, so rounding to a tolerance up front
# would discard the answer.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
a <- readRDS(args[[1]])
b <- readRDS(args[[2]])

cat("A:", a$plant_lib, "\nB:", b$plant_lib, "\n\n")

rows <- list()

walk <- function(x, y, path) {
  if (is.list(x) && is.list(y)) {
    kx <- names(x); ky <- names(y)
    if (!identical(kx, ky)) {
      rows[[length(rows) + 1L]] <<- data.frame(
        node = path, status = "STRUCTURE DIFFERS", n_diff = NA_integer_,
        max_rel = NA_real_)
      return(invisible())
    }
    if (is.null(kx)) kx <- seq_along(x)
    for (k in kx) walk(x[[k]], y[[k]], paste0(path, "$", k))
    return(invisible())
  }
  if (is.numeric(x) && is.numeric(y)) {
    if (length(x) != length(y)) {
      rows[[length(rows) + 1L]] <<- data.frame(
        node = path, status = sprintf("LENGTH %d vs %d", length(x), length(y)),
        n_diff = NA_integer_, max_rel = NA_real_)
      return(invisible())
    }
    if (!length(x)) return(invisible())
    xv <- as.numeric(x); yv <- as.numeric(y)
    both_nan <- is.na(xv) & is.na(yv)
    d <- which(!both_nan & !(xv == yv))
    rel <- if (length(d)) {
      max(abs(xv[d] - yv[d]) / pmax(abs(yv[d]), .Machine$double.xmin))
    } else 0
    rows[[length(rows) + 1L]] <<- data.frame(
      node = path, status = if (length(d)) "differs" else "identical",
      n_diff = length(d), max_rel = rel)
    return(invisible())
  }
  rows[[length(rows) + 1L]] <<- data.frame(
    node = path, status = if (identical(x, y)) "identical" else "DIFFERS",
    n_diff = NA_integer_, max_rel = NA_real_)
}

# data.frames are lists, and species output is one, so walking them column-wise is
# what we want -- it localises any difference to a named variable.
walk(a[names(a) != "plant_lib"], b[names(b) != "plant_lib"], "")

tab <- do.call(rbind, rows)
tab <- tab[order(-ifelse(is.na(tab$max_rel), Inf, tab$max_rel)), ]

n_bad <- sum(tab$status != "identical", na.rm = TRUE)
cat(sprintf("%d of %d numeric nodes differ\n\n", n_bad, nrow(tab)))
print(head(tab, 25), row.names = FALSE)

worst <- suppressWarnings(max(tab$max_rel, na.rm = TRUE))
if (!is.finite(worst)) worst <- 0
struct <- any(grepl("STRUCTURE|LENGTH", tab$status))

cat("\n")
if (struct) {
  cat("FAIL: the two runs differ structurally -- a different number of ODE steps or\n")
  cat("nodes. That would mean the leaf difference changed a discrete decision in the\n")
  cat("solver, which is the outcome this check exists to detect.\n")
  quit(status = 1)
}
if (n_bad == 0L) {
  cat("PASS: the two SCM runs are bit-identical throughout.\n")
  quit(status = 0)
}
cat(sprintf("worst relative difference: %.3e (%.0f ULP)\n",
            worst, worst / .Machine$double.eps))
if (worst <= 1e-10) {
  cat("PASS: differences are at rounding level and no discrete decision changed.\n")
  quit(status = 0)
}
cat("FAIL: too large to be rounding. Investigate before trusting the swap.\n")
quit(status = 1)
