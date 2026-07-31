# PLAN.md item 1: cross-check this package against plant's compiled build.
#
#   Rscript tests/validate/compare_with_plant.R
#
# Runs plant's own `Leaf` -- the compiled class in the installed plant package,
# not a reimplementation -- over exactly the grid in tests/cpp/test_golden.cpp,
# and compares against the golden file this package generated.
#
# RESULT AS OF 2026-07-31: 585 of 2592 value comparisons differ, every one of them
# at 1-2 ULP (worst relative difference 2.2e-16; double eps is 2.2e-16). The
# extraction is faithful; see PLAN.md item 1 for the full argument, in short:
#
#   * The SOURCE is arithmetically identical. A normalised function-body diff over
#     all 44 shared functions found only three `size_t` -> `int` loop counters,
#     which cannot change floating-point arithmetic. (Plus
#     transpiration_full_integration, adaptive Simpson by design, not on this path.)
#   * The CAUSE of the residual is NOT established. Ruled out by experiment: the
#     plant version/branch, compiler flags (-std=c++20 vs gnu++20, -g, -O0..-O3,
#     -DNDEBUG, -fPIC, -ffp-contract on/off), inlining (-fno-inline), and the odelia
#     header version (local checkout vs installed). An earlier version of this
#     comment blamed translation-unit structure; that is CONTRADICTED by
#     scm_regression.R, which found plant-with-its-own-leaf and
#     plant-with-this-package bit-identical -- and that comparison DOES change TU
#     structure. So the honest statement is: 1 ULP, cause unknown, and it is
#     confined to this standalone harness rather than to the package.
#
# WHAT IS ESTABLISHED, and it is the part that matters: swapping plant's own leaf
# for this package changes NOTHING. plant's test-leaf.r passes unchanged (218
# expectations) and a full SCM regression is bit-identical across 78 of 78 recorded
# nodes, including the whole collected trajectory. See scm_regression.R /
# scm_compare.R. The 1 ULP here is a property of comparing a standalone C++ binary
# against plant's R-bound build, not of the extraction.
#
# So do NOT expect zero. Expect <= a few ULP, and investigate anything larger:
# these nested solvers amplify perturbations up to about GSS_tol_abs (1e-3), so a
# genuine arithmetic difference shows up at 1e-4, not 1e-16. Demonstrated -- forcing
# -ffp-contract=off on one side only moves the disagreement to 3e-4.
#
# Run this against `main`, whose set_physiology signature still matches plant's.
# The `feature/api-cleanup` branch deliberately diverges (different signature), so
# this script does not apply there; its deviations are measured in its own commits.
#
# REFERENCE BUILD. Compare against plant built from the commit this package was
# extracted from, NOT whatever happens to be installed -- the installed plant may
# be a different branch. Recreate it with:
#
#   git -C ../plant worktree add /tmp/plant-ref 76df7169
#   R CMD INSTALL --no-docs -l /tmp/rlib-ref /tmp/plant-ref
#   R_LIBS=/tmp/rlib-ref Rscript tests/validate/compare_with_plant.R
#
# (The first attempt at this compared against an installed plant carrying the ATLS
# thermal-damage layer, which is a different branch. It gave identical numbers,
# because ATLS is default-off and genuinely bit-identical when off -- but that was
# luck, not method.)

suppressMessages(library(plant))

golden_path <- "tests/cpp/golden/operating_points.tsv"
if (!file.exists(golden_path)) {
  stop("Run from the repository root: ", golden_path, " not found")
}

# --- the grid, mirroring tests/cpp/test_golden.cpp exactly --------------------

psi_soils <- c(0.5, 1.0, 2.0, 3.0, 4.0, 6.0)
ppfds <- c(100, 500, 900, 1500)
vpds <- c(0.5, 1.0, 2.0, 4.0)
layer_counts <- c(1L, 3L, 5L)

# Fixed drivers and traits, from plant's tests/testthat/test-leaf.r.
theta <- 0.000157
K_s <- 1.0
h <- 5.0
area_leaf <- 0.05
rho <- 608.0
a_bio <- 0.0245
ca <- 40.0
o2 <- 21.0
tleaf <- 25.0
patm <- 101.3

# Every constructor argument passed explicitly, matching the C++ default
# constructor value for value. Relying on R-side defaults would leave the
# comparison hostage to a defaults mismatch, which is exactly the sort of thing
# this script exists to detect.
new_leaf <- function() {
  plant::Leaf(
    vcmax_25 = 96, c = 2.680147, b = 3.898245, psi_crit = 5.870283,
    root_c = 2.680147, root_b = 3.898245, root_psi_crit = 5.870283,
    beta2 = 1.5, jmax_25 = 157.44, a = 0.30,
    curv_fact_elec_trans = 0.7, curv_fact_colim = 0.99,
    GSS_tol_abs = 1e-3, vulnerability_curve_ncontrol = 100,
    ci_abs_tol = 1e-3, ci_niter = 1000,
    g1_TF24 = 7.5, beta_R_H = 3.4e2, beta_R_V = 9.4e3
  )
}

solve_one <- function(psi_soil, ppfd, vpd, layers) {
  l <- new_leaf()
  i <- seq_len(layers)
  ps <- psi_soil + 0.25 * (i - 1)
  depth <- 1.0 * i
  root <- rep(1.0 / layers, layers)

  l$set_physiology(
    area_leaf = area_leaf, mass_root_prop = root, rho = rho, a_bio = a_bio,
    PPFD = ppfd, psi_soil = ps, soil_depth = depth,
    leaf_specific_conductance_max = K_s * theta / h,
    atm_vpd = vpd, ca = ca, sapwood_volume_per_leaf_area = theta * h,
    leaf_temp = tleaf, atm_o2_kpa = o2, atm_kpa = patm
  )
  l$find_root_collar_psi()

  cons <- l$soil_consumption_
  data.frame(
    psi_soil = psi_soil, ppfd = ppfd, vpd = vpd, layers = layers,
    psi_stem = l$opt_psi_stem_, collar = l$root_collar_psi_,
    ci = l$ci_, assim = l$assim_colimited_,
    transpiration = l$transpiration_, gc = l$stom_cond_CO2_,
    profit = l$profit_, e_up = l$E_up_,
    uptake = sum(cons[is.finite(cons)])
  )
}

grid <- expand.grid(layers = layer_counts, vpd = vpds, ppfd = ppfds,
                    psi_soil = psi_soils)
# expand.grid varies the FIRST column fastest; the C++ loops nest
# psi_soil > ppfd > vpd > layers, with layers innermost. Listing them in
# reverse above therefore reproduces the C++ row order.
grid <- grid[, c("psi_soil", "ppfd", "vpd", "layers")]

cat(sprintf("Running plant's Leaf over %d operating points...\n", nrow(grid)))
plant_rows <- do.call(rbind, Map(solve_one, grid$psi_soil, grid$ppfd,
                                 grid$vpd, grid$layers))

# --- compare ------------------------------------------------------------------

golden <- read.delim(golden_path)
stopifnot(nrow(golden) == nrow(plant_rows))

# The inputs must line up before comparing outputs, or we would be comparing
# different operating points and calling it agreement.
for (key in c("psi_soil", "ppfd", "vpd", "layers")) {
  if (!isTRUE(all.equal(golden[[key]], plant_rows[[key]], tolerance = 0))) {
    stop("grid mismatch in '", key, "': the R grid is not in the same order as ",
         "the C++ one, so the comparison would be meaningless")
  }
}

fields <- c("psi_stem", "collar", "ci", "assim", "transpiration", "gc",
            "profit", "e_up", "uptake")

# NaN == NaN for our purposes: shut-down points legitimately carry the NA
# sentinel on `main` (the shutdown-state leak, PLAN item 2).
differs <- function(a, b) !((is.nan(a) & is.nan(b)) | (!is.nan(a) & !is.nan(b) & a == b))

summary_rows <- list()
total_mismatch <- 0L
for (f in fields) {
  a <- plant_rows[[f]]
  b <- golden[[f]]
  bad <- which(differs(a, b))
  total_mismatch <- total_mismatch + length(bad)
  rel <- if (length(bad)) {
    max(abs(a[bad] - b[bad]) / pmax(abs(b[bad]), .Machine$double.xmin), na.rm = TRUE)
  } else 0
  summary_rows[[f]] <- data.frame(field = f, mismatches = length(bad),
                                  max_rel_diff = rel)
}
summ <- do.call(rbind, summary_rows)

cat("\n")
print(summ, row.names = FALSE)
cat(sprintf("\n%d of %d value comparisons differ.\n",
            total_mismatch, nrow(golden) * length(fields)))

# The bar is a few ULP, not zero -- see the header. Bit-identity is unachievable
# across a header-only single-TU build and plant's separate-TU build, and demanding
# it would just mean a permanently red check.
ulp_tol <- 1e-14   # ~45 ULP of headroom; a real arithmetic difference lands at 1e-4
worst_rel <- max(summ$max_rel_diff)

if (total_mismatch == 0L) {
  cat("\nPASS: bit-identical to plant across the whole grid.\n")
  quit(status = 0)
}
if (worst_rel <= ulp_tol) {
  cat(sprintf(paste0("\nPASS: agreement to %.1e relative (%.0f ULP), within the %.0e bar.\n",
                     "The source is arithmetically identical (see the header); the\n",
                     "cause of the 1-ULP residual is NOT established, but it is\n",
                     "confined to this standalone harness -- swapping the package into\n",
                     "plant is bit-identical (scm_compare.R). Anything above ~1e-12\n",
                     "would be a real difference and should be investigated.\n"),
              worst_rel, worst_rel / .Machine$double.eps, ulp_tol))
  quit(status = 0)
}

# Show the worst offenders so the failure is diagnosable, not just reported.
cat("\nWorst mismatches:\n")
worst <- NULL
for (f in fields) {
  a <- plant_rows[[f]]; b <- golden[[f]]
  bad <- which(differs(a, b))
  if (!length(bad)) next
  rel <- abs(a[bad] - b[bad]) / pmax(abs(b[bad]), .Machine$double.xmin)
  o <- bad[order(rel, decreasing = TRUE)][seq_len(min(3, length(bad)))]
  worst <- rbind(worst, data.frame(
    field = f, psi_soil = golden$psi_soil[o], ppfd = golden$ppfd[o],
    vpd = golden$vpd[o], layers = golden$layers[o],
    plant = a[o], leaf = b[o],
    rel_diff = abs(a[o] - b[o]) / pmax(abs(b[o]), .Machine$double.xmin)))
}
print(head(worst[order(worst$rel_diff, decreasing = TRUE), ], 20), row.names = FALSE)

cat(sprintf("\nFAIL: worst relative difference %.3e exceeds the %.0e bar.\n", worst_rel, ulp_tol))
cat("That is too large to be reassociation. Check, in order: that plant was built\n")
cat("from the commit this package was extracted from (see the header -- a different\n")
cat("branch is the likeliest cause); that the grid and constructor arguments still\n")
cat("match tests/cpp/test_golden.cpp; and only then suspect a conversion bug.\n")
quit(status = 1)
