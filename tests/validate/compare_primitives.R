# Localise the 1-2 ULP disagreement with plant to a single function (issue #13).
#
#   Rscript tests/validate/compare_primitives.R
#
# WHY THIS EXISTS, when compare_with_plant.R already exists. That script compares
# the SOLVE: 288 operating points, end to end, and reports 585 of 2592 values
# differing at 1-2 ULP. It cannot say where the difference comes from, and the
# issue explains why -- the nested solvers amplify perturbations up to about
# GSS_tol_abs (1e-3), so a rounding difference in one primitive arrives at the
# output having been through a golden-section search and a ci iteration, by which
# point it points nowhere.
#
# So this compares the PRIMITIVES, called directly at inputs fixed here rather
# than found by a solver. The functions are emitted in call-tree order:
#
#   tier 1  arrh_curve, peak_arrh_curve          pure arithmetic + exp/pow
#   tier 2  proportion_of_conductivity           the Weibull vulnerability curve
#   tier 3  assim_{rubisco,electron}_limited,    rational arithmetic, then a
#           assim_colimited                        quadratic root
#   tier 4  transpiration, stom_cond_CO2         first use of an odelia spline
#   tier 5  psi_stem_to_ci,                      contain their own iteration
#           transpiration_to_psi_stem
#
# The FIRST tier that disagrees is the answer. Clean tiers 1-4 with a dirty tier
# 5 would mean the arithmetic agrees and an iteration terminates differently,
# which is a different bug with a different fix.
#
# REFERENCE BUILD. This must run against plant built with plant's OWN leaf, from
# the commit this package was extracted from -- not the installed plant, which on
# this machine is the `feature/consume-leaf-package` build and therefore contains
# *this package*, making the comparison vacuous. The script refuses to run if it
# detects that. Build the reference with:
#
#   git -C ../plant worktree add /tmp/plant-ref 76df7169
#   R CMD INSTALL --no-docs -l /tmp/rlib-ref /tmp/plant-ref
#   R_LIBS=/tmp/rlib-ref Rscript tests/validate/compare_primitives.R

suppressMessages(library(plant))

# --- Refuse to compare the package against itself ----------------------------

leaf_model_h <- system.file("include", "plant", "leaf_model.h", package = "plant")
if (!nzchar(leaf_model_h)) {
  stop("plant is installed without its headers; cannot verify which leaf it uses")
}
if (any(grepl("leaf.hpp", readLines(leaf_model_h, warn = FALSE), fixed = TRUE))) {
  stop("The plant on this library path is the consume-leaf-package build: its\n",
       "  leaf_model.h includes <leaf.hpp>, so it IS this package and comparing\n",
       "  against it would be vacuous. Build the reference described above and\n",
       "  re-run with R_LIBS pointed at it.")
}

# --- Build the C++ side ------------------------------------------------------

inc <- function(p) system.file("include", package = p)
bin <- file.path(tempdir(), "leaf_primitives")
src <- "tests/validate/primitives.cpp"
if (!file.exists(src)) stop("run from the repository root: ", src, " not found")

message("Compiling ", src, " ...")
cxx <- paste(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX20"),
                     stdout = TRUE), collapse = " ")
std <- paste(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX20STD"),
                     stdout = TRUE), collapse = " ")
ok <- system2(cxx, c(std, "-O2", "-I", "inst/include",
                     "-isystem", shQuote(inc("odelia")),
                     "-isystem", shQuote(inc("BH")),
                     "-o", shQuote(bin), shQuote(src)))
if (ok != 0) stop("the primitives harness did not compile")

cpp <- read.delim(text = paste(system2(bin, stdout = TRUE), collapse = "
"),
                  stringsAsFactors = FALSE, colClasses = "character")

# Hex in, hex out. See the comment above row() in primitives.cpp: R reads its own
# decimal output one ULP off for about 18% of values, which is enough on its own
# to manufacture the disagreement this script exists to explain.
for (col in c("arg1", "arg2", "value")) cpp[[col]] <- as.numeric(cpp[[col]])

# --- The same calls, through plant's bindings --------------------------------

theta <- 0.000157; K_s <- 1.0; h <- 5.0
area_leaf <- 0.05; rho <- 608.0; a_bio <- 0.0245
ca <- 40.0; o2 <- 21.0; tleaf <- 25.0; patm <- 101.3

l <- plant::Leaf(
  vcmax_25 = 96, c = 2.680147, b = 3.898245, psi_crit = 5.870283,
  root_c = 2.680147, root_b = 3.898245, root_psi_crit = 5.870283,
  beta2 = 1.5, jmax_25 = 157.44, a = 0.30,
  curv_fact_elec_trans = 0.7, curv_fact_colim = 0.99,
  GSS_tol_abs = 1e-3, vulnerability_curve_ncontrol = 100,
  ci_abs_tol = 1e-3, ci_niter = 1000,
  g1_TF24 = 7.5, beta_R_H = 3.4e2, beta_R_V = 9.4e3
)
l$set_physiology(
  area_leaf = area_leaf, mass_root_prop = 1.0, rho = rho, a_bio = a_bio,
  PPFD = 900, psi_soil = 2.0, soil_depth = 1.0,
  leaf_specific_conductance_max = K_s * theta / h, atm_vpd = 2.0, ca = ca,
  sapwood_volume_per_leaf_area = theta * h, leaf_temp = tleaf,
  atm_o2_kpa = o2, atm_kpa = patm
)

# One evaluator per function, keyed the same way the C++ side prints them.
eval_row <- function(fn, a1, a2) {
  switch(fn,
    arrh_curve                = l$arrh_curve(a1, 96, a2),
    peak_arrh_curve           = l$peak_arrh_curve(43900, a1, a2, 200000, 640),
    proportion_of_conductivity = l$proportion_of_conductivity(a1),
    assim_rubisco_limited     = l$assim_rubisco_limited(a1),
    assim_electron_limited    = l$assim_electron_limited(a1),
    assim_colimited           = l$assim_colimited(a1),
    transpiration             = l$transpiration(a1, a2),
    stom_cond_CO2             = l$stom_cond_CO2(a1, a2),
    psi_stem_to_ci            = l$psi_stem_to_ci(a1, a2),
    transpiration_to_psi_stem = l$transpiration_to_psi_stem(a1, a2),
    stop("no evaluator for ", fn)
  )
}

cpp$plant <- vapply(seq_len(nrow(cpp)),
                    function(i) eval_row(cpp$function.[i], cpp$arg1[i], cpp$arg2[i]),
                    numeric(1))

# --- Compare -----------------------------------------------------------------

# ULP distance, which is the unit the disagreement is actually in. Relative
# difference would compress everything interesting into 1e-16 and hide whether a
# value is off by one representable step or by ten.
ulps <- function(a, b) {
  d <- abs(a - b)
  ifelse(d == 0, 0, d / pmax(.Machine$double.xmin, abs(a) * .Machine$double.eps))
}
cpp$ulp <- ulps(cpp$value, cpp$plant)

tiers <- c(arrh_curve = 1, peak_arrh_curve = 1,
           proportion_of_conductivity = 2,
           assim_rubisco_limited = 3, assim_electron_limited = 3,
           assim_colimited = 3,
           transpiration = 4, stom_cond_CO2 = 4,
           psi_stem_to_ci = 5, transpiration_to_psi_stem = 5)
cpp$tier <- tiers[cpp$function.]

by_fn <- do.call(rbind, lapply(split(cpp, cpp$function.), function(d) {
  data.frame(tier = d$tier[1], fn = d$function.[1], n = nrow(d),
             differing = sum(d$ulp > 0), worst_ulp = max(d$ulp),
             stringsAsFactors = FALSE)
}))
by_fn <- by_fn[order(by_fn$tier, by_fn$fn), ]

cat("\nPrimitive-by-primitive, in call-tree order\n")
cat("(plant's own leaf at 76df7169 vs this package, identical inputs)\n\n")
print(by_fn, row.names = FALSE)

dirty <- by_fn[by_fn$differing > 0, ]
cat("\n")
if (nrow(dirty) == 0) {
  cat("EVERY primitive is bit-identical, which is the expected result.\n\n",
      "The 1-2 ULP that compare_with_plant.R used to report was never in any of\n",
      "these functions, nor anywhere else in the model: it was R's decimal parser\n",
      "reading the golden file, which returns a double one ULP off the correctly\n",
      "rounded value for about 18% of inputs. Both scripts now read through\n",
      "tests/validate/tsv_to_hex.c and both report zero. See issue #13.\n\n",
      "So this is a regression check, not an investigation. If it ever reports a\n",
      "difference, something real has changed -- start at the lowest dirty tier.\n",
      sep = "")
} else {
  cat("First disagreement is in tier ", min(dirty$tier), ": ",
      paste(dirty$fn[dirty$tier == min(dirty$tier)], collapse = ", "), "\n", sep = "")
  worst <- cpp[which.max(cpp$ulp), ]
  cat(sprintf("Worst single value: %s(%.17g, %.17g)\n  cpp   %.17g\n  plant %.17g\n  %.1f ULP\n",
              worst$function., worst$arg1, worst$arg2, worst$value, worst$plant,
              worst$ulp))
}
