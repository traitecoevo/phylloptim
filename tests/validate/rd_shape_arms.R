# The #41 respiration-shape revalidation, ON THE PLANT SIDE.
#
#   Rscript tests/validate/rd_shape_arms.R
#
# The companion to rd_shape_grid.cpp, which does the same A/B on this package's own
# grid. This one runs plant, because plant is the consumer and an integrated
# demographic output is not predictable from a per-leaf one.
#
# ⚠️ DO NOT RUN THIS PINNED AT 25 C AND CALL IT REVALIDATED, which was the obvious
# temptation and is worthless. #41 is inert at 25 C BY CONSTRUCTION -- the
# respiration reference value is defined there -- so a bit-identical result at 25 C
# proves only that the reference point survived, which is a control and not a
# validation. The information is at the temperatures plant actually reaches.
#
# ⚠️ AND PLANT ITSELF PINS 25 C AT ITS DEFAULTS, which is the fact that decides how
# to read all of this. `TF24_Environment` sets `leaf_temp` as a CONSTANT extrinsic
# driver at 25, and `use_energy_balance` defaults to 0 (Tleaf = Tair). So at plant's
# default configuration #41 changes nothing, and that is provable rather than
# observed. The exposure is in the two configurations that leave 25 C: a
# non-default `leaf_temp` driver, and the Penman-Monteith path (plant #523).
#
# --- what was measured, 2026-08-10 --------------------------------------------
#
# TF24 SCM offspring production, one species, max_patch_lifetime = 5:
#
#   plant configuration                    phylloptim master     this branch   change
#   default (leaf_temp = 25, PM off)     81.857087216483691  81.857087216483691  BIT-IDENTICAL
#   leaf_temp driver = 30 C              24.283073719078399  22.369992929645893  -7.9%
#   leaf_temp driver = 35 C              0.00021301630930330433  6.4471696967520308e-07  -99.7% (/331)
#   Penman-Monteith energy balance on     2.2034981948510728   1.1638497249192270  -47.2%
#
# The PM row is the one to carry forward: it is where plant is heading, and #41
# halves offspring production there. The 35 C row is a factor of 331, but on a
# population already four orders below replacement, so it is a large change to a
# number that was already saying "extinct".
#
# Leaf temperature on the PM path, over PAR {400,1000,2000} x Tair {20..40} x
# VPD {1,2,3}: Tleaf runs 19.1 to 62.0 C, i.e. up to +22 K ABOVE air temperature.
# That is why the PM arm moves so much more than its air temperature suggests --
# a PM leaf at Tair 25 sits at 32.9 C, and mean A there falls 23%.
#
# ⚠️ On the PM path `Tleaf` ITSELF moved, on 22 of 90 cells. The energy balance
# couples to transpiration, so a change in R_d feeds back into leaf temperature
# rather than only into carbon. Do not treat leaf temperature as an input on that
# path.
#
# --- how to build the two arms ------------------------------------------------
#
# Four installs: phylloptim twice and plant twice, because plant compiles
# phylloptim's headers in. ⚠️ `R_LIBS` FALLS BACK TO THE SITE LIBRARY SILENTLY if
# the intended one is not ready, which happened once while writing this and
# produced a plausible file describing the wrong package -- hence the hard check
# below, and print find.package() in anything you add.
#
#   SP=/tmp/rd41; W=$(git rev-parse --show-toplevel)
#   SITE=$(Rscript -e 'cat(.Library.site[[1]])')
#   mkdir -p $SP/{lib-pm,lib-branch,lib-plant-master,lib-plant-branch}
#   git -C $W archive origin/master | tar -x -C $SP/pm-src     # mkdir first
#   R_LIBS=$SITE R CMD INSTALL --no-docs -l $SP/lib-pm     $SP/pm-src
#   R_LIBS=$SITE R CMD INSTALL --no-docs -l $SP/lib-branch $W
#   # plant, twice, from a CLEAN source copy each time -- a stale src/*.o is a
#   # segfault that looks like a real bug (and would silently be the other arm)
#   git -C ../plant archive <plant-sha> --prefix=plant/ | tar -x -C $SP
#   R_LIBS="$SP/lib-pm:$SITE"     R CMD INSTALL --no-docs -l $SP/lib-plant-master $SP/plant
#   # ...fresh copy, then...
#   R_LIBS="$SP/lib-branch:$SITE" R CMD INSTALL --no-docs -l $SP/lib-plant-branch $SP/plant2
#
# then run this script once per arm with R_LIBS set to that arm and EXPECT_LIB=$SP.
#
# plant needs no source change for this: it never calls `set_traits` (so the
# fourteenth argument is invisible to it) and binds only `R_d_` of the respiration
# fields. Checked against plant's inst/RcppR6_classes.yml.
#
# --- what this does NOT cover -------------------------------------------------
#
# `compare_with_plant.R` compares against plant's OWN Leaf -- an independent
# implementation -- which is a stronger check than either arm here, and it still
# does not run. Its header blamed the constructor signature; that is no longer the
# obstacle. Against plant at 76df7169 the `Leaf()` and `set_physiology()` calls
# both work, and what actually stops it is that the reference has
# `root_collar_psi_` where this package has `opt_root_psi_` -- the #25 rename,
# which flipped the sign -- on top of results deliberately moved by #15, 11a, 11b,
# #25, #77 and #84. Reviving it is issue #64 and needs a decision about the
# reference, not a rename. Recorded here so the next attempt starts from the real
# obstacle.
suppressMessages(library(plant))

# ⚠️ HARD-FAIL ON THE WRONG LIBRARY, for the reason in the header. A warning would
# not do: every number below is plausible for either arm.
want <- Sys.getenv("EXPECT_LIB")
if (!nzchar(want)) {
  stop("set EXPECT_LIB to the directory both packages must resolve inside")
}
for (pkg in c("plant", "phylloptim")) {
  if (!startsWith(normalizePath(find.package(pkg)), normalizePath(want))) {
    stop(pkg, " resolved outside ", want, ": ", find.package(pkg))
  }
}
cat("plant:     ", find.package("plant"), "\n")
cat("phylloptim:", find.package("phylloptim"), "\n\n")

# One SCM run. The setup is test-strategy-tf24.R's own regression case, so the
# default arm is comparable with a number plant already asserts.
offspring <- function(leaf_temp = NULL, energy_balance = FALSE) {
  p0 <- scm_base_parameters("TF24")
  p0$max_patch_lifetime <- 5
  p1 <- add_strategies(p0, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                       hyperpar = TF24_hyperpar, birth_rate = list(20))
  if (energy_balance) {
    # ⚠️ THE FLAG LIVES IN `$pars`, and assigning it at the strategy's top level
    # is accepted and silently ignored: `p1$strategies[[1]]$use_energy_balance <-
    # 1` adds a list element nothing reads, and the run comes back EXACTLY equal
    # to the default arm -- which reads as "PM changes nothing" rather than as a
    # broken switch. Asserted rather than trusted.
    st <- p1$strategies[[1]]
    st$pars$use_energy_balance <- 1.0
    p1$strategies[[1]] <- st
    stopifnot(p1$strategies[[1]]$pars$use_energy_balance == 1.0)
  }
  env <- Environment("TF24")
  if (!is.null(leaf_temp)) {
    env$extrinsic_drivers_set_constant("leaf_temp", leaf_temp)
  }
  run_scm(p1, env, Control())$offspring_production
}

cases <- list(
  list(label = "default (leaf_temp driver = 25, PM off)", temp = NULL, pm = FALSE),
  list(label = "leaf_temp driver = 30",                   temp = 30,   pm = FALSE),
  list(label = "leaf_temp driver = 35",                   temp = 35,   pm = FALSE),
  list(label = "Penman-Monteith energy balance on",       temp = NULL, pm = TRUE)
)
for (cs in cases) {
  v <- offspring(cs$temp, cs$pm)
  # %.17g, so the two arms can be compared bit-for-bit rather than to the 2e-2
  # the plant test uses -- the point of the default arm is exact equality.
  cat(sprintf("%-40s %s\n", cs$label, paste(sprintf("%.17g", v), collapse = " ")))
}
