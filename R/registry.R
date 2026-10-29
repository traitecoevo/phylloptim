# The model registry's live half: the registry as data frames, its replicas,
# the parameters' values read off the package, and the check that every
# replica agrees with the package. The pure half (reading the YAML, Doxygen
# locations) is R/registry-locations.R; the probes that drive the package for
# each equation are R/registry-probes.R.

#' The model registry
#'
#' phylloptim's temperature paths as one record: every equation (its LaTeX,
#' the C++ symbol that evaluates it, the parameters it reads, the bound fields
#' it writes, a base-R replica) and every parameter (symbol, units, class,
#' meaning, source, placeholder flag, the symbol that holds it). The file is
#' `inst/model-registry.yml`; documents are generated from it rather than
#' written beside it, and CI checks it (`tests/testthat/test-registry.R`,
#' `docs.yml`). Its header states the three rules it is kept by: no stored line
#' numbers or defaults, validated in CI, one registry per package.
#'
#' @param path the YAML file; the package's own by default
#' @return an object of class `model_registry`: `package`, `sections` (a data
#'   frame of `id` and `title`), `parameters` and `equations` (data frames, one
#'   row each, with list-columns for `reads`, `writes`, `also`, `uses` and
#'   `aliases`), `raw` (the validated list) and `path`
#' @seealso [registry_values()] for the values at render time,
#'   [registry_check()] for the replicas against the package,
#'   [registry_resolve()] for file, line and API link from Doxygen's XML
#' @examples
#' reg <- model_registry()
#' reg
#' subset(reg$equations, section == "respiration", c(id, title, code))
#' @export
model_registry <- function(path = system.file("model-registry.yml", package = "phylloptim")) {
  raw <- registry_read(path)
  tabs <- registry_tables(raw)
  structure(list(package = raw$package,
                 sections = data.frame(id = names(raw$sections),
                                       title = unname(unlist(raw$sections)),
                                       stringsAsFactors = FALSE),
                 parameters = tabs$parameters,
                 equations = tabs$equations,
                 raw = raw, path = path),
            class = "model_registry")
}

#' @export
print.model_registry <- function(x, ...) {
  cat("<model_registry> ", x$package, ": ", nrow(x$parameters), " parameters, ",
      nrow(x$equations), " equations\n", sep = "")
  for (i in seq_len(nrow(x$sections))) {
    s <- x$sections$id[i]
    cat(sprintf("  %-16s %3d parameters  %3d equations\n", s,
                sum(x$parameters$section == s), sum(x$equations$section == s)))
  }
  invisible(x)
}

#' The replicas, as functions
#'
#' Every equation's `replica` evaluated in one environment whose parent is the
#' base environment, so a replica can call another by its id
#' (`damage_rates` calls `impairment_rate`) and nothing else: no package
#' internals, which is the point of a replica.
#'
#' @param registry a [model_registry()]
#' @param id an equation id
#' @return `registry_replicas()`: an environment holding one function per
#'   equation id. `registry_replica()`: that one function.
#' @export
registry_replicas <- function(registry) {
  env <- new.env(parent = baseenv())
  for (i in seq_len(nrow(registry$equations))) {
    assign(registry$equations$id[i],
           eval(parse(text = registry$equations$replica[i]), env), envir = env)
  }
  env
}

#' @rdname registry_replicas
#' @export
registry_replica <- function(registry, id) {
  env <- registry_replicas(registry)
  if (!exists(id, envir = env, inherits = FALSE)) {
    stop("no equation `", id, "` in the registry", call. = FALSE)
  }
  get(id, envir = env, inherits = FALSE)
}

# The bound field a `phylloptim::Leaf::` symbol names, or NA.
.registry_field <- function(code) {
  ifelse(grepl("^phylloptim::Leaf::", code), sub("^phylloptim::Leaf::", "", code), NA_character_)
}

# The constexprs' values from Doxygen's initializers, which may name one another
# (`kg_to_mol_h2o = 1.0 / molar_mass_h2o`): evaluate in passes until nothing new
# resolves. Returns a named numeric vector keyed by qualified symbol.
.registry_initializer_values <- function(locations) {
  if (is.null(locations)) return(numeric())
  vars <- locations[locations$kind == "variable" & !is.na(locations$initializer), , drop = FALSE]
  env <- new.env(parent = baseenv())
  out <- numeric()
  pending <- seq_len(nrow(vars))
  for (pass in 1:4) {
    if (!length(pending)) break
    still <- integer()
    for (i in pending) {
      expr <- gsub("phylloptim::", "", vars$initializer[i], fixed = TRUE)
      v <- tryCatch(eval(parse(text = expr), env), error = function(e) NULL)
      if (is.numeric(v) && length(v) == 1L && is.finite(v)) {
        out[vars$symbol[i]] <- v
        assign(sub("^.*::", "", vars$symbol[i]), v, envir = env)
      } else {
        still <- c(still, i)
      }
    }
    pending <- still
  }
  out
}

#' The parameters' values, read at render time
#'
#' Adds `value` and `value_source` to the registry's parameters: a bound field
#' is read off `leaf`, a driver off [set_drivers()]'s default, a constexpr off
#' its initializer in Doxygen's XML when `locations` is given, and a literal
#' off the registry (the one kind of number it stores). Anything else is `NA`.
#' Nothing is read from a comment or a document, which is rule 1 of the
#' registry's header.
#'
#' @param registry a [model_registry()]
#' @param leaf the `Leaf` to read bound fields off; `NULL` reads none
#' @param locations the output of [registry_locations()], for the constexprs
#' @return the registry's `parameters` data frame with `value` and
#'   `value_source` (`"binding"`, `"set_drivers()"`, `"initializer"`,
#'   `"literal"` or `NA`)
#' @examples
#' v <- registry_values(model_registry())
#' v[v$section == "respiration", c("id", "symbol", "value", "value_source")]
#' @export
registry_values <- function(registry, leaf = leaf_model(), locations = NULL) {
  p <- registry$parameters
  consts <- .registry_initializer_values(locations)
  drivers <- formals(set_drivers)
  bound <- if (is.null(leaf)) character() else names(leaf)
  value <- rep(NA_real_, nrow(p))
  src <- rep(NA_character_, nrow(p))
  for (i in seq_len(nrow(p))) {
    if (!is.na(p$literal[i])) {
      value[i] <- p$literal[i]; src[i] <- "literal"
    } else if (!is.na(p$driver[i]) && p$driver[i] %in% names(drivers)) {
      value[i] <- eval(drivers[[p$driver[i]]]); src[i] <- "set_drivers()"
    } else if (!is.na(f <- .registry_field(p$code[i])) && f %in% bound) {
      v <- leaf[[f]]
      value[i] <- if (is.logical(v)) as.numeric(v) else v
      src[i] <- "binding"
    } else if (p$code[i] %in% names(consts)) {
      value[i] <- consts[[p$code[i]]]; src[i] <- "initializer"
    }
  }
  p$value <- value
  p$value_source <- src
  p
}

#' Check every replica against the package
#'
#' Runs each equation's probe (`R/registry-probes.R`), which drives a `Leaf`
#' and reads the quantity back, and compares it with the registry's replica of
#' the printed equation on a grid. One row per equation: the worst absolute and
#' relative disagreement, where relative is to the package value floored at 1 so
#' a quantity passing through zero is judged absolutely. An equation passes if
#' either is within its tolerance: `tol` unless the entry sets its own (an
#' identity or a derivative, which a finite difference cannot pin to 1e-10).
#'
#' @param registry a [model_registry()]
#' @param temps the leaf temperatures of the grid, deg C
#' @param tol the tolerance for an entry that sets none
#' @param assert stop, naming the failures, if any entry fails
#' @return a data frame: `id`, `kind`, `n` (comparisons made), `max_abs`,
#'   `max_rel`, `tol`, `pass`
#' @examples
#' \donttest{
#' registry_check()
#' }
#' @export
registry_check <- function(registry = model_registry(), temps = seq(0, 50, by = 2.5),
                           tol = 1e-10, assert = FALSE) {
  rep <- registry_replicas(registry)
  eq <- registry$equations
  rows <- lapply(seq_len(nrow(eq)), function(i) {
    id <- eq$id[i]
    probe <- .registry_probes[[id]]
    if (is.null(probe)) {
      stop("registry_check(): equation `", id, "` has no probe in R/registry-probes.R",
           call. = FALSE)
    }
    cmp <- probe(rep, temps)
    pkg <- as.numeric(cmp$pkg)
    mine <- as.numeric(cmp$mine)
    if (!length(pkg) || length(pkg) != length(mine)) {
      stop("registry_check(): the probe for `", id, "` returned ", length(pkg),
           " package values against ", length(mine), " replica values", call. = FALSE)
    }
    if (!all(is.finite(pkg)) || !all(is.finite(mine))) {
      stop("registry_check(): a non-finite value in the probe for `", id, "`", call. = FALSE)
    }
    a <- abs(pkg - mine)
    r <- a / pmax(abs(pkg), 1)
    t <- if (is.na(eq$tol[i])) tol else eq$tol[i]
    data.frame(id = id, kind = eq$kind[i], n = length(pkg), max_abs = max(a),
               max_rel = max(r), tol = t, pass = max(r) <= t || max(a) <= t,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (assert && !all(out$pass)) {
    stop("registry_check(): replicas disagree with the package beyond their tolerance:\n",
         paste(utils::capture.output(print(out[!out$pass, ])), collapse = "\n"),
         call. = FALSE)
  }
  out
}
