#!/usr/bin/env Rscript
# Resolve the model registry against Doxygen's XML, and fail if it cannot.
#
# Run by .github/workflows/docs.yml after the render, where the package is not
# installed: only R/registry-locations.R is sourced, so that file must stay
# free of anything else in the package (test-registry.R asserts it).
#
#   Rscript tools/registry_symbols.R docs/xml docs/html/model-registry.yml
#
# Three things are checked. Every `code` (and `also`) symbol the registry names
# resolves to a member in the XML, so a rename in the headers that the registry
# did not follow fails here. Every resolved line maps back to the source
# through the filter's line map, which refuses a header whose code lines the
# filter altered. And every parameter an equation `reads` is referenced by the
# equation's function (Doxygen's REFERENCES_RELATION), so an equation that
# quietly stopped reading a parameter, or started reading one not registered,
# fails here. The second argument, if given, receives the resolved registry:
# every parameter and equation with its file, line, body span, API anchor and
# initializer, for a document that wants to link without a Doxygen run of its
# own.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: registry_symbols.R <doxygen xml dir> [<output yml>]")
xml_dir <- args[1]
out_path <- if (length(args) >= 2L) args[2] else NA_character_

root <- normalizePath(file.path(dirname(xml_dir), ".."))
source(file.path(root, "R", "registry-locations.R"))

raw <- registry_read(file.path(root, "inst", "model-registry.yml"))
locations <- registry_locations(xml_dir,
                                include_dir = file.path(root, "inst", "include"),
                                filter = file.path(root, "tools", "doxygen_filter.awk"))
res <- registry_resolve(raw, locations = locations, unreferenced = "error")

cat(sprintf("model registry: %d parameters and %d equations resolved against %d symbols\n",
            nrow(res$parameters), nrow(res$equations), nrow(locations)))

if (!is.na(out_path)) {
  strip <- function(df) {
    df$replica <- NULL
    lapply(seq_len(nrow(df)), function(i) {
      row <- lapply(df[i, , drop = FALSE], function(col) {
        v <- if (is.list(col)) col[[1]] else col
        if (is.factor(v)) as.character(v) else v
      })
      row <- row[!vapply(row, function(v) length(v) == 1L && is.na(v), TRUE)]
      names(row) <- names(df)[names(df) %in% names(row)]
      row
    })
  }
  out <- list(package = raw$package,
              generated_by = "tools/registry_symbols.R",
              parameters = stats::setNames(strip(res$parameters), res$parameters$id),
              equations = stats::setNames(strip(res$equations), res$equations$id))
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  yaml::write_yaml(out, out_path)
  cat("wrote", out_path, "\n")
}
