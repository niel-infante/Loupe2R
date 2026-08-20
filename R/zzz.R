.pkg_env <- new.env(parent = emptyenv())
.pkg_env$cloupe_pkg <- "~/Tools/cloupe"

#' Set the path to the local cellgeni/cloupe Python package
#'
#' The default is `~/Tools/cloupe`. Use this function to point Loupe2R at a
#' different installation. Changes persist for the R session only.
#'
#' To make the setting permanent, add to your `.Rprofile`:
#'
#' ```r
#' options(loupe2r.cloupe_pkg = "/path/to/cloupe")
#' ```
#'
#' @param path Character. Path to the root of the cloupe Python package.
#' @return The path, invisibly.
#' @export
set_cloupe_path <- function(path) {
  .pkg_env$cloupe_pkg <- path
  invisible(path)
}

#' Resolve and validate the configured cloupe package path
#'
#' Internal. Expands and checks `.pkg_env$cloupe_pkg`, failing fast with an
#' actionable message rather than letting reticulate's import fail deep
#' inside Python with a cryptic traceback.
#' @noRd
.resolve_cloupe_pkg <- function() {
  path   <- path.expand(.pkg_env$cloupe_pkg)
  marker <- file.path(path, "cloupe", "__init__.py")
  if (!dir.exists(path) || !file.exists(marker)) {
    stop(
      "Cannot find the cellgeni/cloupe Python package at:\n  ", path, "\n\n",
      "Loupe2R requires a local installation of cellgeni/cloupe (AGPL-3.0,\n",
      "https://github.com/cellgeni/cloupe), which is NOT bundled with this\n",
      "R package. Install it, e.g.:\n",
      "  pip install git+https://github.com/cellgeni/cloupe.git\n",
      "into the Python environment reticulate will use, then either:\n",
      "  - call set_cloupe_path(\"/path/to/cloupe\"), or\n",
      "  - set options(loupe2r.cloupe_pkg = \"/path/to/cloupe\") in .Rprofile.",
      call. = FALSE
    )
  }
  normalizePath(path, mustWork = TRUE)
}

.onLoad <- function(libname, pkgname) {
  opt <- getOption("loupe2r.cloupe_pkg")
  if (!is.null(opt)) .pkg_env$cloupe_pkg <- opt
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Loupe2R imports .cloupe files via an unofficial, reverse-engineered\n",
    "parser (cellgeni/cloupe, AGPL-3.0, described upstream as a 'PoC'). A\n",
    ".cloupe file captures only ONE bin resolution (not the multi-resolution\n",
    "structure Seurat::Load10X_Spatial(bin.size=c(8,16)) gives you), and\n",
    "array_row/array_col are a rounded approximation, not read from an\n",
    "authoritative source. See the package README/vignette before trusting\n",
    "results beyond exploration. (Silence with suppressPackageStartupMessages().)"
  )
}
