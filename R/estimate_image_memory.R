#' Estimate tissue-image memory usage across candidate resolutions
#'
#' Scopes the memory tissue-image reconstruction and reload will need at a
#' range of \code{px_per_bin} choices, without doing a real extraction.
#' Calls \code{cloupe_extract}'s own \code{estimate_image_memory()} (Python)
#' via reticulate, which opens the file with the count matrix itself
#' skipped entirely and reads only the tile-pyramid's index metadata --
#' never actual tile pixel bytes -- so this is cheap enough to run before
#' committing to a real \code{cloupe_to_seurat()} call.
#'
#' This is scoped to the tissue image specifically, not total memory use.
#' The non-image part (count matrix, barcodes, coordinates) doesn't vary
#' with resolution choice, and for real Visium HD datasets in the
#' few-hundred-thousand-barcode range has been measured at roughly 6-10 GB
#' peak RSS regardless (see \code{cloupe_to_seurat()}'s own documentation).
#' The image is the part a caller actually has a lever over via
#' \code{px_per_bin}, which is what this function reports on.
#'
#' Two memory numbers are reported per resolution, since they come from two
#' different, separately crash-prone steps: \code{python_peak_gb} is
#' \code{cloupe_extract}'s own stitching peak (Pillow), while
#' \code{r_reload_peak_gb} models \code{png::readPNG()}'s double-precision
#' expansion when the resulting image is loaded back into this package's
#' Seurat object -- the tighter constraint in practice, and the step that
#' has actually crashed on a real oversized image (see this package's own
#' \code{px_per_bin} documentation on \code{\link{cloupe_to_seurat}}).
#' Both are approximations from a formula, not profiled measurements.
#'
#' @param cloupe_path   Path to the .cloupe file.
#' @param px_per_bin    Candidate resolutions to report on, as a list (not a
#'                      plain numeric vector -- one candidate can be
#'                      \code{NULL} for native resolution, which a numeric
#'                      vector can't hold alongside numbers). Default
#'                      \code{list(NULL, 16, 8, 4, 2, 1)}, in the same units
#'                      as \code{cloupe_to_seurat()}'s \code{px_per_bin}.
#' @param version_check If TRUE (default), abort with an error when the
#'                      file's internal .cloupe format version(s) haven't
#'                      been validated against real paired SpaceRanger
#'                      output. Same meaning as \code{cloupe_to_seurat()}'s.
#' @param condaenv      Name of conda environment to activate before the
#'                      estimate (NULL = use current reticulate Python).
#'
#' @return A data.frame with one row per requested \code{px_per_bin}:
#'   \describe{
#'     \item{px_per_bin}{The requested value (\code{NA} for native).}
#'     \item{image_width, image_height}{The resulting image dimensions.}
#'     \item{clamped_to_native}{TRUE if this request was coarser than the
#'       file's native resolution already is, so no downsampling actually
#'       happens (see \code{target_image_size()} in \code{cloupe_extract}).}
#'     \item{python_peak_gb}{Approximate Pillow-side peak, in GB.}
#'     \item{r_reload_peak_gb}{Approximate \code{png::readPNG()} reload
#'       peak, in GB -- the tighter constraint in practice.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(reticulate)
#' use_condaenv("loupe2py", required = TRUE)
#'
#' estimate_image_memory("path/to/sample.cloupe")
#' }
#'
#' @export
estimate_image_memory <- function(
  cloupe_path,
  px_per_bin    = list(NULL, 16, 8, 4, 2, 1),
  version_check = TRUE,
  condaenv      = NULL
) {
  if (!is.null(condaenv))
    reticulate::use_condaenv(condaenv, required = TRUE)

  cloupe_path <- normalizePath(path.expand(cloupe_path), mustWork = TRUE)

  cloupe_extract <- tryCatch(
    reticulate::import("cloupe_extract"),
    error = function(e) stop(
      "Failed to import cloupe_extract. Install it with:\n",
      "  pip install \"cloupe_extract @ git+https://github.com/niel-infante/Loupe2Py.git#subdirectory=cloupe_extract\"\n",
      "into the Python environment reticulate will use.\n\n",
      "Original error: ", conditionMessage(e)
    )
  )

  # NULL elements of px_per_bin convert to Python None automatically via
  # reticulate's list conversion; version_check=TRUE raises a Python
  # UnvalidatedFormatVersionError here, which reticulate propagates as an
  # R error, same as cloupe_to_seurat()'s own version_check.
  results <- cloupe_extract$estimate_image_memory(
    cloupe_path, px_per_bin = px_per_bin, version_check = version_check
  )

  .image_memory_table(results)
}
