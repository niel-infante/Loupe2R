#' Detect the mitochondrial gene prefix for a set of gene names
#'
#' Internal. Returns "^MT-" if any human-style mitochondrial genes are
#' present, otherwise "^mt-" (mouse-style), matching the fallback the rest of
#' the package assumes even when neither pattern actually matches (in which
#' case percent.mt will simply be 0 for every cell).
#'
#' @param gene_names Character vector of gene symbols.
#' @return A regex string, either "^MT-" or "^mt-".
#' @noRd
.detect_mt_pattern <- function(gene_names) {
  if (any(grepl("^MT-", gene_names))) "^MT-" else "^mt-"
}

#' Derive a sample name from a .cloupe file path
#'
#' Internal. Strips the directory and the final extension from `cloupe_path`.
#'
#' @param cloupe_path Character. Path to a .cloupe file.
#' @return Character. The basename with its final extension removed.
#' @noRd
.derive_sample_name <- function(cloupe_path) {
  sub("\\.[^.]+$", "", basename(cloupe_path))
}

#' Align a tissue-positions data.frame to a set of cell names
#'
#' Internal. Sets `pos`'s rownames from its `barcode` column, then reindexes
#' it to `cell_names` (the count matrix's column order), so downstream
#' metadata assignment lines up one-to-one with the assay's cells.
#'
#' @param pos A data.frame with a `barcode` column (e.g. read from
#'   tissue_positions.csv).
#' @param cell_names Character vector of cell/barcode names to align to,
#'   typically `colnames(mat)`.
#' @return `pos` reindexed to `cell_names`, with `cell_names` as rownames.
#' @noRd
.align_positions <- function(pos, cell_names) {
  rownames(pos) <- pos$barcode
  pos[cell_names, ]
}

#' Reshape estimate_image_memory()'s raw Python result into a data.frame
#'
#' Internal. `results` is what reticulate hands back from
#' `cloupe_extract$estimate_image_memory()`: a list of named lists, one per
#' requested px_per_bin, each with `px_per_bin` (NULL for native),
#' `image_size` (a 2-element width/height list), `clamped_to_native`,
#' `python_peak_gb`, `r_reload_peak_gb`. Kept separate from
#' `estimate_image_memory()` itself so this reshape step is testable without
#' a real Python call.
#'
#' @param results List of named lists, as returned by
#'   `cloupe_extract$estimate_image_memory()`.
#' @return A data.frame with columns `px_per_bin` (NA = native),
#'   `image_width`, `image_height`, `clamped_to_native`, `python_peak_gb`,
#'   `r_reload_peak_gb`.
#' @noRd
.image_memory_table <- function(results) {
  data.frame(
    px_per_bin = vapply(
      results,
      function(r) if (is.null(r$px_per_bin)) NA_real_ else as.numeric(r$px_per_bin),
      numeric(1)
    ),
    image_width       = vapply(results, function(r) r$image_size[[1]], numeric(1)),
    image_height      = vapply(results, function(r) r$image_size[[2]], numeric(1)),
    clamped_to_native = vapply(results, function(r) r$clamped_to_native, logical(1)),
    python_peak_gb    = vapply(results, function(r) r$python_peak_gb, numeric(1)),
    r_reload_peak_gb  = vapply(results, function(r) r$r_reload_peak_gb, numeric(1))
  )
}
