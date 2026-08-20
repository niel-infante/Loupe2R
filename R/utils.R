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
