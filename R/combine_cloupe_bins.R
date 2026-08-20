#' Combine single-resolution cloupe imports into one multi-resolution object
#'
#' A single \code{.cloupe} file only ever captures the one bin resolution
#' that was open in Loupe Browser when it was saved. If you have separate
#' \code{.cloupe} files for the same Visium HD capture area at different bin
#' sizes (e.g. 8\eqn{\mu}{u}m and 16\eqn{\mu}{u}m), import each one
#' separately with \code{\link{cloupe_to_seurat}()} and combine the results
#' with this function into a single Seurat object with one assay and one
#' image per resolution -- matching the naming convention
#' \code{Seurat::Load10X_Spatial(bin.size = c(8, 16))} uses
#' (\code{"Spatial.008um"}, \code{"Spatial.016um"}, ...).
#'
#' Because Seurat v5 assays are not required to share an identical set of
#' cells, this is a supported structure, not a workaround -- it is the same
#' mechanism Seurat's own official Visium HD loader relies on. One
#' consequence carries over here too: the combined object's own cell
#' identity -- \code{Cells()}/\code{ncol()}, \code{meta.data}, idents --
#' stays pinned to the *base* object's cells (whichever object is passed
#' first) even after switching \code{DefaultAssay()} to a different
#' resolution; only \code{ncol(object[[assay]])} reflects that assay's own
#' cell count. Index into a specific assay/image explicitly for
#' resolution-specific work, exactly as in Seurat's own Visium HD workflows.
#'
#' Each input's bin size is read from \code{Misc(x, "bin_size_um")}, which
#' \code{\link{cloupe_to_seurat}()} sets automatically -- it is not guessed
#' from a filename.
#'
#' @param ... Two or more Seurat objects, each returned by
#'   \code{\link{cloupe_to_seurat}()} on a \code{.cloupe} file for a
#'   different bin resolution of the *same* tissue section. The first object
#'   determines the combined object's default assay/image and its
#'   \code{meta.data} baseline; 10x recommends 8\eqn{\mu}{u}m bins for
#'   general analysis, so pass that one first if in doubt.
#'
#' @return A single Seurat object with one assay/image per input object,
#'   named \code{"<assay>.<NNN>um"} / \code{"<slice>.<NNN>um"}.
#'
#' @examples
#' \dontrun{
#' s008 <- cloupe_to_seurat("sample_008um.cloupe")
#' s016 <- cloupe_to_seurat("sample_016um.cloupe")
#' srt  <- combine_cloupe_bins(s008, s016)
#' Assays(srt)  # "Spatial.008um" "Spatial.016um"
#' }
#'
#' @export
combine_cloupe_bins <- function(...) {
  objs <- list(...)
  if (length(objs) < 2)
    stop("combine_cloupe_bins() needs at least two Seurat objects (one per bin resolution).")
  if (!all(vapply(objs, inherits, logical(1), what = "Seurat")))
    stop("All arguments to combine_cloupe_bins() must be Seurat objects.")

  bin_um <- vapply(objs, function(x) {
    b <- Seurat::Misc(x, "bin_size_um")
    if (is.null(b) || !is.numeric(b) || length(b) != 1)
      stop(
        "Every object passed to combine_cloupe_bins() must carry ",
        "Misc(x, \"bin_size_um\"), as set automatically by cloupe_to_seurat()."
      )
    b
  }, numeric(1))

  if (anyDuplicated(round(bin_um, 3)))
    stop(
      "Two or more objects resolve to the same bin size (~",
      paste(unique(bin_um[duplicated(round(bin_um, 3))]), collapse = ", "),
      " um) -- combine_cloupe_bins() expects one object per distinct resolution."
    )

  .warn_if_bin_counts_inconsistent(objs, bin_um)

  suffix <- function(um) sprintf("%03dum", round(um))

  # Base object: rename its default assay + image(s) in place.
  base            <- objs[[1]]
  base_suffix     <- suffix(bin_um[1])
  base_assay_orig <- Seurat::DefaultAssay(base)
  base_assay_new  <- paste0(base_assay_orig, ".", base_suffix)

  # RenameAssays() (not re-exported by Seurat itself, hence SeuratObject::)
  # handles reassigning the default-assay pointer atomically; a manual
  # remove/reassign errors with "Cannot delete default assay".
  base <- SeuratObject::RenameAssays(
    base,
    assay.name     = base_assay_orig,
    new.assay.name = base_assay_new
  )

  for (img_name in Seurat::Images(base)) {
    img_obj <- base[[img_name]]
    base[[img_name]] <- NULL
    base[[paste0(img_name, ".", base_suffix)]] <- img_obj
  }

  Seurat::DefaultAssay(base) <- base_assay_new
  # No explicit "default image" to update here either (see the comment in
  # cloupe_to_seurat.R) -- Images(object, assay = DefaultAssay(object)) is
  # how modern Seurat/SeuratObject (>= 5.4) resolves this automatically.

  # Additional objects: attach their assay + image(s) into base under the
  # same "<name>.<NNN>um" convention, without touching base's meta.data.
  for (i in seq_along(objs)[-1]) {
    obj <- objs[[i]]
    suf <- suffix(bin_um[i])

    add_assay_name <- paste0(Seurat::DefaultAssay(obj), ".", suf)
    base[[add_assay_name]] <- obj[[Seurat::DefaultAssay(obj)]]

    for (img_name in Seurat::Images(obj)) {
      base[[paste0(img_name, ".", suf)]] <- obj[[img_name]]
    }
  }

  base
}

#' Warn if bin counts across resolutions look inconsistent for one tissue area
#'
#' Internal. Nothing in a single .cloupe file identifies which capture area
#' it came from, so this is the best available guard against accidentally
#' combining resolutions from different samples: for the same tissue area,
#' bin counts should scale roughly with (coarse_bin/fine_bin)^2.
#'
#' @param objs List of Seurat objects.
#' @param bin_um Numeric vector of bin sizes (microns), same length as objs.
#' @return Invisible NULL; called for its warning() side effect.
#' @noRd
.warn_if_bin_counts_inconsistent <- function(objs, bin_um) {
  ncells <- vapply(objs, ncol, numeric(1))
  ord    <- order(bin_um)
  for (i in seq_len(length(ord) - 1)) {
    a <- ord[i]; b <- ord[i + 1]
    expected_ratio <- (bin_um[b] / bin_um[a])^2
    observed_ratio <- ncells[a] / ncells[b]
    if (is.finite(expected_ratio) &&
        (observed_ratio < expected_ratio / 4 || observed_ratio > expected_ratio * 4)) {
      warning(
        sprintf(
          paste(
            "Bin counts for %gum (%d bins) and %gum (%d bins) look",
            "inconsistent for the same tissue area (expected roughly %.1fx",
            "more %gum bins; observed %.1fx). Confirm these .cloupe files",
            "are the same capture area."
          ),
          bin_um[a], ncells[a], bin_um[b], ncells[b],
          expected_ratio, bin_um[a], observed_ratio
        ),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}
