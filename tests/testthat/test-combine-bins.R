# Synthetic single-resolution Seurat objects shaped like cloupe_to_seurat()
# output, without needing any real .cloupe file or Python dependency.
.make_synthetic_bin_obj <- function(n_cells, bin_um, seed = 1) {
  set.seed(seed)
  mat <- matrix(rpois(20 * n_cells, 3), nrow = 20, ncol = n_cells)
  rownames(mat) <- paste0("gene", seq_len(20))
  colnames(mat) <- paste0("bc", seq_len(n_cells))
  srt <- suppressWarnings(Seurat::CreateSeuratObject(counts = mat, assay = "Spatial"))

  img_arr <- array(0.5, dim = c(10, 10, 3))
  coords_df <- data.frame(
    tissue = 1L,
    row = seq_len(n_cells) %% 5,
    col = seq_len(n_cells) %/% 5,
    imagerow = stats::runif(n_cells, 0, 10),
    imagecol = stats::runif(n_cells, 0, 10),
    row.names = colnames(mat)
  )
  sfacs <- Seurat::scalefactors(spot = 1, fiducial = 1, hires = 1, lowres = 1)
  img <- tryCatch(
    methods::new("VisiumV2", image = img_arr, scale.factors = sfacs,
                 coordinates = coords_df, assay = "Spatial", key = "slice1_"),
    error = function(e)
      methods::new("VisiumV1", image = img_arr, scale.factors = sfacs,
                   coordinates = coords_df, assay = "Spatial", key = "slice1_")
  )
  srt[["slice1"]] <- img
  Seurat::Misc(srt, "bin_size_um") <- bin_um
  srt
}

test_that("combine_cloupe_bins merges two resolutions into one object", {
  s8  <- .make_synthetic_bin_obj(40, 8, seed = 1)
  s16 <- .make_synthetic_bin_obj(10, 16, seed = 2)

  srt <- combine_cloupe_bins(s8, s16)

  expect_setequal(Seurat::Assays(srt), c("Spatial.008um", "Spatial.016um"))
  expect_setequal(Seurat::Images(srt), c("slice1.008um", "slice1.016um"))
  expect_equal(Seurat::DefaultAssay(srt), "Spatial.008um")
  expect_equal(ncol(srt[["Spatial.008um"]]), 40)
  expect_equal(ncol(srt[["Spatial.016um"]]), 10)
  # object-level cell identity tracks the base (first) object
  expect_equal(ncol(srt), 40)
})

test_that("combine_cloupe_bins supports three or more resolutions", {
  s2  <- .make_synthetic_bin_obj(160, 2, seed = 3)
  s8  <- .make_synthetic_bin_obj(10, 8, seed = 4)
  s16 <- .make_synthetic_bin_obj(4, 16, seed = 5)

  srt <- combine_cloupe_bins(s2, s8, s16)
  expect_setequal(Seurat::Assays(srt), c("Spatial.002um", "Spatial.008um", "Spatial.016um"))
})

test_that("combine_cloupe_bins requires at least two objects", {
  s8 <- .make_synthetic_bin_obj(10, 8)
  expect_error(combine_cloupe_bins(s8), "at least two Seurat objects")
})

test_that("combine_cloupe_bins requires Seurat objects", {
  s8 <- .make_synthetic_bin_obj(10, 8)
  expect_error(combine_cloupe_bins(s8, "not a seurat object"), "must be Seurat objects")
})

test_that("combine_cloupe_bins requires Misc bin_size_um on every object", {
  s8 <- .make_synthetic_bin_obj(10, 8)
  s_no_bin <- .make_synthetic_bin_obj(10, 16)
  Seurat::Misc(s_no_bin, "bin_size_um") <- NULL
  expect_error(combine_cloupe_bins(s8, s_no_bin), "bin_size_um")
})

test_that("combine_cloupe_bins errors on duplicate bin sizes", {
  s8a <- .make_synthetic_bin_obj(40, 8, seed = 1)
  s8b <- .make_synthetic_bin_obj(40, 8, seed = 6)
  expect_error(combine_cloupe_bins(s8a, s8b), "same bin size")
})

test_that("combine_cloupe_bins warns when bin counts look inconsistent", {
  s8    <- .make_synthetic_bin_obj(400, 8, seed = 1)
  s16_bad <- .make_synthetic_bin_obj(5, 16, seed = 7)  # way too few for the same tissue
  expect_warning(combine_cloupe_bins(s8, s16_bad), "inconsistent for the same tissue area")
})

test_that("combine_cloupe_bins does not warn about bin-count consistency when counts are fine", {
  # Note: Seurat/SeuratObject's own assignment machinery can still emit its
  # own benign warnings here (e.g. "Key 'spatial_' taken, using
  # 'spatial008um_' instead", "Adding image data that isn't associated with
  # any assays") as a side effect of the remove/reassign rename pattern --
  # those aren't ours to silence and aren't what this test is checking.
  s8  <- .make_synthetic_bin_obj(400, 8, seed = 1)
  s16 <- .make_synthetic_bin_obj(100, 16, seed = 2)  # ~4x fewer, matches (16/8)^2
  warnings_caught <- character(0)
  withCallingHandlers(
    combine_cloupe_bins(s8, s16),
    warning = function(w) {
      warnings_caught <<- c(warnings_caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("inconsistent for the same tissue area", warnings_caught)))
})
