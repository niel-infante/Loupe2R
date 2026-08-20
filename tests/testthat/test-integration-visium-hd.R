# Opt-in: cross-validates cloupe_to_seurat() against the official SpaceRanger
# output for the same sample. Not run by default (no real .cloupe file is
# committed to this repo) -- set LOUPE2R_TEST_DIR to a SpaceRanger Visium HD
# `outs/` directory containing cloupe_008um.cloupe and
# binned_outputs/square_008um/ to run it, e.g.:
#   LOUPE2R_TEST_DIR=/path/to/outs Rscript -e "testthat::test_local()"

test_dir <- Sys.getenv("LOUPE2R_TEST_DIR", unset = "")
skip_if_not(nzchar(test_dir), "Set LOUPE2R_TEST_DIR to run the local Visium HD integration test")

cloupe_file  <- file.path(test_dir, "cloupe_008um.cloupe")
official_dir <- file.path(test_dir, "binned_outputs", "square_008um")
skip_if_not(file.exists(cloupe_file), paste("cloupe_008um.cloupe not found in", test_dir))
skip_if_not(dir.exists(official_dir), paste("binned_outputs/square_008um not found in", test_dir))
skip_if_not_installed("arrow")
skip_if_not_installed("digest")

fixture <- jsonlite::read_json(
  test_path("fixtures", "cloupe_008um_expected.json"),
  simplifyVector = TRUE
)

extract_dir <- withr::local_tempdir()
srt_cloupe <- suppressWarnings(cloupe_to_seurat(cloupe_file, outdir = extract_dir, keep_files = TRUE))

test_that("extraction matches the checked-in regression fixture", {
  # Real .cloupe files are too large to commit; these summary stats were
  # computed once from this same reference file and cross-validated against
  # its paired official SpaceRanger output (see the vignette). If a future
  # SpaceRanger/Loupe Browser release changes the format in a way this
  # parser silently misreads, this is far more likely to catch it than
  # eyeballing a plot.
  expect_equal(nrow(srt_cloupe), fixture$n_features)
  expect_equal(ncol(srt_cloupe), fixture$n_barcodes)
  expect_equal(sum(srt_cloupe$nCount_Spatial), fixture$total_umi)
  expect_equal(Seurat::Misc(srt_cloupe, "bin_size_um"), fixture$bin_size_um)
  expect_equal(
    srt_cloupe[["slice1"]]@scale.factors$spot,
    fixture$spot_diameter_fullres,
    tolerance = 1e-5
  )
  expect_equal(
    digest::digest(sort(colnames(srt_cloupe)), algo = "sha256"),
    fixture$barcode_set_sha256
  )
})

test_that("no unvalidated .cloupe format version is reported for this file", {
  fmt <- Seurat::Misc(srt_cloupe, "cloupe_format_info")
  expect_equal(fmt$container, fixture$cloupe_format_info$container)
  expect_equal(fmt$matrix, fixture$cloupe_format_info$matrix)
})

test_that("cloupe_to_seurat agrees with official SpaceRanger output", {
  srt_official  <- Seurat::Load10X_Spatial(test_dir, bin.size = 8)
  raw_official  <- Seurat::Read10X_h5(file.path(official_dir, "raw_feature_bc_matrix.h5"))

  cloupe_bc   <- colnames(srt_cloupe)
  filtered_bc <- colnames(srt_official)
  raw_bc      <- colnames(raw_official)

  overlap_filtered <- length(intersect(cloupe_bc, filtered_bc)) / length(cloupe_bc)
  overlap_raw      <- length(intersect(cloupe_bc, raw_bc)) / length(cloupe_bc)
  message(sprintf(
    "Barcode overlap vs filtered official: %.4f, vs raw official: %.4f",
    overlap_filtered, overlap_raw
  ))
  expect_gt(max(overlap_filtered, overlap_raw), 0.95)

  shared <- intersect(cloupe_bc, filtered_bc)
  expect_gt(length(shared), 0)

  nc_cloupe   <- srt_cloupe$nCount_Spatial[shared]
  nc_official <- srt_official[[paste0("nCount_", Seurat::DefaultAssay(srt_official))]][shared, 1]
  expect_gt(cor(nc_cloupe, nc_official), 0.99)

  pos_official <- as.data.frame(arrow::read_parquet(
    file.path(official_dir, "spatial", "tissue_positions.parquet")
  ))
  rownames(pos_official) <- pos_official$barcode
  pos_official <- pos_official[shared, ]

  pos_cloupe <- read.csv(file.path(extract_dir, "tissue_positions.csv"), stringsAsFactors = FALSE)
  rownames(pos_cloupe) <- pos_cloupe$barcode
  pos_cloupe <- pos_cloupe[shared, ]

  row_diff <- abs(pos_cloupe$array_row - pos_official$array_row)
  col_diff <- abs(pos_cloupe$array_col - pos_official$array_col)
  message(sprintf(
    "array_row exact match: %.4f (max diff %d), array_col exact match: %.4f (max diff %d)",
    mean(row_diff == 0), max(row_diff), mean(col_diff == 0), max(col_diff)
  ))
  # array_row/array_col are parsed directly from Visium HD barcode names
  # (e.g. "s_008um_00269_00526-1"), which is exact, not an approximation --
  # see the comment above parse_array_position() in cloupe_extract.py.
  expect_gt(mean(row_diff == 0), 0.99)
  expect_gt(mean(col_diff == 0), 0.99)
  expect_true(all(row_diff <= 2))
  expect_true(all(col_diff <= 2))

  cor_row <- cor(pos_cloupe$pxl_row_in_fullres, pos_official$pxl_row_in_fullres)
  cor_col <- cor(pos_cloupe$pxl_col_in_fullres, pos_official$pxl_col_in_fullres)
  cor_row_vs_col <- cor(pos_cloupe$pxl_row_in_fullres, pos_official$pxl_col_in_fullres)
  message(sprintf("pxl_row correlation: %.6f, pxl_col correlation: %.6f", cor_row, cor_col))
  expect_gt(cor_row, 0.999)
  expect_gt(cor_col, 0.999)
  # sanity check against an axis swap: row shouldn't correlate with col
  expect_lt(abs(cor_row_vs_col), 0.5)
})
