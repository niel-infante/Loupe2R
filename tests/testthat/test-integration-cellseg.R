# Opt-in: cross-validates cloupe_to_seurat() against official SpaceRanger
# cell-segmentation output for the same sample. Not run by default (no real
# .cloupe file is committed to this repo) -- set LOUPE2R_CELLSEG_TEST_DIR to
# a directory containing cloupe_cell.cloupe and segmented_outputs/
# (filtered_feature_cell_matrix/, cell_segmentations.geojson) to run it, e.g.:
#   LOUPE2R_CELLSEG_TEST_DIR=/path/to/dir Rscript -e "testthat::test_local()"

test_dir <- Sys.getenv("LOUPE2R_CELLSEG_TEST_DIR", unset = "")
skip_if_not(nzchar(test_dir), "Set LOUPE2R_CELLSEG_TEST_DIR to run the local cell-segmentation integration test")

cloupe_file  <- file.path(test_dir, "cloupe_cell.cloupe")
official_dir <- file.path(test_dir, "segmented_outputs")
skip_if_not(file.exists(cloupe_file), paste("cloupe_cell.cloupe not found in", test_dir))
skip_if_not(dir.exists(official_dir), paste("segmented_outputs not found in", test_dir))
skip_if_not_installed("sf")

test_that("cloupe_to_seurat agrees with official SpaceRanger cell-segmentation output", {
  extract_dir <- withr::local_tempdir()
  srt <- suppressWarnings(cloupe_to_seurat(
    cloupe_file, outdir = extract_dir, keep_files = TRUE, include_image = FALSE
  ))

  cloupe_bc <- colnames(srt)

  # ---- Count concordance vs. official filtered_feature_cell_matrix ----
  off_mat <- Seurat::ReadMtx(
    mtx            = file.path(official_dir, "filtered_feature_cell_matrix", "matrix.mtx.gz"),
    cells          = file.path(official_dir, "filtered_feature_cell_matrix", "barcodes.tsv.gz"),
    features       = file.path(official_dir, "filtered_feature_cell_matrix", "features.tsv.gz"),
    feature.column = 2,
    cell.column    = 1,
    skip.cell      = 0,
    skip.feature   = 0
  )
  official_bc <- colnames(off_mat)

  overlap_counts <- length(intersect(cloupe_bc, official_bc)) / length(cloupe_bc)
  message(sprintf(
    "Barcode overlap vs official filtered_feature_cell_matrix: %.4f%% (%d/%d)",
    overlap_counts * 100, length(intersect(cloupe_bc, official_bc)), length(cloupe_bc)
  ))
  expect_gt(overlap_counts, 0.99)

  shared_counts <- intersect(cloupe_bc, official_bc)
  nc_cloupe   <- srt$nCount_Spatial[shared_counts]
  nc_official <- Matrix::colSums(off_mat[, shared_counts, drop = FALSE])
  message(sprintf(
    "Count exact match: %.4f%%   correlation: %.6f",
    mean(nc_cloupe == nc_official) * 100, cor(nc_cloupe, nc_official)
  ))
  expect_gt(mean(nc_cloupe == nc_official), 0.99)
  expect_gt(cor(nc_cloupe, nc_official), 0.999)

  # ---- Position concordance vs. official cell_segmentations.geojson ----
  # True area centroid of each cell's official polygon boundary -- same
  # definition (and, for cells with GeoJSON coverage, the same underlying
  # polygon) our own extraction now reads directly. Uses sf::st_centroid()
  # rather than a hand-rolled shoelace formula; both compute the same
  # geometric quantity, this just avoids re-deriving it in R.
  geo <- sf::st_read(file.path(official_dir, "cell_segmentations.geojson"), quiet = TRUE)
  geo <- geo[sf::st_geometry_type(geo) == "POLYGON", ]

  # These are planar pixel coordinates, not geographic lon/lat -- GeoJSON
  # carries no CRS, and sf defaults to treating that as longlat and routing
  # centroid computation through spherical (s2) geometry, which is both the
  # wrong geometric model here and chokes on real-image polygon rings that
  # are perfectly valid as plane polygons. Force planar (GEOS) computation.
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)

  centroids <- sf::st_coordinates(sf::st_centroid(geo))

  official_pos <- data.frame(
    barcode    = sprintf("cellid_%09d-1", as.integer(geo$cell_id)),
    official_x = centroids[, "X"],
    official_y = centroids[, "Y"]
  )
  rownames(official_pos) <- official_pos$barcode

  shared_pos <- intersect(cloupe_bc, official_pos$barcode)
  message(sprintf(
    "Barcode overlap vs official cell_segmentations.geojson: %.4f%% (%d/%d)",
    length(shared_pos) / length(cloupe_bc) * 100, length(shared_pos), length(cloupe_bc)
  ))
  expect_gt(length(shared_pos) / length(cloupe_bc), 0.99)

  cx <- srt$pxl_col_in_fullres[shared_pos]
  cy <- srt$pxl_row_in_fullres[shared_pos]
  ox <- official_pos[shared_pos, "official_x"]
  oy <- official_pos[shared_pos, "official_y"]
  dist <- sqrt((cx - ox)^2 + (cy - oy)^2)
  message(sprintf(
    "Centroid distance (px): mean=%.4f  median=%.4f  max=%.4f",
    mean(dist), stats::median(dist), max(dist)
  ))
  expect_lt(mean(dist), 0.01)
  expect_lt(max(dist), 1)
})
