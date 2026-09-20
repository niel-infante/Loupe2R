# Opt-in: real-data regression guard for the R-side crash that motivated the
# px_per_bin fix -- png::readPNG() expanding a native-resolution tissue image
# to double precision and hitting R's 30 GB vector memory limit. Not run by
# default (the ~7.6 GB .cloupe file isn't committed to this repo) -- set
# LOUPE2R_LARGE_IMAGE_TEST_DIR to a directory containing cloupe_008um.cloupe
# to run it, e.g.:
#   LOUPE2R_LARGE_IMAGE_TEST_DIR=/path/to/crc_dataset Rscript -e "testthat::test_local()"

test_dir <- Sys.getenv("LOUPE2R_LARGE_IMAGE_TEST_DIR", unset = "")
skip_if_not(nzchar(test_dir), "Set LOUPE2R_LARGE_IMAGE_TEST_DIR to run the large-image integration test")

cloupe_file <- file.path(test_dir, "cloupe_008um.cloupe")
skip_if_not(file.exists(cloupe_file), paste("cloupe_008um.cloupe not found in", test_dir))

test_that("cloupe_to_seurat with default px_per_bin completes and embeds a small image", {
  # Default px_per_bin -- caller passes nothing. This is the direct
  # regression test for the R-side png::readPNG() crash already reproduced
  # against this exact file at native resolution: proves the new safe
  # default alone prevents it, with zero caller changes required.
  srt <- suppressWarnings(cloupe_to_seurat(cloupe_file))

  img <- srt[["slice1"]]@image
  expect_lt(dim(img)[1] * dim(img)[2], 100e6)
  # Native resolution for this sample is 75250 x 48740 -- confirm this is
  # genuinely smaller, not a coincidentally-small crop.
  expect_lt(dim(img)[1], 48740)
  expect_lt(dim(img)[2], 75250)

  scale_facs <- srt[["slice1"]]@scale.factors
  expect_lt(scale_facs$hires, 1.0)
  expect_equal(scale_facs$hires, scale_facs$lowres)
})
