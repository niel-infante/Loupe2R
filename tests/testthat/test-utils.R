test_that(".detect_mt_pattern finds human mitochondrial genes", {
  expect_equal(.detect_mt_pattern(c("MT-ND1", "GAPDH", "MT-CO1")), "^MT-")
})

test_that(".detect_mt_pattern finds mouse mitochondrial genes", {
  expect_equal(.detect_mt_pattern(c("mt-Nd1", "Gapdh", "mt-Co1")), "^mt-")
})

test_that(".detect_mt_pattern falls back to mouse pattern when neither matches", {
  # Matches the existing behaviour of the code this was extracted from:
  # no MT- genes -> falls through to "^mt-", even if that also matches
  # nothing (percent.mt ends up 0 either way).
  expect_equal(.detect_mt_pattern(c("GAPDH", "ACTB")), "^mt-")
})

test_that(".detect_mt_pattern prefers human when both conventions are present", {
  expect_equal(.detect_mt_pattern(c("MT-ND1", "mt-Nd1")), "^MT-")
})

test_that(".derive_sample_name strips directory and final extension", {
  expect_equal(.derive_sample_name("/path/to/sample1.cloupe"), "sample1")
})

test_that(".derive_sample_name only strips the final extension", {
  expect_equal(.derive_sample_name("sample.v2.cloupe"), "sample.v2")
})

test_that(".derive_sample_name handles a path with no extension", {
  expect_equal(.derive_sample_name("sample_no_ext"), "sample_no_ext")
})

test_that(".align_positions reindexes to the requested cell order", {
  pos <- data.frame(
    barcode = c("bc2", "bc1", "bc3"),
    value   = c(20, 10, 30),
    stringsAsFactors = FALSE
  )
  aligned <- .align_positions(pos, c("bc1", "bc2", "bc3"))
  expect_equal(rownames(aligned), c("bc1", "bc2", "bc3"))
  expect_equal(aligned$value, c(10, 20, 30))
})

test_that(".align_positions produces NA rows for cell names absent from pos", {
  pos <- data.frame(
    barcode = c("bc1", "bc2"),
    value   = c(10, 20),
    stringsAsFactors = FALSE
  )
  aligned <- .align_positions(pos, c("bc1", "bc_missing"))
  expect_equal(aligned$value, c(10, NA))
})

# .image_memory_table -- fake results shaped exactly as reticulate hands
# back cloupe_extract$estimate_image_memory()'s Python list-of-dicts: NULL
# for native px_per_bin, image_size as a 2-element list (Python tuples
# convert to R lists, not vectors).
.fake_image_memory_result <- function(px_per_bin, w, h, clamped, py_gb, r_gb) {
  list(
    px_per_bin = px_per_bin,
    image_size = list(w, h),
    clamped_to_native = clamped,
    python_peak_gb = py_gb,
    r_reload_peak_gb = r_gb
  )
}

test_that(".image_memory_table converts NULL px_per_bin to NA", {
  results <- list(
    .fake_image_memory_result(NULL, 75250, 48740, FALSE, 22.01, 88.02),
    .fake_image_memory_result(4, 10298, 6670, FALSE, 1.38, 1.65)
  )
  tbl <- .image_memory_table(results)
  expect_equal(tbl$px_per_bin, c(NA_real_, 4))
})

test_that(".image_memory_table unpacks image_size into width/height columns", {
  results <- list(.fake_image_memory_result(4, 10298, 6670, FALSE, 1.38, 1.65))
  tbl <- .image_memory_table(results)
  expect_equal(tbl$image_width, 10298)
  expect_equal(tbl$image_height, 6670)
})

test_that(".image_memory_table preserves row order and all columns", {
  results <- list(
    .fake_image_memory_result(NULL, 100, 100, FALSE, 1, 2),
    .fake_image_memory_result(16, 90, 90, FALSE, 0.9, 1.8),
    .fake_image_memory_result(1, 10, 10, TRUE, 0.01, 0.02)
  )
  tbl <- .image_memory_table(results)
  expect_equal(nrow(tbl), 3)
  expect_equal(tbl$px_per_bin, c(NA_real_, 16, 1))
  expect_equal(tbl$clamped_to_native, c(FALSE, FALSE, TRUE))
  expect_named(tbl, c("px_per_bin", "image_width", "image_height",
                       "clamped_to_native", "python_peak_gb", "r_reload_peak_gb"))
})
