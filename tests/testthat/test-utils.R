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
