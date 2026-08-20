test_that("set_cloupe_path updates the resolved path", {
  old <- .pkg_env$cloupe_pkg
  on.exit(.pkg_env$cloupe_pkg <- old, add = TRUE)

  set_cloupe_path("/some/custom/path")
  expect_equal(.pkg_env$cloupe_pkg, "/some/custom/path")
})

test_that(".resolve_cloupe_pkg errors when the directory does not exist", {
  old <- .pkg_env$cloupe_pkg
  on.exit(.pkg_env$cloupe_pkg <- old, add = TRUE)

  set_cloupe_path(file.path(tempdir(), "definitely_does_not_exist_12345"))
  expect_error(.resolve_cloupe_pkg(), "Cannot find the cellgeni/cloupe")
})

test_that(".resolve_cloupe_pkg errors when the marker file is missing", {
  old <- .pkg_env$cloupe_pkg
  on.exit(.pkg_env$cloupe_pkg <- old, add = TRUE)

  bad_dir <- withr::local_tempdir()
  # directory exists but has no cloupe/__init__.py inside it
  set_cloupe_path(bad_dir)
  expect_error(.resolve_cloupe_pkg(), "Cannot find the cellgeni/cloupe")
})

test_that(".resolve_cloupe_pkg succeeds when the marker file is present", {
  old <- .pkg_env$cloupe_pkg
  on.exit(.pkg_env$cloupe_pkg <- old, add = TRUE)

  good_dir <- withr::local_tempdir()
  dir.create(file.path(good_dir, "cloupe"))
  file.create(file.path(good_dir, "cloupe", "__init__.py"))

  set_cloupe_path(good_dir)
  resolved <- .resolve_cloupe_pkg()
  expect_equal(normalizePath(resolved), normalizePath(good_dir))
})
