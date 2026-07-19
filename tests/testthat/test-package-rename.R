test_that("package metadata and native registration use fastINFLECT", {
  root <- find_package_root()
  desc <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]
  namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  rcpp_r <- readLines(file.path(root, "R", "RcppExports.R"), warn = FALSE)
  rcpp_cpp <- readLines(file.path(root, "src", "RcppExports.cpp"), warn = FALSE)

  expect_identical(unname(desc[["Package"]]), "fastINFLECT")
  expect_true(any(grepl("useDynLib\\(fastINFLECT", namespace)))
  expect_true(any(grepl("_fastINFLECT_inflect_iqr_cpp", rcpp_r, fixed = TRUE)))
  expect_true(any(grepl("R_init_fastINFLECT", rcpp_cpp, fixed = TRUE)))
  expect_false(any(grepl("useDynLib\\(INFLECT", namespace)))
})
