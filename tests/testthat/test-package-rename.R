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

test_that("release metadata and rendered docs point at fastINFLECT", {
  root <- find_package_root()
  desc <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]

  expect_true(grepl("github.com/mdmanurung/fastINFLECT", desc[["URL"]], fixed = TRUE))
  expect_identical(unname(desc[["BugReports"]]), "https://github.com/mdmanurung/fastINFLECT/issues")
  pkgdown_path <- file.path(root, "_pkgdown.yml")
  if (file.exists(pkgdown_path)) {
    pkgdown <- readLines(pkgdown_path, warn = FALSE)
    expect_true(any(grepl("mdmanurung.github.io/fastINFLECT", pkgdown, fixed = TRUE)))
  }

  docs_to_check <- file.path(root, c(
    "vignettes/comparing-metaclustering.html",
    "docs/articles/comparing-metaclustering.html",
    "docs/search.json"
  ))
  docs_to_check <- docs_to_check[file.exists(docs_to_check)]
  if (length(docs_to_check)) {
    stale <- unlist(lapply(docs_to_check, readLines, warn = FALSE), use.names = FALSE)
    expect_false(any(grepl("INFLECT_0.2.1", stale, fixed = TRUE)))
  }

  search_path <- file.path(root, "docs", "search.json")
  if (file.exists(search_path)) {
    search_json <- paste(readLines(search_path, warn = FALSE), collapse = "\n")
    expect_false(grepl('"path":\\[\\]', search_json))
    expect_true(grepl('"path":"https://mdmanurung.github.io/fastINFLECT/', search_json))
  }
})
