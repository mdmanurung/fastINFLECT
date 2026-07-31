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

test_that("release metadata points at fastINFLECT", {
  root <- find_package_root()
  desc <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]

  expect_true(grepl("github.com/mdmanurung/fastINFLECT", desc[["URL"]], fixed = TRUE))
  expect_identical(unname(desc[["BugReports"]]), "https://github.com/mdmanurung/fastINFLECT/issues")
})

test_that("generated documentation is complete and route-clean in a checkout", {
  root <- find_package_root()
  if (!file.exists(file.path(root, ".git"))) {
    skip("Generated pkgdown files are excluded from the source tarball by .Rbuildignore.")
  }

  required_relative <- c(
    "_pkgdown.yml",
    "vignettes/fastINFLECT.html",
    "docs/articles/fastINFLECT.html",
    "docs/articles/fastINFLECT.md",
    "docs/articles/fastINFLECT_files/figure-html/marker-performance-1.png",
    "docs/articles/fastINFLECT_files/figure-html/overview-1.png",
    "docs/articles/index.html",
    "docs/articles/index.md",
    "docs/index.html",
    "docs/index.md",
    "docs/search.json",
    "docs/sitemap.xml",
    "docs/llms.txt",
    "docs/pkgdown.yml"
  )
  required_paths <- file.path(root, required_relative)
  missing_relative <- required_relative[!file.exists(required_paths)]
  expect_true(
    length(missing_relative) == 0L,
    info = paste("Missing required generated files:", paste(missing_relative, collapse = ", "))
  )
  if (length(missing_relative) != 0L) {
    return(invisible(NULL))
  }

  retired_slugs <- c(
    "using-fastINFLECT-1",
    "benchmark-inflect-vs-consensus",
    "comparing-metaclustering",
    "migrating-to-fastINFLECT-2"
  )
  retired_relative <- unlist(lapply(retired_slugs, function(slug) {
    c(
      file.path("vignettes", paste0(slug, c(".Rmd", ".html", "_files"))),
      file.path("docs", "articles", paste0(slug, c(".html", ".md", "_files")))
    )
  }), use.names = FALSE)
  unexpected_retired <- retired_relative[file.exists(file.path(root, retired_relative))]
  expect_true(
    length(unexpected_retired) == 0L,
    info = paste("Retired documentation paths remain:", paste(unexpected_retired, collapse = ", "))
  )

  route_relative <- c(
    "_pkgdown.yml",
    "docs/articles/fastINFLECT.html",
    "docs/articles/fastINFLECT.md",
    "docs/articles/index.html",
    "docs/articles/index.md",
    "docs/index.html",
    "docs/index.md",
    "docs/search.json",
    "docs/sitemap.xml",
    "docs/llms.txt",
    "docs/pkgdown.yml"
  )
  route_text <- setNames(lapply(route_relative, function(path) {
    paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
  }), route_relative)
  for (slug in retired_slugs) {
    expect_false(
      any(vapply(route_text, grepl, logical(1), pattern = slug, fixed = TRUE)),
      info = slug
    )
  }

  pkgdown <- readLines(file.path(root, "_pkgdown.yml"), warn = FALSE)
  get_started_line <- grep(
    "^[[:space:]]*text:[[:space:]]*Get started[[:space:]]*$",
    pkgdown
  )
  expect_length(get_started_line, 1L)
  if (length(get_started_line) == 1L) {
    expect_identical(
      trimws(pkgdown[[get_started_line + 1L]]),
      "href: articles/fastINFLECT.html"
    )
  }

  expect_true(grepl(
    'href="fastINFLECT.html"',
    route_text[["docs/articles/index.html"]],
    fixed = TRUE
  ))
  expect_true(grepl(
    'href="articles/fastINFLECT.html"',
    route_text[["docs/index.html"]],
    fixed = TRUE
  ))
  canonical_url <- "https://mdmanurung.github.io/fastINFLECT/articles/fastINFLECT.html"
  expect_true(grepl(
    paste0('"path":"', canonical_url, '"'),
    route_text[["docs/search.json"]],
    fixed = TRUE
  ))
  expect_true(grepl(
    paste0("<loc>", canonical_url, "</loc>"),
    route_text[["docs/sitemap.xml"]],
    fixed = TRUE
  ))

  rendered_relative <- c(
    "vignettes/fastINFLECT.html",
    "docs/articles/fastINFLECT.html",
    "docs/articles/fastINFLECT.md"
  )
  rendered_text <- lapply(rendered_relative, function(path) {
    paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
  })
  stale <- unlist(rendered_text, use.names = FALSE)
  expect_false(any(grepl("INFLECT_0.2.1", stale, fixed = TRUE)))

  search_lines <- readLines(file.path(root, "docs", "search.json"), warn = FALSE)
  expect_length(search_lines, 1L)
  search_json <- paste(search_lines, collapse = "")
  expect_true(startsWith(search_json, "[{"))
  expect_true(endsWith(search_json, "}]"))
  expect_false(grepl('"path":[]', search_json, fixed = TRUE))
  expect_false(grepl('"path":null', search_json, fixed = TRUE))
  search_records <- strsplit(
    substring(search_json, 2L, nchar(search_json) - 1L),
    "},{",
    fixed = TRUE
  )[[1L]]
  expect_gt(length(search_records), 0L)
  expect_true(all(grepl(
    '^\\{?"path":"[^"[:space:]][^"]*"',
    search_records,
    perl = TRUE
  )))
})
