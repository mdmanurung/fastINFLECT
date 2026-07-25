test_that("release sources consistently identify fastINFLECT 1.0.0", {
  root <- find_package_root()
  desc <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]

  expect_identical(unname(desc[["Version"]]), "1.0.0")
  expect_true(file.exists(file.path(
    root,
    "vignettes",
    "using-fastINFLECT-1.Rmd"
  )))
  expect_false(file.exists(file.path(
    root,
    "vignettes",
    "migrating-to-fastINFLECT-2.Rmd"
  )))

  source_paths <- file.path(root, c(
    "README.md",
    "NEWS.md",
    "R/INFLECT.R",
    "R/INFLECT3-package.R",
    "R/QC-to-curve.R",
    "R/inflect-provenance.R",
    "R/iteration-QC.R",
    "R/plot-markerperformance.R",
    "vignettes/benchmark-inflect-vs-consensus.Rmd",
    "vignettes/comparing-metaclustering.Rmd",
    "vignettes/using-fastINFLECT-1.Rmd"
  ))
  source_paths <- source_paths[file.exists(source_paths)]
  release_text <- paste(
    unlist(lapply(source_paths, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  stale_patterns <- c(
    "fastINFLECT 2.0",
    "Version 2.0",
    "version 2.0",
    "pre-2.0",
    "2.0 migration",
    "migrating-to-fastINFLECT-2"
  )
  for (pattern in stale_patterns) {
    expect_false(grepl(pattern, release_text, fixed = TRUE), info = pattern)
  }

  pkgdown_path <- file.path(root, "_pkgdown.yml")
  if (file.exists(pkgdown_path)) {
    pkgdown <- paste(readLines(pkgdown_path, warn = FALSE), collapse = "\n")
    expect_true(grepl("using-fastINFLECT-1", pkgdown, fixed = TRUE))
    expect_false(grepl("migrating-to-fastINFLECT-2", pkgdown, fixed = TRUE))
  }
})
