test_that("release sources consistently identify fastINFLECT 1.0.0", {
  root <- find_package_root()
  source_paths <- file.path(root, c(
    "DESCRIPTION",
    "README.md",
    "NEWS.md",
    "R/INFLECT.R",
    "R/INFLECT3-package.R",
    "R/QC-to-curve.R",
    "R/inflect-provenance.R",
    "R/iteration-QC.R",
    "R/plot-markerperformance.R",
    "vignettes/fastINFLECT.Rmd"
  ))
  missing_sources <- source_paths[!file.exists(source_paths)]
  expect_true(
    length(missing_sources) == 0L,
    info = paste("Missing required release sources:", paste(missing_sources, collapse = ", "))
  )
  if (length(missing_sources) != 0L) {
    return(invisible(NULL))
  }

  desc <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]
  expect_identical(unname(desc[["Version"]]), "1.0.0")

  retired_vignette_slugs <- c(
    "using-fastINFLECT-1",
    "benchmark-inflect-vs-consensus",
    "comparing-metaclustering",
    "migrating-to-fastINFLECT-2"
  )
  retired_vignette_paths <- file.path(
    root,
    "vignettes",
    paste0(retired_vignette_slugs, ".Rmd")
  )
  expect_false(
    any(file.exists(retired_vignette_paths)),
    info = paste(
      "Retired vignette sources must remain absent:",
      paste(retired_vignette_paths[file.exists(retired_vignette_paths)], collapse = ", ")
    )
  )

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
    "migrating-to-fastINFLECT-2",
    "using-fastINFLECT-1",
    "benchmark-inflect-vs-consensus",
    "comparing-metaclustering"
  )
  for (pattern in stale_patterns) {
    expect_false(grepl(pattern, release_text, fixed = TRUE), info = pattern)
  }

  pkgdown_path <- file.path(root, "_pkgdown.yml")
  if (file.exists(file.path(root, ".git"))) {
    expect_true(file.exists(pkgdown_path))
    pkgdown <- readLines(pkgdown_path, warn = FALSE)
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
    pkgdown_text <- paste(pkgdown, collapse = "\n")
    for (slug in retired_vignette_slugs) {
      expect_false(grepl(slug, pkgdown_text, fixed = TRUE), info = slug)
    }
  }
})
