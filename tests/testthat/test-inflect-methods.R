make_minimal_inflect_result <- function(provenance = NULL) {
  accuracy <- list(
    "5" = matrix(
      TRUE,
      nrow = 2,
      ncol = 3,
      dimnames = list(paste0("cluster", 1:2), c("CD3", "CD4", "CD8"))
    )
  )
  result <- list(
    collection.U = data.frame(
      i = c(5, 10, 15),
      Unimodality = c(50, 65, 72)
    ),
    fittedcurve = data.frame(x = 5:15, y = seq(50, 72, length.out = 11)),
    lfunction = list(knee = 10, range = 15, angle = 42),
    ggplot = structure(list(name = "stored plot"), class = "ggplot"),
    metaclustering.list = list(),
    Accuracy.sets = accuracy,
    accuracy.sets = accuracy
  )
  if (!is.null(provenance)) {
    result$provenance <- provenance
  }
  structure(result, class = c("inflect.results", "list"))
}

test_that("print.inflect.results emits a compact summary and returns invisibly", {
  source_pkg_file("inflect-results-methods.R")

  result <- make_minimal_inflect_result()
  output <- capture.output(returned <- withVisible(print(result)))

  expect_false(returned$visible)
  expect_identical(returned$value, result)
  expect_true(any(grepl("fastINFLECT result", output, fixed = TRUE)))
  expect_true(any(grepl("knee", output, fixed = TRUE)))
})

test_that("summary.inflect.results returns one row with stable base columns", {
  source_pkg_file("inflect-results-methods.R")

  result <- make_minimal_inflect_result()
  result_summary <- summary(result)

  expected_columns <- c(
    "knee",
    "range",
    "angle",
    "n_points",
    "min_i",
    "max_i",
    "min_unimodality",
    "max_unimodality",
    "n_markers"
  )

  expect_s3_class(result_summary, "data.frame")
  expect_equal(nrow(result_summary), 1)
  expect_identical(names(result_summary), expected_columns)
})

test_that("summary.inflect.results returns NA for missing lfunction fields", {
  source_pkg_file("inflect-results-methods.R")

  result <- make_minimal_inflect_result()
  result$lfunction <- list(knee = 10)
  result_summary <- summary(result)

  expect_equal(result_summary$knee, 10)
  expect_true(is.na(result_summary$range))
  expect_true(is.na(result_summary$angle))
})

test_that("plot.inflect.results returns the stored ggplot object", {
  source_pkg_file("inflect-results-methods.R")

  result <- make_minimal_inflect_result()

  expect_identical(plot(result), result$ggplot)
})

test_that("as.data.frame.inflect.results returns collection.U", {
  source_pkg_file("inflect-results-methods.R")

  result <- make_minimal_inflect_result()

  expect_identical(as.data.frame(result), result$collection.U)
})

test_that("summary/print ignore an empty-list provenance (9-column base contract holds)", {
  source_pkg_file("inflect-results-methods.R")

  result <- make_minimal_inflect_result()
  result$provenance <- list()   # as produced by new_inflect_results() defaults

  result_summary <- summary(result)
  expect_equal(ncol(result_summary), 9L)
  expect_false("uniform.test" %in% names(result_summary))

  output <- capture.output(print(result))
  expect_false(any(grepl("provenance", output, fixed = TRUE)))
})

test_that("summary.inflect.results includes provenance columns when present", {
  source_pkg_file("inflect-results-methods.R")

  result <- make_minimal_inflect_result(provenance = list(
    uniform.test = "both",
    th.pvalue = 0.05,
    th.IQR = 2,
    zeroes.in = FALSE,
    basedata = "Curve",
    package_version = "0.2.1"
  ))
  result_summary <- summary(result)

  expected_columns <- c(
    "uniform.test",
    "th.pvalue",
    "th.IQR",
    "zeroes.in",
    "basedata",
    "package_version"
  )

  expect_true(all(expected_columns %in% names(result_summary)))
})
