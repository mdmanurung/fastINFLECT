test_that("Lfunction returns scalar knee and range without plotting", {
  source_pkg_file("leastError.R")
  source_pkg_file("angleplot.R")
  source_pkg_file("Lfunction.R")

  assign("LinesAngles", function(...) 42, envir = globalenv())
  on.exit(rm("LinesAngles", envir = globalenv()), add = TRUE)

  dataframe <- data.frame(
    x = 1:25,
    y = c(seq(0, 20, length.out = 12), seq(20, 24, length.out = 13))
  )

  full <- Lfunction(dataframe, cutoff = 1000, plot = FALSE)
  refined <- Lfunction(dataframe, cutoff = 2, plot = FALSE)

  expect_length(full$knee, 1)
  expect_length(full$range, 1)
  expect_false(is.data.frame(full$range))
  expect_length(refined$knee, 1)
  expect_length(refined$range, 1)
  expect_false(is.data.frame(refined$range))
})

test_that("Lfunction handles the smallest supported curve without a one-point second line", {
  testthat::skip_if_not_installed("LearnGeom")
  source_pkg_file("leastError.R")
  source_pkg_file("angleplot.R")
  source_pkg_file("Lfunction.R")

  assign("LinesAngles", LearnGeom::LinesAngles, envir = globalenv())
  on.exit(rm("LinesAngles", envir = globalenv()), add = TRUE)

  result <- Lfunction(
    data.frame(
      x = 5:9,
      y = c(10, 25, 45, 58, 64)
    ),
    cutoff = 1000,
    plot = FALSE
  )

  expect_length(result$knee, 1)
  expect_length(result$range, 1)
  expect_false(is.na(result$angle))
})

test_that("leastError rejects too-short input clearly", {
  source_pkg_file("leastError.R")

  expect_error(
    leastError(data.frame(x = 1:4, y = 1:4)),
    "at least five rows"
  )
})

test_that("QC.to.curve validates basedata before fitting", {
  source_pkg_file("QC-to-curve.R")

  collection <- list(
    U.set = data.frame(i = 1:6, Unimodality = c(1, 2, 3, 4, 4.5, 4.7))
  )

  expect_error(QC.to.curve(collection, basedata = "bad"), "basedata")
})
