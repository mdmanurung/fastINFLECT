test_that("bimodality.coefficient separates unimodal from bimodal samples", {
  source_pkg_file("bimodality.R")

  set.seed(1)
  unimodal <- bimodality.coefficient(rnorm(2000))
  bimodal <- bimodality.coefficient(c(rnorm(1000), rnorm(1000, 8)))

  expect_true(is.finite(unimodal) && is.finite(bimodal))
  expect_lt(unimodal, 5 / 9)   # below the uniform benchmark
  expect_gt(bimodal, 5 / 9)    # clearly bimodal
  expect_gt(bimodal, unimodal)
})

test_that("bimodality.coefficient handles degenerate input", {
  source_pkg_file("bimodality.R")

  expect_true(is.na(bimodality.coefficient(c(1, 2, 3))))   # n < 4
  expect_true(is.na(bimodality.coefficient(rep(2, 10))))   # zero variance
  expect_true(is.na(bimodality.coefficient(c(1, NA, 2, NA), na.rm = TRUE)))
})

test_that("bimodality.coefficient matches its closed-form definition", {
  source_pkg_file("bimodality.R")

  set.seed(2)
  x <- rgamma(500, shape = 2)
  n <- length(x)
  m <- mean(x); d <- x - m
  m2 <- mean(d^2); m3 <- mean(d^3); m4 <- mean(d^4)
  G1 <- (m3 / m2^1.5) * sqrt(n * (n - 1)) / (n - 2)
  G2 <- ((n + 1) * (m4 / m2^2 - 3) + 6) * (n - 1) / ((n - 2) * (n - 3))
  expected <- (G1^2 + 1) / (G2 + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))

  expect_equal(bimodality.coefficient(x), expected)
})
