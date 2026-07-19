test_that("inflect_kneedle finds the knee of a concave increasing curve", {
  source_pkg_file("select-k.R")

  ## Saturating curve: steep rise then plateau; knee should sit near the bend.
  x <- 1:40
  y <- 100 * x / (x + 8)
  knee <- inflect_kneedle(x, y)
  expect_true(is.finite(knee))
  expect_true(knee >= 3 && knee <= 20)

  ## Degenerate inputs return NA rather than erroring.
  expect_true(is.na(inflect_kneedle(1:3, rep(5, 3))))
  expect_true(is.na(inflect_kneedle(numeric(0), numeric(0))))
})

test_that("inflect_kneedle orients correctly for all four curve shapes", {
  source_pkg_file("select-k.R")

  x <- 1:40
  base <- 100 * x / (x + 8)   # concave, increasing; corner sits at a low x

  k_ci <- inflect_kneedle(x, base, concave = TRUE, increasing = TRUE)
  expect_true(k_ci >= 3 && k_ci <= 20)

  ## Vertical reflection (y -> max-y) turns concave-increasing into convex-decreasing;
  ## the corner keeps its x location.
  k_vd <- inflect_kneedle(x, max(base) - base, concave = FALSE, increasing = FALSE)
  expect_equal(k_vd, k_ci)

  ## Horizontal reflection (reverse x) turns it into concave-decreasing;
  ## the corner moves to the mirrored x position.
  k_cd <- inflect_kneedle(x, rev(base), concave = TRUE, increasing = FALSE)
  expect_equal(k_cd, max(x) + min(x) - k_ci)

  ## Convex-increasing (reflect the decreasing one vertically) is finite and sensible.
  k_vi <- inflect_kneedle(x, max(base) - rev(base), concave = FALSE, increasing = TRUE)
  expect_equal(k_vi, k_cd)
})

test_that("inflect_threshold_k returns the smallest k meeting the target", {
  source_pkg_file("select-k.R")

  df <- data.frame(i = 5:12, Unimodality = c(80, 88, 93, 95, 96, 98, 99, 99.5))
  expect_equal(inflect_threshold_k(df, 0.95), 8L)   # first i with >= 95%
  expect_equal(inflect_threshold_k(df, 95), 8L)     # percentage form
  expect_equal(inflect_threshold_k(df, 0.99), 11L)
  expect_true(is.na(inflect_threshold_k(df, 0.999))) # target never reached

  expect_error(inflect_threshold_k(data.frame(a = 1), 0.95), "collection.U")
  expect_error(inflect_threshold_k(df, -1), "target")
})

test_that(".inflect_selection assembles all three criteria", {
  source_pkg_file("select-k.R")

  df <- data.frame(i = 5:20, Unimodality = 100 * (5:20) / (5:20 + 4))
  sel <- .inflect_selection(df, lfunction = list(knee = 11), target = 0.9)

  expect_s3_class(sel, "data.frame")
  expect_equal(sel$method, c("inflection", "kneedle", "threshold"))
  expect_equal(sel$k[sel$method == "inflection"], 11)
  expect_true(is.finite(sel$k[sel$method == "kneedle"]))
  ## threshold: smallest i with Unimodality >= 90
  expect_equal(sel$k[sel$method == "threshold"], inflect_threshold_k(df, 0.9))
})
