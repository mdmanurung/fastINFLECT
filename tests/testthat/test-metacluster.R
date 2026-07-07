test_that("iteration.metacluster matches repeated hclust metaclustering calls", {
  source_pkg_file("metaClustering-hclust.R")
  source_pkg_file("iteration-metacluster.R")

  codes <- matrix(
    c(
      0, 0,
      0, 1,
      4, 4,
      4, 5,
      8, 8
    ),
    ncol = 2,
    byrow = TRUE
  )
  flowsom <- list(map = list(codes = codes))
  set.i <- 2:4

  expected <- stats::setNames(
    lapply(set.i, function(k) metaClusteringhclust(codes, nClus = k)),
    as.character(set.i)
  )
  observed <- iteration.metacluster(flowsom, set.i = set.i, multicore = FALSE)

  expect_named(observed, as.character(set.i))
  expect_equal(observed, expected)
})
