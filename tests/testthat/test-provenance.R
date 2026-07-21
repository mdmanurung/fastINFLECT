test_that("new_inflect_results stores provenance", {
  env <- new.env(parent = globalenv())
  source_pkg_file("inflect-results-class.R", envir = env)

  provenance <- list(
    set.i = as.integer(c(5, 10)),
    uniform.test = "both",
    markers = c("CD3", "CD4")
  )
  result <- env$new_inflect_results(
    collection.U = data.frame(i = 5, Unimodality = 100),
    fittedcurve = data.frame(x = 1:5, y = 1:5),
    lfunction = list(knee = 1, range = 5, angle = 0),
    ggplot = "plot",
    metaclustering.list = list("5" = as.integer(1)),
    accuracy.sets = list("5" = matrix(TRUE, nrow = 1, ncol = 1)),
    provenance = provenance
  )

  expect_identical(result$provenance, provenance)
})

test_that("new_inflect_results validates provenance", {
  env <- new.env(parent = globalenv())
  source_pkg_file("inflect-results-class.R", envir = env)

  expect_error(
    env$new_inflect_results(
      collection.U = data.frame(i = 5, Unimodality = 100),
      fittedcurve = data.frame(x = 1:5, y = 1:5),
      lfunction = list(knee = 1, range = 5, angle = 0),
      ggplot = "plot",
      provenance = "not a list"
    ),
    "`provenance` must be a list"
  )
})

test_that("normalize_set_i expands length-two inputs using current INFLECT rules", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-provenance.R", envir = env)

  fake <- make_fake_flowsom(
    matrix(rep(1, 20), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    n_nodes = 100L
  )

  expect_identical(
    env$normalize_set_i(c(10, 20), fake),
    as.integer(c(5:10, 15, 20, seq(30, 90, 10)))
  )
})

test_that("normalize_set_i keeps default k values within small kohonen maps", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-provenance.R", envir = env)

  fake <- make_fake_kohonen(
    data = matrix(rep(1:18, length.out = 36), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    codes = matrix(rep(1:18, length.out = 18), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    unit.classif = rep(1:9, length.out = 18)
  )

  expect_identical(
    env$normalize_set_i(c(150, 200), fake),
    as.integer(5:9)
  )
})

test_that("normalize_set_i rejects SOMs too small for curve fitting", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-provenance.R", envir = env)

  fake <- make_fake_kohonen(
    data = matrix(rep(1:8, length.out = 16), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    codes = matrix(rep(1:8, length.out = 8), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    unit.classif = rep(1:4, length.out = 8)
  )

  expect_error(
    env$normalize_set_i(c(150, 200), fake),
    "at least five SOM nodes"
  )
})

test_that("resolve_inflect_markers returns evaluated markers using FlowSOMQC rules", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-provenance.R", envir = env)

  fake <- make_fake_flowsom(
    matrix(
      rep(1, 40),
      ncol = 4,
      dimnames = list(NULL, c("CD10", "CD2", "CD1", "CD3"))
    ),
    cols_used = c(2, 4)
  )

  expect_equal(env$resolve_inflect_markers(fake), c("CD2", "CD3"))
  expect_equal(
    env$resolve_inflect_markers(
      fake,
      only.clustering.markers = FALSE,
      acquired_markers = c("CD10", "CD1", "CD3")
    ),
    c("CD3", "CD1", "CD10")
  )
})
