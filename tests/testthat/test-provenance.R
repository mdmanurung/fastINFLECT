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

test_that("normalize_set_i treats valid schedules literally", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-qc-core.R", envir = env)
  source_pkg_file("inflect-provenance.R", envir = env)

  fake <- make_fake_flowsom(
    matrix(rep(1, 20), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    n_nodes = 100L
  )

  schedule <- seq.int(25L, 100L, by = 5L)
  expect_identical(env$normalize_set_i(schedule, fake), schedule)
})

test_that("normalize_set_i rejects ambiguous or invalid schedules early", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-qc-core.R", envir = env)
  source_pkg_file("inflect-provenance.R", envir = env)

  fake <- make_fake_flowsom(
    matrix(rep(1, 20), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    n_nodes = 100L
  )

  expect_error(env$normalize_set_i(c(25, 100), fake), "at least five")
  expect_error(
    env$normalize_set_i(c(5, 10, 10, 15, 20), fake),
    "unique"
  )
  expect_error(
    env$normalize_set_i(c(5, 10, 15.5, 20, 25), fake),
    "integer"
  )
  expect_error(
    env$normalize_set_i(c(5, 15, 10, 20, 25), fake),
    "strictly increasing"
  )
  expect_error(
    env$normalize_set_i(c(5, 10, 15, 20, 101), fake),
    "between 1 and the number of SOM nodes"
  )
})

test_that("inflect_adaptive_set_i is explicit and bounded", {
  env <- new.env(parent = globalenv())
  source_pkg_file("inflect-provenance.R", envir = env)

  expected <- c(seq.int(5L, 25L), seq.int(30L, 100L, by = 5L))
  expect_identical(
    env$inflect_adaptive_set_i(n_nodes = 900L, max_k = 100L),
    expected
  )
  expect_equal(max(env$inflect_adaptive_set_i(900L, 103L)), 103L)
  expect_error(env$inflect_adaptive_set_i(100L, 101L), "`max_k`")
  expect_no_warning(
    overflow_safe <- env$inflect_adaptive_set_i(
      n_nodes = .Machine$integer.max,
      max_k = 100L,
      dense_until = .Machine$integer.max,
      medium_until = .Machine$integer.max
    )
  )
  expect_identical(overflow_safe, 5:100)
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
