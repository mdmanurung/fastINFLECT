test_that("INFLECT forwards QC arguments and returns an inflect.results object", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-provenance.R", envir = env)
  source_pkg_file("inflect-results-class.R", envir = env)
  source_pkg_file("select-k.R", envir = env)
  source_pkg_file("INFLECT.R", envir = env)

  captured <- new.env(parent = emptyenv())
  env$iteration.metacluster <- function(FlowSOM.results, set.i, multicore, cores) {
    stats::setNames(
      rep(list(as.integer(1)), length(set.i)),
      as.character(set.i)
    )
  }
  env$iteration.QC <- function(FlowSOM.results,
                               metaclustering.list,
                               set.i,
                               multicore,
                               cores,
                               zeroes.in,
                               only.clustering.markers,
                               acquired_markers,
                               uniform.test,
                               th.pvalue,
                               th.IQR,
                               verbose,
                               max.n.diptest = NULL,
                               seed = 1L,
                               ...) {
    captured$set.i <- set.i
    captured$zeroes.in <- zeroes.in
    captured$only.clustering.markers <- only.clustering.markers
    captured$acquired_markers <- acquired_markers
    captured$uniform.test <- uniform.test
    captured$th.pvalue <- th.pvalue
    captured$th.IQR <- th.IQR
    captured$verbose <- verbose
    captured$max.n.diptest <- max.n.diptest
    captured$seed <- seed
    captured$dots <- list(...)

    list(
      U.set = data.frame(i = set.i, Unimodality = 100),
      Accuracy.matrixes = stats::setNames(
        rep(
          list(matrix(TRUE, nrow = 1, ncol = 1, dimnames = list("1", "CD3"))),
          length(set.i)
        ),
        as.character(set.i)
      )
    )
  }
  env$QC.to.curve <- function(collection.U, basedata, ggtitle) {
    list(
      collection.U = collection.U[[1]],
      fittedcurve = data.frame(x = 1:5, y = 1:5),
      lfunction = list(knee = 1, range = 5, angle = 0),
      ggplot = "plot"
    )
  }

  fake <- make_fake_flowsom(
    matrix(rep(1:10, 2), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    cols_used = 1
  )
  fake$map$nNodes <- 40L
  expected_set_i <- as.integer(c(5:7, 12, 17, 27))
  result <- suppressWarnings(env$INFLECT(
    fake,
    set.i = c(7, 17),
    multicore = FALSE,
    zeroes.in = TRUE,
    only.clustering.markers = FALSE,
    acquired_markers = "CD4",
    uniform.test = "spread",
    th.pvalue = 0.2,
    th.IQR = 4,
    verbose = TRUE,
    simulate.p.value = TRUE
  ))

  expect_identical(captured$set.i, expected_set_i)
  expect_true(captured$zeroes.in)
  expect_false(captured$only.clustering.markers)
  expect_equal(captured$acquired_markers, "CD4")
  expect_equal(captured$uniform.test, "spread")
  expect_equal(captured$th.pvalue, 0.2)
  expect_equal(captured$th.IQR, 4)
  expect_true(captured$verbose)
  expect_null(captured$max.n.diptest)
  expect_equal(captured$seed, 1L)
  expect_true(captured$dots$simulate.p.value)
  expect_s3_class(result, "inflect.results")
  expect_identical(result$Accuracy.sets, result$accuracy.sets)
  expect_identical(result$provenance$set.i, expected_set_i)
  expect_equal(result$provenance$zeroes.in, TRUE)
  expect_equal(result$provenance$basedata, "Curve")
  expect_equal(result$provenance$only.clustering.markers, FALSE)
  expect_equal(result$provenance$acquired_markers, "CD4")
  expect_equal(result$provenance$markers, "CD4")
  expect_equal(result$provenance$uniform.test, "spread")
  expect_equal(result$provenance$th.pvalue, 0.2)
  expect_equal(result$provenance$th.IQR, 4)
  expect_true(result$provenance$diptest_args$simulate.p.value)
})

test_that("new_inflect_results creates S3 result objects with both accuracy aliases", {
  env <- new.env(parent = globalenv())
  source_pkg_file("inflect-results-class.R", envir = env)

  accuracy <- list("5" = matrix(TRUE, nrow = 1, ncol = 1))
  result <- env$new_inflect_results(
    collection.U = data.frame(i = 5, Unimodality = 100),
    fittedcurve = data.frame(x = 1:5, y = 1:5),
    lfunction = list(knee = 1, range = 5, angle = 0),
    ggplot = "plot",
    metaclustering.list = list("5" = as.integer(1)),
    accuracy.sets = accuracy
  )

  expect_s3_class(result, "inflect.results")
  expect_identical(result$Accuracy.sets, accuracy)
  expect_identical(result$accuracy.sets, accuracy)
})

test_that("INFLECT normalizes kohonen inputs before running pipeline stages", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-provenance.R", envir = env)
  source_pkg_file("inflect-results-class.R", envir = env)
  source_pkg_file("select-k.R", envir = env)
  source_pkg_file("INFLECT.R", envir = env)

  captured <- new.env(parent = emptyenv())
  env$iteration.metacluster <- function(FlowSOM.results, set.i, multicore, cores) {
    captured$metacluster_view <- inherits(FlowSOM.results, "inflect_som_view")
    captured$source_type <- FlowSOM.results$inflect_source$type
    stats::setNames(
      rep(list(as.integer(c(1, 2, 3, 4))), length(set.i)),
      as.character(set.i)
    )
  }
  env$iteration.QC <- function(FlowSOM.results,
                               metaclustering.list,
                               set.i,
                               multicore,
                               cores,
                               zeroes.in,
                               only.clustering.markers,
                               acquired_markers,
                               uniform.test,
                               th.pvalue,
                               th.IQR,
                               verbose,
                               ...) {
    captured$qc_view <- inherits(FlowSOM.results, "inflect_som_view")
    list(
      U.set = data.frame(i = set.i, Unimodality = seq_along(set.i)),
      Accuracy.matrixes = stats::setNames(
        rep(
          list(matrix(TRUE, nrow = 1, ncol = 2, dimnames = list("1", c("CD3", "CD4")))),
          length(set.i)
        ),
        as.character(set.i)
      )
    )
  }
  env$QC.to.curve <- function(collection.U, basedata, ggtitle) {
    list(
      collection.U = collection.U[[1]],
      fittedcurve = data.frame(x = collection.U[[1]]$i, y = collection.U[[1]]$Unimodality),
      lfunction = list(knee = collection.U[[1]]$i[[1]], range = collection.U[[1]]$i[[5]], angle = 0),
      ggplot = "plot"
    )
  }

  object <- make_fake_kohonen(
    data = matrix(rep(1:10, 2), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    codes = matrix(rep(1:8, 2), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    unit.classif = rep(1:4, length.out = 10)
  )

  result <- suppressWarnings(env$INFLECT(
    object,
    set.i = 2:6,
    multicore = FALSE,
    uniform.test = "spread"
  ))

  expect_true(captured$metacluster_view)
  expect_true(captured$qc_view)
  expect_equal(captured$source_type, "kohonen")
  expect_equal(result$provenance$som_type, "kohonen")
  expect_equal(result$provenance$som_data_layer, "X")
  expect_equal(result$provenance$markers, c("CD3", "CD4"))
})
