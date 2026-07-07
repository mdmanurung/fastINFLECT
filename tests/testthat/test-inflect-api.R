test_that("INFLECT forwards QC arguments and returns an inflect.results object", {
  env <- new.env(parent = globalenv())
  source_pkg_file("inflect-results-class.R", envir = env)
  source_pkg_file("INFLECT.R", envir = env)

  captured <- new.env(parent = emptyenv())
  env$iteration.metacluster <- function(FlowSOM.results, set.i, multicore, cores) {
    stats::setNames(list(as.integer(1)), as.character(set.i))
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
    captured$uniform.test <- uniform.test
    captured$th.pvalue <- th.pvalue
    captured$th.IQR <- th.IQR
    captured$verbose <- verbose
    captured$dots <- list(...)

    list(
      U.set = data.frame(i = set.i, Unimodality = 100),
      Accuracy.matrixes = stats::setNames(
        list(matrix(TRUE, nrow = 1, ncol = 1, dimnames = list("1", "CD3"))),
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

  fake <- make_fake_flowsom(matrix(1:10, ncol = 1, dimnames = list(NULL, "CD3")))
  result <- suppressWarnings(env$INFLECT(
    fake,
    set.i = 5,
    multicore = FALSE,
    uniform.test = "spread",
    th.pvalue = 0.2,
    th.IQR = 4,
    verbose = TRUE,
    simulate.p.value = TRUE
  ))

  expect_equal(captured$uniform.test, "spread")
  expect_equal(captured$th.pvalue, 0.2)
  expect_equal(captured$th.IQR, 4)
  expect_true(captured$verbose)
  expect_true(captured$dots$simulate.p.value)
  expect_s3_class(result, "inflect.results")
  expect_identical(result$Accuracy.sets, result$accuracy.sets)
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
