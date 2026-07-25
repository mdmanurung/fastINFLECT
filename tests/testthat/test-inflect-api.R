test_that("INFLECT forwards QC arguments and returns an inflect.results object", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-qc-core.R", envir = env)
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
                               max.events.per.node = NULL,
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
    captured$max.events.per.node <- max.events.per.node
    captured$seed <- seed
    captured$dots <- list(...)

    matrices <- stats::setNames(
      rep(
        list(matrix(TRUE, nrow = 1, ncol = 1, dimnames = list("1", "CD3"))),
        length(set.i)
      ),
      as.character(set.i)
    )
    details <- lapply(matrices, function(x) {
      list(
        dip_pass = x,
        iqr_pass = x,
        combined_pass = x,
        criterion_pass = x
      )
    })
    list(
      scores = data.frame(
        k = set.i,
        qc_pass_rate = 100,
        criterion = "iqr"
      ),
      criterion_pass = matrices,
      dip_pass = matrices,
      iqr_pass = matrices,
      combined_pass = matrices,
      qc.details = details,
      provenance = list(
        criterion = "iqr",
        criterion_label = "IQR-spread QC pass rate",
        scoring_mode = "indexed_node_capped",
        event_sampling = list(
          mode = "capped_per_node",
          max_events_per_node = max.events.per.node,
          original_events = 10L,
          retained_events = 10L,
          per_node = data.frame(
            node = seq_len(40L),
            original_events = c(10L, rep(0L, 39L)),
            retained_events = c(10L, rep(0L, 39L))
          )
        )
      )
    )
  }
  env$QC.to.curve <- function(collection.U, basedata, ggtitle) {
    list(
      scores = collection.U$scores,
      collection.U = data.frame(
        i = collection.U$scores$k,
        Unimodality = collection.U$scores$qc_pass_rate
      ),
      fittedcurve = data.frame(
        k = 1:5,
        qc_pass_rate = 1:5,
        x = 1:5,
        y = 1:5
      ),
      lfunction = list(knee = 1, range = 5, angle = 0),
      ggplot = "plot"
    )
  }

  fake <- make_fake_flowsom(
    matrix(rep(1:10, 2), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
    cols_used = 1,
    n_nodes = 40L
  )
  expected_set_i <- 7:11
  expect_identical(eval(formals(env$INFLECT)$zeroes.in), TRUE)
  expect_error(env$INFLECT(fake), "`set.i` is required")
  expect_error(
    env$INFLECT(fake, set.i = 5:9, multicore = TRUE, cores = 1L),
    "`cores`.*>= 2"
  )
  result <- suppressWarnings(env$INFLECT(
    fake,
    set.i = expected_set_i,
    multicore = FALSE,
    zeroes.in = TRUE,
    only.clustering.markers = FALSE,
    acquired_markers = "CD4",
    uniform.test = "spread",
    th.pvalue = 0.2,
    th.IQR = 4,
    max.events.per.node = 1000L,
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
  expect_equal(captured$max.events.per.node, 1000L)
  expect_equal(captured$seed, 1L)
  expect_true(captured$dots$simulate.p.value)
  expect_s3_class(result, "inflect.results")
  expect_named(result$scores, c("k", "qc_pass_rate", "criterion"))
  expect_identical(result$criterion_pass, result$iqr_pass)
  expect_length(result$qc.details, length(expected_set_i))
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
  expect_equal(result$provenance$criterion, "iqr")
  expect_equal(result$provenance$max.events.per.node, 1000L)
  expect_equal(result$provenance$scoring_mode, "indexed_node_capped")
  expect_equal(result$provenance$n_events_retained, 10L)
  expect_true(result$provenance$diptest_args$simulate.p.value)
  inflection <- result$selection[result$selection$method == "inflection", ]
  expect_false(inflection$partition_available)
  expect_equal(inflection$k_status, "fitted_estimate_no_partition")
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
  expect_equal(result$scores$qc_pass_rate, 100)
})

test_that("INFLECT normalizes kohonen inputs before running pipeline stages", {
  env <- new.env(parent = globalenv())
  source_pkg_file("som-adapter.R", envir = env)
  source_pkg_file("inflect-qc-core.R", envir = env)
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
    matrices <- stats::setNames(
      rep(
        list(matrix(TRUE, nrow = 1, ncol = 2, dimnames = list("1", c("CD3", "CD4")))),
        length(set.i)
      ),
      as.character(set.i)
    )
    list(
      scores = data.frame(
        k = set.i,
        qc_pass_rate = seq_along(set.i),
        criterion = "iqr"
      ),
      criterion_pass = matrices,
      dip_pass = matrices,
      iqr_pass = matrices,
      combined_pass = matrices,
      qc.details = rep(list(list()), length(set.i)),
      provenance = list()
    )
  }
  env$QC.to.curve <- function(collection.U, basedata, ggtitle) {
    list(
      scores = collection.U$scores,
      collection.U = data.frame(
        i = collection.U$scores$k,
        Unimodality = collection.U$scores$qc_pass_rate
      ),
      fittedcurve = data.frame(
        k = collection.U$scores$k,
        qc_pass_rate = collection.U$scores$qc_pass_rate,
        x = collection.U$scores$k,
        y = collection.U$scores$qc_pass_rate
      ),
      lfunction = list(
        knee = collection.U$scores$k[[1]],
        range = collection.U$scores$k[[5]],
        angle = 0
      ),
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
