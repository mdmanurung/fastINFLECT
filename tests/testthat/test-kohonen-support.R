test_that("as_inflect_som normalizes kohonen SOM objects", {
  source_pkg_file("som-adapter.R")

  values <- matrix(
    c(1, 2, 3, 4, 2, 3, 4, 5),
    ncol = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  codes <- matrix(
    c(1, 2, 3, 4, 2, 3, 4, 5),
    ncol = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  object <- make_fake_kohonen(
    data = values,
    codes = codes,
    unit.classif = c(1L, 2L, 3L, 4L)
  )

  view <- as_inflect_som(object)

  expect_s3_class(view, "inflect_som_view")
  expect_equal(view$data, values)
  expect_equal(view$prettyColnames, c("CD3", "CD4"))
  expect_identical(view$map$colsUsed, 1:2)
  expect_equal(view$map$mapping[, 1], c(1L, 2L, 3L, 4L))
  expect_equal(view$map$nNodes, 4L)
  expect_equal(view$map$codes, codes)
  expect_equal(view$inflect_source$type, "kohonen")
})

test_that("as_inflect_som validates FlowSOM mapping and marker indices", {
  source_pkg_file("som-adapter.R")

  base <- make_fake_flowsom(
    matrix(
      rep(1:20, length.out = 20),
      ncol = 2,
      dimnames = list(NULL, c("CD3", "CD4"))
    )
  )
  base$map$nNodes <- 2L
  base$map$codes <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  base$map$mapping <- matrix(rep(1:2, length.out = nrow(base$data)), ncol = 1)

  missing_mapping <- base
  missing_mapping$map$mapping[1, 1] <- NA_integer_
  expect_error(as_inflect_som(missing_mapping), "must not contain missing values")

  out_of_range <- base
  out_of_range$map$mapping[1, 1] <- 3L
  expect_error(as_inflect_som(out_of_range), "outside `seq_len\\(map\\$nNodes\\)`")

  short_mapping <- base
  short_mapping$map$mapping <- matrix(short_mapping$map$mapping[-1, 1], ncol = 1)
  expect_error(as_inflect_som(short_mapping), "length must match")

  bad_cols <- base
  bad_cols$map$colsUsed <- c(1L, 3L)
  expect_error(as_inflect_som(bad_cols), "valid data column indices")
})

test_that("as_inflect_som uses fused numeric code layers for kohonen XYF objects", {
  source_pkg_file("som-adapter.R")

  x <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8),
    ncol = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  y <- matrix(
    c(1, 0, 0, 1),
    ncol = 1,
    dimnames = list(NULL, "response")
  )
  x.codes <- matrix(
    c(1, 2, 3, 4, 2, 3, 4, 5),
    ncol = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  y.codes <- matrix(
    c(0, 1, 0, 1),
    ncol = 1,
    dimnames = list(NULL, "response")
  )
  object <- make_fake_kohonen(
    data = list(X = x, Y = y),
    codes = list(X = x.codes, Y = y.codes),
    unit.classif = c(1L, 2L, 3L, 4L),
    distance.weights = c(1, 1),
    user.weights = c(1, 1)
  )

  view <- as_inflect_som(object)

  expect_equal(view$data, x)
  expect_equal(view$prettyColnames, c("CD3", "CD4"))
  expect_equal(dim(view$map$codes), c(4L, 3L))
  expect_equal(colnames(view$map$codes), c("X.CD3", "X.CD4", "Y.response"))
  expect_equal(view$inflect_source$code_layers, c("X", "Y"))
  expect_equal(view$inflect_source$data_layer, "X")
})

test_that("FlowSOMQC accepts kohonen SOM objects", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")

  values <- matrix(
    rep(c(1, 2), each = 6),
    ncol = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  codes <- matrix(
    c(1, 2, 1, 2),
    ncol = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  object <- make_fake_kohonen(
    data = values,
    codes = codes,
    unit.classif = rep(1:2, each = 3)
  )

  result <- FlowSOMQC(
    FlowSOM.results = object,
    metaclustering = as.integer(c(1, 2)),
    uniform.test = "spread",
    progress = FALSE
  )

  expect_true(all(result, na.rm = TRUE))
  expect_equal(colnames(result), c("CD3", "CD4"))
})

test_that("iteration.QC accepts kohonen SOM objects", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")
  source_pkg_file("iteration-QC.R")

  values <- matrix(
    rep(c(1, 2), each = 6),
    ncol = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  codes <- matrix(
    c(1, 2, 1, 2),
    ncol = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  object <- make_fake_kohonen(
    data = values,
    codes = codes,
    unit.classif = rep(1:2, each = 3)
  )
  metaclustering.list <- list("2" = as.integer(c(1, 2)))

  result <- iteration.QC(
    FlowSOM.results = object,
    metaclustering.list = metaclustering.list,
    set.i = 2,
    workers = 1L,
    uniform.test = "spread",
    progress = FALSE
  )

  expect_equal(result$U.set$i, 2)
  expect_equal(names(result$Accuracy.matrixes), "2")
  expect_equal(colnames(result$Accuracy.matrixes[[1]]), c("CD3", "CD4"))
})

test_that("iteration.metacluster accepts kohonen XYF objects", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("iteration-metacluster.R")

  x.codes <- matrix(
    c(1, 2, 8, 9, 1, 2, 8, 9),
    ncol = 2,
    dimnames = list(NULL, c("CD3", "CD4"))
  )
  y.codes <- matrix(
    c(0, 0, 1, 1),
    ncol = 1,
    dimnames = list(NULL, "response")
  )
  object <- make_fake_kohonen(
    data = list(
      X = matrix(rep(1:8, each = 1), ncol = 2, dimnames = list(NULL, c("CD3", "CD4"))),
      Y = matrix(c(0, 0, 1, 1), ncol = 1, dimnames = list(NULL, "response"))
    ),
    codes = list(X = x.codes, Y = y.codes),
    unit.classif = c(1L, 2L, 3L, 4L)
  )

  result <- iteration.metacluster(object, set.i = c(2, 3))

  expect_equal(names(result), c("2", "3"))
  expect_true(all(vapply(result, is.integer, logical(1))))
  expect_equal(unname(lengths(result)), c(4L, 4L))
})

test_that("as_inflect_som accepts real kohonen SOM and XYF objects", {
  testthat::skip_if_not_installed("kohonen")
  source_pkg_file("som-adapter.R")

  set.seed(1)
  rows <- c(1:10, 51:60, 101:110)
  x <- scale(as.matrix(iris[rows, 1:4]))
  y <- stats::model.matrix(~ Species - 1, iris[rows, ])

  som <- kohonen::som(
    x,
    grid = kohonen::somgrid(3, 3, "hexagonal"),
    rlen = 5
  )
  xyf <- kohonen::xyf(
    x,
    y,
    grid = kohonen::somgrid(3, 3, "hexagonal"),
    rlen = 5
  )

  som.view <- as_inflect_som(som)
  xyf.view <- as_inflect_som(xyf)

  expect_equal(nrow(som.view$map$codes), 9L)
  expect_equal(ncol(som.view$data), 4L)
  expect_equal(nrow(xyf.view$map$codes), 9L)
  expect_equal(
    ncol(xyf.view$map$codes),
    ncol(xyf$codes[[1]]) + ncol(xyf$codes[[2]])
  )
  expect_equal(xyf.view$prettyColnames, colnames(x))
  expect_equal(xyf.view$inflect_source$type, "kohonen")
  expect_length(xyf.view$inflect_source$code_layers, 2L)
})

test_that("INFLECT accepts real kohonen SOM and XYF objects with literal set.i", {
  testthat::skip_if_not_installed("kohonen")
  testthat::skip_if_not_installed("LearnGeom")
  source_pkg_file("som-adapter.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")
  source_pkg_file("iteration-QC.R")
  source_pkg_file("leastError.R")
  source_pkg_file("angleplot.R")
  source_pkg_file("Lfunction.R")
  source_pkg_file("QC-to-curve.R")
  source_pkg_file("inflect-provenance.R")
  source_pkg_file("inflect-results-class.R")
  source_pkg_file("select-k.R")
  source_pkg_file("INFLECT.R")

  assign("LinesAngles", LearnGeom::LinesAngles, envir = globalenv())
  on.exit(rm("LinesAngles", envir = globalenv()), add = TRUE)

  set.seed(1)
  rows <- c(1:10, 51:60, 101:110)
  x <- scale(as.matrix(iris[rows, 1:4]))
  y <- stats::model.matrix(~ Species - 1, iris[rows, ])

  som <- kohonen::som(
    x,
    grid = kohonen::somgrid(3, 3, "hexagonal"),
    rlen = 5
  )
  xyf <- kohonen::xyf(
    x,
    y,
    grid = kohonen::somgrid(3, 3, "hexagonal"),
    rlen = 5
  )

  som.result <- INFLECT(
    som,
    set.i = 5:9,
    workers = 1L,
    uniform.test = "spread",
    progress = FALSE
  )
  xyf.result <- INFLECT(
    xyf,
    set.i = 5:9,
    workers = 1L,
    uniform.test = "spread",
    progress = FALSE
  )

  expect_s3_class(som.result, "inflect.results")
  expect_s3_class(xyf.result, "inflect.results")
  expect_equal(som.result$provenance$set.i, 5:9)
  expect_equal(xyf.result$provenance$set.i, 5:9)
  expect_equal(som.result$provenance$som_type, "kohonen")
  expect_equal(xyf.result$provenance$som_type, "kohonen")
})
