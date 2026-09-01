## Build a small but non-trivial FlowSOM: several SOM nodes, several events each,
## drawn so that some (cluster, marker) combinations are bimodal and some are not.
make_multinode_flowsom <- function(seed = 1L) {
  set.seed(seed)
  n_nodes <- 12L
  per_node <- 40L
  markers <- c("CD3", "CD4", "CD8", "CD19")
  centres <- matrix(runif(n_nodes * length(markers), 0, 6), nrow = n_nodes,
                    dimnames = list(NULL, markers))
  data <- do.call(rbind, lapply(seq_len(n_nodes), function(i) {
    m <- matrix(rnorm(per_node * length(markers), 0, 0.4), ncol = length(markers))
    sweep(m, 2, centres[i, ], `+`)
  }))
  colnames(data) <- markers
  mapping <- matrix(rep(seq_len(n_nodes), each = per_node), ncol = 1)
  obj <- list(
    data = data,
    scale = FALSE,
    prettyColnames = markers,
    map = list(
      colsUsed = seq_along(markers),
      mapping = mapping,
      nNodes = n_nodes,
      codes = centres
    )
  )
  class(obj) <- "FlowSOM"
  obj
}

legacy_accuracy_row <- function(expr,
                                uniform.test,
                                th.pvalue,
                                th.IQR,
                                p_of) {
  markers <- colnames(expr)
  out <- stats::setNames(rep(NA, length(markers)), markers)
  if (nrow(expr) <= 1L) {
    return(out)
  }
  for (j in seq_along(markers)) {
    me <- expr[, j]
    uniform <- TRUE
    if (uniform.test != "spread") {
      uniform <- p_of(me) >= th.pvalue
    }
    if (uniform && uniform.test != "unimodality") {
      uniform <- .inflect_iqr(me) < th.IQR
    }
    out[[j]] <- uniform
  }
  out
}

test_that("indexed scoring exactly matches the former full-matrix row path", {
  source_pkg_file("inflect-qc-core.R")

  fs <- make_multinode_flowsom(11L)
  p_of <- .inflect_make_p_of()
  rows_to_test <- list(1L, 1:40, c(1:17, 101:160), seq_len(nrow(fs$data)))
  for (rows in rows_to_test) {
    expr <- fs$data[rows, , drop = FALSE]
    for (uniform.test in c("both", "spread", "unimodality")) {
      expected <- legacy_accuracy_row(
        expr = expr,
        uniform.test = uniform.test,
        th.pvalue = 0.05,
        th.IQR = 2,
        p_of = p_of
      )
      actual <- .inflect_accuracy_row_indexed(
        data = fs$data,
        rows = rows,
        uniform.test = uniform.test,
        th.pvalue = 0.05,
        th.IQR = 2,
        p_of = p_of
      )
      expect_identical(
        actual,
        expected,
        info = paste(uniform.test, length(rows))
      )
    }
  }
})

test_that("memoised iteration.QC reproduces per-k FlowSOMQC exactly", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")
  source_pkg_file("metaClustering-hclust.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom()
  set.i <- 3:10
  ml <- iteration.metacluster(fs, set.i = set.i)

  reference <- lapply(set.i, function(k) {
    FlowSOMQC(fs, ml[[as.character(k)]], progress = FALSE)
  })
  ref_scores <- vapply(reference, function(m) sum(m, na.rm = TRUE) * 100 / prod(dim(m)), numeric(1))

  qc <- iteration.QC(fs, ml, set.i, workers = 1L, progress = FALSE)

  expect_equal(qc$U.set$Unimodality, ref_scores)
  expect_equal(qc$scores$qc_pass_rate, ref_scores)
  for (idx in seq_along(set.i)) {
    expect_identical(
      unname(qc$Accuracy.matrixes[[idx]]),
      unname(qc$qc.details[[idx]]$criterion_pass)
    )
    expect_equal(
      unname(qc$Accuracy.matrixes[[idx]]),
      unname(reference[[idx]]),
      ignore_attr = TRUE
    )
  }
})

test_that("memoisation holds across uniform.test modes", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")
  source_pkg_file("metaClustering-hclust.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(2L)
  set.i <- 3:8
  ml <- iteration.metacluster(fs, set.i = set.i)

  for (ut in c("both", "spread", "unimodality")) {
    ref <- lapply(set.i, function(k) {
      FlowSOMQC(fs, ml[[as.character(k)]], uniform.test = ut, progress = FALSE)
    })
    qc <- iteration.QC(
      fs, ml, set.i, workers = 1L,
      uniform.test = ut, progress = FALSE
    )
    ok <- all(mapply(function(a, b) isTRUE(all.equal(
      unname(a),
      unname(b),
      check.attributes = FALSE
    )),
                     qc$Accuracy.matrixes, ref))
    expect_true(ok, info = ut)
  }
})

test_that("serial and two-core QC results are identical", {
  skip_if(.Platform$OS.type != "unix")
  skip_if(parallel::detectCores() < 2L)
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("metaClustering-hclust.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(13L)
  set.i <- 3:8
  ml <- iteration.metacluster(fs, set.i = set.i)
  serial <- suppressWarnings(iteration.QC(
    fs, ml, set.i,
    workers = 1L,
    uniform.test = "both"
  ))
  parallel <- suppressWarnings(iteration.QC(
    fs, ml, set.i,
    workers = 2L,
    uniform.test = "both"
  ))

  expect_identical(parallel$U.set, serial$U.set)
  expect_identical(parallel$Accuracy.matrixes, serial$Accuracy.matrixes)
})

test_that("invalid worker results report task and subtree identifiers", {
  source_pkg_file("inflect-qc-core.R")

  markers <- c("CD3", "CD4")
  good <- stats::setNames(c(TRUE, FALSE), markers)
  expect_error(
    .inflect_validate_worker_results(
      list(good, NULL),
      keys = c("1,2", "3,4"),
      marker_names = markers,
      parallel = TRUE
    ),
    "Parallel QC worker failure.*task 2 \\(subtree 3,4\\).*NULL"
  )

  failed <- structure("injected worker error", class = "try-error")
  expect_error(
    .inflect_validate_worker_results(
      list(failed),
      keys = "5,6",
      marker_names = markers,
      parallel = TRUE
    ),
    "task 1 \\(subtree 5,6\\).*try-error.*injected worker error"
  )
  expect_error(
    .inflect_validate_worker_results(
      list(good),
      keys = c("7,8", "9,10"),
      marker_names = markers,
      parallel = TRUE
    ),
    "task 2 \\(subtree 9,10\\): missing result"
  )

  details <- .inflect_empty_qc_row(markers)
  details$dip_p_value <- NULL
  expect_error(
    .inflect_validate_worker_results(
      list(details),
      keys = "11,12",
      marker_names = markers,
      parallel = FALSE
    ),
    "missing field\\(s\\): dip_p_value"
  )
})

test_that("fast dip p-value matches diptest::dip.test p-values", {
  skip_if_not_installed("diptest")
  source_pkg_file("inflect-qc-core.R")

  set.seed(42)
  max_diff <- 0
  for (n in c(4:12, 20, 60, 250)) {
    for (gen in list(function(m) rnorm(m),
                     function(m) c(rnorm(m %/% 2), rnorm(m - m %/% 2, 6)),
                     function(m) rexp(m))) {
      for (rep in 1:6) {
        x <- gen(n)
        p_ref <- suppressWarnings(diptest::dip.test(x)$p.value)
        p_fast <- suppressWarnings(.inflect_dip_pvalue(diptest::dip(x), length(x)))
        max_diff <- max(max_diff, abs(p_ref - p_fast))
      }
    }
  }
  expect_lt(max_diff, 1e-12)

  x <- rnorm(20)
  D <- rep(suppressWarnings(diptest::dip(x)), 3L)
  expect_equal(.inflect_dip_pvalue(D, length(x)),
               rep(.inflect_dip_pvalue(D[[1L]], length(x)), 3L))
  expect_error(.inflect_dip_pvalue(D, c(length(x), length(x) + 1L)),
               "same length")
})

test_that("small-sample dip interpolation is warning-free and exact", {
  skip_if_not_installed("diptest")
  source_pkg_file("inflect-qc-core.R")

  set.seed(142)
  for (n in 4:8) {
    x <- rnorm(n)
    expect_no_warning(
      p_fast <- .inflect_dip_pvalue(diptest::dip(x), n)
    )
    expect_equal(
      p_fast,
      suppressWarnings(diptest::dip.test(x)$p.value),
      tolerance = 1e-12
    )
  }
})

test_that("non-finite marker values fail before QC scoring", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")
  source_pkg_file("metaClustering-hclust.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(4L)
  fs$data[c(3, 40), "CD3"] <- NA_real_
  fs$data[81, "CD3"] <- Inf
  mc <- as.integer(rep(1:3, length.out = fs$map$nNodes))

  expect_error(
    FlowSOMQC(fs, mc, uniform.test = "unimodality", progress = FALSE),
    "only finite values.*CD3 \\(3\\)"
  )
})

test_that("spread scoring uses full marker expression when dip subsampling is enabled", {
  source_pkg_file("inflect-qc-core.R")

  expr <- matrix(
    c(rep(0, 50), rep(100, 50)),
    ncol = 1,
    dimnames = list(NULL, "CD3")
  )
  full <- .inflect_accuracy_row(
    expr = expr,
    uniform.test = "spread",
    th.pvalue = 0.05,
    th.IQR = 50,
    p_of = function(x) 1,
    subsample = NULL,
    seed = 1L
  )
  capped <- .inflect_accuracy_row(
    expr = expr,
    uniform.test = "spread",
    th.pvalue = 0.05,
    th.IQR = 50,
    p_of = function(x) 1,
    subsample = 4L,
    seed = 1L
  )

  expect_false(unname(full[["CD3"]]))
  expect_identical(capped, full)
})

test_that("metaclustering inputs fail early instead of dropping events", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(5L)
  n_nodes <- fs$map$nNodes

  expect_error(
    FlowSOMQC(fs, as.integer(rep(1L, n_nodes - 1L)), progress = FALSE),
    "one entry per SOM node"
  )
  expect_error(
    FlowSOMQC(fs, as.integer(rep(1L, n_nodes + 1L)), progress = FALSE),
    "one entry per SOM node"
  )
  expect_error(
    FlowSOMQC(fs, as.integer(c(rep(1L, n_nodes - 1L), NA)), progress = FALSE),
    "must not contain missing values"
  )
  expect_error(
    FlowSOMQC(fs, as.integer(c(rep(1L, n_nodes - 1L), 0L)), progress = FALSE),
    "positive integers"
  )

  gapped <- FlowSOMQC(
    fs,
    as.integer(rep(c(2L, 4L), length.out = n_nodes)),
    uniform.test = "spread",
    progress = FALSE
  )
  expect_equal(nrow(gapped), 2L)
  expect_equal(rownames(gapped), c("1", "2"))

  ml <- list("2" = as.integer(rep(1L, n_nodes - 1L)))
  expect_error(
    iteration.QC(fs, ml, set.i = 2L, workers = 1L, progress = FALSE),
    "one entry per SOM node"
  )
})

test_that("max.n.diptest is validated once at the QC boundary", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(6L)
  ml <- list("2" = as.integer(rep(1:2, length.out = fs$map$nNodes)))
  bad_caps <- list(0L, 3L, NA_integer_, -1L, c(4L, 5L), 4.5)

  for (cap in bad_caps) {
    expect_error(
      iteration.QC(fs, ml, set.i = 2L, workers = 1L,
                   max.n.diptest = cap, progress = FALSE),
      "`max.n.diptest`"
    )
  }
  expect_identical(.inflect_validate_max_n_diptest(4L), 4L)
  expect_no_error(
    suppressWarnings(iteration.QC(fs, ml, set.i = 2L, workers = 1L,
                 max.n.diptest = 10L, progress = FALSE))
  )
})

test_that("size-robust subsampling is deterministic and bounded", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")
  source_pkg_file("metaClustering-hclust.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(3L)
  set.i <- 3:8
  ml <- iteration.metacluster(fs, set.i = set.i)

  a <- suppressWarnings(iteration.QC(
    fs, ml, set.i, workers = 1L,
    max.n.diptest = 30L, seed = 7L, progress = FALSE
  ))
  b <- suppressWarnings(iteration.QC(
    fs, ml, set.i, workers = 1L,
    max.n.diptest = 30L, seed = 7L, progress = FALSE
  ))
  expect_equal(a$U.set, b$U.set)

  ## the RNG state must be restored: an outer sequence is unaffected
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  suppressWarnings(invisible(iteration.QC(
    fs, ml, set.i, workers = 1L,
    max.n.diptest = 30L, seed = 1L, progress = FALSE
  )))
  after <- runif(1)
  expect_equal(before, after)
})

test_that("node-capped sampling is deterministic, bounded, and worker invariant", {
  skip_if(.Platform$OS.type != "unix")
  skip_if(parallel::detectCores() < 2L)
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("metaClustering-hclust.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(17L)
  set.i <- 3:8
  ml <- iteration.metacluster(fs, set.i = set.i)
  expect_warning(serial_a <- iteration.QC(
    fs, ml, set.i,
    workers = 1L,
    uniform.test = "both",
    max.events.per.node = 15L,
    seed = 42L,
    progress = FALSE
  ), "sensitivity analysis")
  serial_b <- suppressWarnings(iteration.QC(
    fs, ml, set.i,
    workers = 1L,
    uniform.test = "both",
    max.events.per.node = 15L,
    seed = 42L,
    progress = FALSE
  ))
  two_core <- suppressWarnings(iteration.QC(
    fs, ml, set.i,
    workers = 2L,
    uniform.test = "both",
    max.events.per.node = 15L,
    seed = 42L,
    progress = FALSE
  ))

  expect_identical(serial_a$U.set, serial_b$U.set)
  expect_identical(serial_a$Accuracy.matrixes, serial_b$Accuracy.matrixes)
  expect_identical(serial_a$U.set, two_core$U.set)
  expect_identical(serial_a$Accuracy.matrixes, two_core$Accuracy.matrixes)

  sampling <- serial_a$provenance$event_sampling
  expect_equal(sampling$mode, "capped_per_node")
  expect_equal(sampling$original_events, nrow(fs$data))
  expect_equal(sampling$retained_events, fs$map$nNodes * 15L)
  expect_true(all(sampling$per_node$retained_events <= 15L))
  expect_equal(
    sampling$per_node$original_events,
    tabulate(fs$map$mapping[, 1], nbins = fs$map$nNodes)
  )

  set.seed(99)
  before <- runif(1)
  set.seed(99)
  suppressWarnings(invisible(iteration.QC(
    fs, ml, set.i,
    workers = 1L,
    uniform.test = "spread",
    max.events.per.node = 15L,
    seed = 42L,
    progress = FALSE
  )))
  after <- runif(1)
  expect_equal(before, after)
})

test_that("node caps and worker counts are validated at the QC boundary", {
  source_pkg_file("inflect-qc-core.R")

  for (cap in list(0L, NA_integer_, -1L, c(5L, 10L), 2.5)) {
    expect_error(
      .inflect_validate_max_events_per_node(cap),
      "`max.events.per.node`"
    )
  }
  expect_identical(.inflect_validate_max_events_per_node(1000L), 1000L)
  expect_error(.inflect_validate_workers(0L), "`workers`")
  expect_identical(.inflect_validate_workers(2L), 2L)
  windows <- .inflect_parallel_plan(
    workers = 2L,
    n_tasks = 10L,
    os_type = "windows"
  )
  expect_false(windows$use_parallel)
  expect_identical(windows$effective_workers, 1L)
  expect_identical(windows$backend, "serial_platform_fallback")
})

test_that("QC exposes neutral criterion-specific evidence and deprecated aliases", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(31L)
  set.i <- 3:6
  ml <- iteration.metacluster(fs, set.i)
  qc <- iteration.QC(
    fs,
    ml,
    set.i,
    uniform.test = "spread"
  )

  expect_named(
    qc,
    c(
      "scores",
      "dip_pass",
      "iqr_pass",
      "combined_pass",
      "criterion_pass",
      "qc.details",
      "provenance",
      "U.set",
      "Accuracy.matrixes"
    )
  )
  expect_named(qc$scores, c("k", "qc_pass_rate", "criterion"))
  expect_true(all(qc$scores$criterion == "iqr"))
  expect_identical(qc$criterion_pass, qc$iqr_pass)
  expect_identical(qc$Accuracy.matrixes, qc$criterion_pass)
  expect_equal(qc$U.set$Unimodality, qc$scores$qc_pass_rate)

  expected_fields <- c(
    "dip_pass",
    "iqr_pass",
    "combined_pass",
    "criterion_pass",
    "dip_p_value",
    "iqr",
    "event_count",
    "test_event_count",
    "failure_reason"
  )
  expect_named(qc$qc.details[[1]], expected_fields)
  expect_equal(dim(qc$qc.details[[1]]$dip_pass), c(set.i[[1]], 4L))
  expect_equal(
    qc$qc.details[[1]]$combined_pass,
    qc$qc.details[[1]]$dip_pass & qc$qc.details[[1]]$iqr_pass
  )
  expect_identical(qc$provenance$criterion, "iqr")
  expect_match(qc$provenance$criterion_label, "IQR")
  expect_equal(
    qc$provenance$criterion_summary$qc_pass_rate[
      qc$provenance$criterion_summary$selected
    ],
    qc$scores$qc_pass_rate
  )
})

test_that("complete finite transformed distributions are always retained", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")

  fs <- make_multinode_flowsom(32L)
  fs$data[1:20, "CD3"] <- seq(-5, 0, length.out = 20L)
  mc <- as.integer(rep(1:3, length.out = fs$map$nNodes))

  expect_no_warning(
    complete <- FlowSOMQC(fs, mc, uniform.test = "both", progress = FALSE)
  )
  complete_details <- attr(complete, "qc.details")
  expect_identical(
    complete_details$test_event_count,
    complete_details$event_count
  )
  expect_identical(
    attr(complete, "provenance")$value_handling,
    list(
      rule = "require finite QC marker data and retain all values unchanged",
      validation = "passed"
    )
  )
  expect_error(
    FlowSOMQC(fs, mc, zeroes.in = FALSE, progress = FALSE),
    "unused argument.*zeroes.in"
  )
})

test_that("iteration.QC validates requested names and exact cluster counts", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(33L)
  two <- as.integer(rep(1:2, length.out = fs$map$nNodes))
  three <- as.integer(rep(1:3, length.out = fs$map$nNodes))

  expect_error(
    iteration.QC(fs, list(two), set.i = 2L),
    "unique names"
  )
  expect_error(
    iteration.QC(
      fs,
      list("2" = two, "3" = three),
      set.i = 2L
    ),
    "unexpected: 3"
  )
  expect_error(
    iteration.QC(fs, list("2" = rep(1L, fs$map$nNodes)), set.i = 2L),
    "exactly 2 unique cluster labels"
  )
})

test_that("compiled IQR and p-value match the pure-R paths when available", {
  skip_if_not_installed("fastINFLECT")
  ns <- asNamespace("fastINFLECT")
  skip_if_not(exists(".inflect_iqr_cpp", where = ns, inherits = FALSE))
  iqr_cpp <- get(".inflect_iqr_cpp", ns)
  dp_cpp <- get(".inflect_dip_pvalue_cpp", ns)

  set.seed(8)
  for (i in 1:200) {
    x <- c(rnorm(sample(5:400, 1)), rexp(sample(0:40, 1)))
    q <- stats::quantile(x, names = FALSE)
    expect_identical(iqr_cpp(x), q[[4L]] - q[[2L]])
  }

  data("qDiptab", package = "diptest", envir = environment())
  nn <- as.integer(dimnames(qDiptab)[["n"]])
  Ps <- as.numeric(dimnames(qDiptab)[["Pr"]])
  set.seed(9)
  samples <- lapply(1:300, function(i) rnorm(sample(6:300, 1)))
  D <- vapply(samples, diptest::dip, numeric(1))
  n <- vapply(samples, length, integer(1))
  p_ref <- vapply(samples, function(x) suppressWarnings(diptest::dip.test(x))$p.value, numeric(1))
  p_cpp <- dp_cpp(D, as.integer(n), qDiptab, nn, Ps)
  ## compiled table interpolation must not change any pass/fail decision
  expect_identical(sum((p_ref >= 0.05) != (p_cpp >= 0.05)), 0L)

  ## NA and NaN statistics must both map to p = 1, exactly like R's is.na() guard
  expect_equal(dp_cpp(c(NA_real_, NaN), c(50L, 50L), qDiptab, nn, Ps), c(1, 1))
})
