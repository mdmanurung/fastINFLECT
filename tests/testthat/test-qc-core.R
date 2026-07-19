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

test_that("memoised iteration.QC reproduces per-k FlowSOMQC exactly", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")
  source_pkg_file("metaClustering-hclust.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom()
  set.i <- 3:10
  ml <- iteration.metacluster(fs, set.i = set.i, multicore = FALSE)

  reference <- lapply(set.i, function(k) {
    FlowSOMQC(fs, ml[[as.character(k)]], zeroes.in = FALSE, verbose = FALSE)
  })
  ref_scores <- vapply(reference, function(m) sum(m, na.rm = TRUE) * 100 / prod(dim(m)), numeric(1))

  qc <- iteration.QC(fs, ml, set.i, multicore = FALSE, zeroes.in = FALSE, verbose = FALSE)

  expect_equal(qc$U.set$Unimodality, ref_scores)
  for (idx in seq_along(set.i)) {
    expect_identical(unname(qc$Accuracy.matrixes[[idx]]), unname(reference[[idx]]))
  }
})

test_that("memoisation holds across uniform.test modes and zeroes.in", {
  source_pkg_file("som-adapter.R")
  source_pkg_file("inflect-qc-core.R")
  source_pkg_file("FlowSOM-QC.R")
  source_pkg_file("metaClustering-hclust.R")
  source_pkg_file("iteration-metacluster.R")
  source_pkg_file("iteration-QC.R")

  fs <- make_multinode_flowsom(2L)
  set.i <- 3:8
  ml <- iteration.metacluster(fs, set.i = set.i, multicore = FALSE)

  for (ut in c("both", "spread", "unimodality")) {
    for (zi in c(TRUE, FALSE)) {
      ref <- lapply(set.i, function(k) {
        FlowSOMQC(fs, ml[[as.character(k)]], zeroes.in = zi, uniform.test = ut, verbose = FALSE)
      })
      qc <- iteration.QC(fs, ml, set.i, multicore = FALSE, zeroes.in = zi,
                         uniform.test = ut, verbose = FALSE)
      ok <- all(mapply(function(a, b) identical(unname(a), unname(b)),
                       qc$Accuracy.matrixes, ref))
      expect_true(ok, info = paste(ut, zi))
    }
  }
})

test_that("fast dip p-value matches diptest::dip.test decisions", {
  skip_if_not_installed("diptest")
  source_pkg_file("inflect-qc-core.R")

  set.seed(42)
  disagreements <- 0L
  for (n in c(6, 20, 60, 250)) {
    for (gen in list(function(m) rnorm(m),
                     function(m) c(rnorm(m %/% 2), rnorm(m - m %/% 2, 6)),
                     function(m) rexp(m))) {
      for (rep in 1:6) {
        x <- gen(n)
        p_ref <- suppressWarnings(diptest::dip.test(x))$p.value
        p_fast <- .inflect_dip_pvalue(diptest::dip(x), length(x))
        disagreements <- disagreements + as.integer((p_ref >= 0.05) != (p_fast >= 0.05))
      }
    }
  }
  expect_identical(disagreements, 0L)
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
  ml <- iteration.metacluster(fs, set.i = set.i, multicore = FALSE)

  a <- iteration.QC(fs, ml, set.i, multicore = FALSE, zeroes.in = FALSE,
                    max.n.diptest = 30L, seed = 7L, verbose = FALSE)
  b <- iteration.QC(fs, ml, set.i, multicore = FALSE, zeroes.in = FALSE,
                    max.n.diptest = 30L, seed = 7L, verbose = FALSE)
  expect_equal(a$U.set, b$U.set)

  ## the RNG state must be restored: an outer sequence is unaffected
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  invisible(iteration.QC(fs, ml, set.i, multicore = FALSE, zeroes.in = FALSE,
                         max.n.diptest = 30L, seed = 1L, verbose = FALSE))
  after <- runif(1)
  expect_equal(before, after)
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
