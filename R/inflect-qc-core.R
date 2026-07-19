## Internal fast QC core shared by FlowSOMQC() and iteration.QC().
##
## The two computational levers implemented here are:
##   1. `.inflect_dip_pvalue()` reproduces the table-based p-value of
##      `diptest::dip.test()` from the raw dip statistic returned by the
##      ~12x-cheaper `diptest::dip()`. The public p-value closure also mirrors
##      `dip.test()`'s incomplete-case handling before calling the fast path.
##   2. `.inflect_prepare_qc()` / `.inflect_accuracy_matrix()` factor the QC of a
##      single metaclustering so that `iteration.QC()` can memoise the accuracy
##      of each distinct dendrogram subtree (SOM-node set) across all k.
##
## Nothing here is exported; it is an implementation detail of the public API.

.inflect_env <- new.env(parent = emptyenv())

## Lazily load and cache diptest's tabulated null quantiles of the dip statistic.
.inflect_qdiptab <- function() {
  if (is.null(.inflect_env$qDiptab)) {
    e <- new.env()
    utils::data("qDiptab", package = "diptest", envir = e)
    .inflect_env$qDiptab <- e$qDiptab
  }
  .inflect_env$qDiptab
}

## TRUE when the package's compiled accelerators are linked and callable. False in
## sourced-file unit tests or an R-only install, where the pure-R paths are used.
.inflect_have_cpp <- function() {
  if (is.null(.inflect_env$have_cpp)) {
    .inflect_env$have_cpp <- tryCatch(
      is.function(get0(".inflect_iqr_cpp", mode = "function")) &&
        isTRUE(.inflect_iqr_cpp(c(1, 2, 3, 4, 5)) == 2),
      error = function(e) FALSE
    )
  }
  .inflect_env$have_cpp
}

## Type-7 inter-quartile range, via the O(n) compiled path when available and
## stats::quantile() otherwise. Missing values deliberately fall through to
## stats::quantile() so the legacy error behaviour is preserved.
.inflect_iqr <- function(x) {
  if (anyNA(x)) {
    qs <- stats::quantile(x, names = FALSE)
    return(qs[[4L]] - qs[[2L]])
  }
  if (.inflect_have_cpp()) {
    return(.inflect_iqr_cpp(x))
  }
  qs <- stats::quantile(x, names = FALSE)
  qs[[4L]] - qs[[2L]]
}

## Exact reproduction of the table branch of diptest::dip.test(): from the raw
## dip statistic `D` and sample size `n` it returns the interpolated p-value.
##
## For n <= 8, ties in the tabulated grid affect stats::approx()'s regularisation,
## so those few cases use the exact approx() route. Larger n use the compiled
## direct interpolation path when available; this avoids approx() setup overhead in
## the hot path while keeping p-values matched to dip.test().
.inflect_dip_pvalue_approx <- function(D, n, qd, nn, P.s) {
  max.n <- max(nn)
  if (is.na(D) || n <= 3L) {
    return(1)
  }
  if (n >= max.n) {
    n0 <- n1 <- max.n
    i.n <- i2 <- length(nn)
    f.n <- 0
  } else {
    i.n <- findInterval(n, nn)
    n0 <- nn[i.n]
    i2 <- i.n + 1L
    n1 <- nn[i2]
    f.n <- (n - n0) / (n1 - n0)
  }
  y.0 <- sqrt(n0) * qd[i.n, ]
  y.1 <- sqrt(n1) * qd[i2, ]
  sD <- sqrt(n) * D
  1 - stats::approx(y.0 + f.n * (y.1 - y.0), P.s, rule = 2, xout = sD)[["y"]]
}

.inflect_dip_pvalue <- function(D, n) {
  qd <- .inflect_qdiptab()
  nn <- as.integer(dimnames(qd)[["n"]])
  P.s <- as.numeric(dimnames(qd)[["Pr"]])
  max.n <- max(nn)
  L <- length(nn)
  M <- length(P.s)

  D <- as.numeric(D)
  n <- as.integer(n)
  if (length(n) == 1L && length(D) != 1L) {
    n <- rep.int(n, length(D))
  }
  if (length(n) != length(D)) {
    stop("`D` and `n` must have the same length, or `n` must be length 1.", call. = FALSE)
  }
  p <- numeric(length(D))
  simple <- is.na(D) | n <= 3L
  p[simple] <- 1

  exact <- !simple & n <= 8L
  if (any(exact)) {
    p[exact] <- vapply(which(exact), function(idx) {
      .inflect_dip_pvalue_approx(D[[idx]], n[[idx]], qd, nn, P.s)
    }, numeric(1))
  }

  fast <- !simple & !exact
  if (!any(fast)) {
    return(p)
  }
  if (.inflect_have_cpp()) {
    p[fast] <- .inflect_dip_pvalue_cpp(D[fast], n[fast], qd, nn, P.s)
    return(p)
  }

  for (idx in which(fast)) {
    ni <- n[idx]
    Di <- D[idx]
    if (ni >= max.n) {
      n0 <- n1 <- max.n
      i.n <- i2 <- L
      f.n <- 0
    } else {
      i.n <- findInterval(ni, nn)
      n0 <- nn[i.n]
      i2 <- i.n + 1L
      n1 <- nn[i2]
      f.n <- (ni - n0) / (n1 - n0)
    }
    y.0 <- sqrt(n0) * qd[i.n, ]
    grid <- y.0 + f.n * (sqrt(n1) * qd[i2, ] - y.0)
    sD <- sqrt(ni) * Di

    if (sD <= grid[1L]) {
      pv <- P.s[1L]
    } else if (sD >= grid[M]) {
      pv <- P.s[M]
    } else {
      j <- findInterval(sD, grid)
      gj <- grid[j]
      t <- (sD - gj) / (grid[j + 1L] - gj)
      pv <- P.s[j] + t * (P.s[j + 1L] - P.s[j])
    }
    p[idx] <- 1 - pv
  }
  p
}

## Apply the zeroes.in rule exactly as FlowSOMQC did: drop non-positive values and,
## if fewer than five positives remain, left-pad with zeros back up to five.
.inflect_marker_expression <- function(values, zeroes.in) {
  if (isFALSE(zeroes.in)) {
    me <- values[values > 0]
    if (length(me) < 5L) {
      me <- c(rep(0, 5L - length(me)), me)
    }
    return(me)
  }
  values
}

## Build the p-value closure. With simulate.p.value / B in `dots` we cannot use the
## table and fall back to diptest::dip.test()'s Monte-Carlo branch; otherwise we use
## the fast statistic-plus-table path.
.inflect_make_p_of <- function(dots = list()) {
  if (isTRUE(dots$simulate.p.value)) {
    force(dots)
    function(me) do.call(diptest::dip.test, c(list(me), dots))$p.value
  } else {
    function(me) {
      ## dip.test() drops incomplete cases before computing dip(); mirror that
      ## cleanup before the faster statistic-plus-table route.
      me <- me[stats::complete.cases(me)]
      .inflect_dip_pvalue(diptest::dip(me), length(me))
    }
  }
}

## Optional seeded subsampling cap that removes the sample-size confound of the dip
## test (large clusters make it over-powered). NULL leaves the sample untouched.
.inflect_maybe_subsample <- function(me, max.n, seed_key) {
  if (is.null(max.n) || length(me) <= max.n) {
    return(me)
  }
  ## deterministic per-cell RNG so the recommended k stays reproducible
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed_key)
  me[sort(sample.int(length(me), max.n))]
}

## Compute one accuracy row (logical, one entry per ordered marker) for the events
## in `expr` (a rows x markers numeric matrix). Reproduces the FlowSOMQC inner loop:
## clusters with <= 1 event yield all-NA; otherwise each marker passes when it is
## unimodal (dip) and/or tight (IQR) as selected by uniform.test.
.inflect_accuracy_row <- function(expr,
                                  zeroes.in,
                                  uniform.test,
                                  th.pvalue,
                                  th.IQR,
                                  p_of,
                                  subsample = NULL,
                                  seed = 1L) {
  markers <- colnames(expr)
  out <- stats::setNames(rep(NA, length(markers)), markers)
  if (nrow(expr) <= 1L) {
    return(out)
  }
  do_dip <- uniform.test != "spread"
  do_iqr <- uniform.test != "unimodality"
  for (j in seq_along(markers)) {
    me <- .inflect_marker_expression(expr[, j], zeroes.in)
    if (!is.null(subsample)) {
      me <- .inflect_maybe_subsample(me, subsample, seed + j)
    }
    uniform <- TRUE
    if (do_dip) {
      uniform <- p_of(me) >= th.pvalue
    }
    ## `uniform && (iqr < th)` is FALSE whenever `uniform` is already FALSE, so
    ## skipping the (sorting) IQR call in that case is exact, not an approximation.
    if (uniform && do_iqr) {
      uniform <- .inflect_iqr(me) < th.IQR
    }
    out[j] <- uniform
  }
  out
}

## Prepare the (unscaled, marker-ordered) event matrix and metadata once. Mirrors the
## data handling in FlowSOMQC so downstream QC works on identical inputs.
.inflect_prepare_qc <- function(view,
                                only.clustering.markers = TRUE,
                                acquired_markers = NULL) {
  data <- view$data
  if (isTRUE(view$scale)) {
    for (j in seq_len(ncol(data))) {
      data[, j] <- data[, j] * view$scaled.scale[j] + view$scaled.center[j]
    }
  }
  colnames(data) <- view$prettyColnames

  clustering.markers <- view$prettyColnames[view$map$colsUsed]
  if (only.clustering.markers) {
    markers <- clustering.markers
  } else {
    if (!is.null(acquired_markers) && all(acquired_markers %in% view$prettyColnames)) {
      markers <- acquired_markers
    } else {
      stop("Error in acquired_markers: The 'acquired_markers' vector must match names in 'FlowSOM.result$prettyColnames' ")
    }
  }

  ordered.markers <- c(
    gtools::mixedsort(intersect(markers, clustering.markers)),
    gtools::mixedsort(setdiff(markers, clustering.markers))
  )

  list(
    data = data[, ordered.markers, drop = FALSE],
    ordered.markers = ordered.markers,
    mapping = view$map$mapping[, 1]
  )
}

## Full accuracy matrix for one metaclustering, given a prepared QC context. Rows are
## clusters seq_len(max(metaclustering)); columns are the ordered markers. An optional
## `cache` environment memoises rows by the metacluster's sorted SOM-node set, so
## identical node sets appearing at different k are computed only once.
.inflect_accuracy_matrix <- function(prep,
                                     metaclustering,
                                     zeroes.in,
                                     uniform.test,
                                     th.pvalue,
                                     th.IQR,
                                     p_of,
                                     node_events = NULL,
                                     cache = NULL,
                                     subsample = NULL,
                                     seed = 1L,
                                     verbose = FALSE) {
  ordered.markers <- prep$ordered.markers
  clusters <- seq_len(max(metaclustering))
  accuracy.matrix <- matrix(
    nrow = length(clusters), ncol = length(ordered.markers),
    dimnames = list(clusters, ordered.markers)
  )

  event_cluster <- metaclustering[prep$mapping]
  cluster.rows <- split(seq_len(nrow(prep$data)), event_cluster)

  for (cl in clusters) {
    if (verbose) {
      message("Cluster: ", cl, " on ", length(clusters))
    }
    key <- NULL
    if (!is.null(cache)) {
      nodes <- which(metaclustering == cl)
      key <- paste0(nodes, collapse = ",")
      cached <- cache[[key]]
      if (!is.null(cached)) {
        accuracy.matrix[cl, ] <- cached
        next
      }
    }

    rows <- cluster.rows[[as.character(cl)]]
    expr <- if (is.null(rows)) {
      prep$data[integer(0), , drop = FALSE]
    } else {
      prep$data[rows, , drop = FALSE]
    }
    row <- .inflect_accuracy_row(
      expr = expr,
      zeroes.in = zeroes.in,
      uniform.test = uniform.test,
      th.pvalue = th.pvalue,
      th.IQR = th.IQR,
      p_of = p_of,
      subsample = subsample,
      seed = seed
    )
    accuracy.matrix[cl, ] <- row
    if (!is.null(cache)) {
      assign(key, row, envir = cache)
    }
  }
  accuracy.matrix
}
