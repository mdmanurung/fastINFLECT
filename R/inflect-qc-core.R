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
  1 - stats::approx(
    y.0 + f.n * (y.1 - y.0),
    P.s,
    rule = 2,
    ties = mean,
    xout = sD
  )[["y"]]
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

.inflect_validate_max_n_diptest <- function(max.n) {
  if (is.null(max.n)) {
    return(NULL)
  }
  if (length(max.n) != 1L ||
      !is.numeric(max.n) ||
      is.na(max.n) ||
      !is.finite(max.n) ||
      max.n < 4L ||
      max.n > .Machine$integer.max ||
      max.n != floor(max.n)) {
    stop("`max.n.diptest` must be NULL or a single integer >= 4.", call. = FALSE)
  }
  as.integer(max.n)
}

.inflect_validate_max_events_per_node <- function(max.events) {
  if (is.null(max.events)) {
    return(NULL)
  }
  if (length(max.events) != 1L ||
      !is.numeric(max.events) ||
      is.na(max.events) ||
      !is.finite(max.events) ||
      max.events < 1L ||
      max.events > .Machine$integer.max ||
      max.events != floor(max.events)) {
    stop(
      "`max.events.per.node` must be NULL or a single positive integer.",
      call. = FALSE
    )
  }
  as.integer(max.events)
}

.inflect_validate_seed <- function(seed) {
  if (length(seed) != 1L ||
      !is.numeric(seed) ||
      is.na(seed) ||
      !is.finite(seed) ||
      seed < 0 ||
      seed > .Machine$integer.max ||
      seed != floor(seed)) {
    stop(
      "`seed` must be a single non-negative integer no greater than `.Machine$integer.max`.",
      call. = FALSE
    )
  }
  as.integer(seed)
}

.inflect_validate_qc_arguments <- function(zeroes.in,
                                           th.pvalue,
                                           th.IQR) {
  if (length(zeroes.in) != 1L ||
      !is.logical(zeroes.in) ||
      is.na(zeroes.in)) {
    stop("`zeroes.in` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(th.pvalue) != 1L ||
      !is.numeric(th.pvalue) ||
      is.na(th.pvalue) ||
      !is.finite(th.pvalue) ||
      th.pvalue < 0 ||
      th.pvalue > 1) {
    stop("`th.pvalue` must be a finite number in [0, 1].", call. = FALSE)
  }
  if (length(th.IQR) != 1L ||
      !is.numeric(th.IQR) ||
      is.na(th.IQR) ||
      !is.finite(th.IQR) ||
      th.IQR < 0) {
    stop("`th.IQR` must be a finite non-negative number.", call. = FALSE)
  }
  invisible(TRUE)
}

.inflect_criterion <- function(uniform.test) {
  switch(
    uniform.test,
    both = "combined",
    spread = "iqr",
    unimodality = "dip",
    stop("Unknown QC criterion: ", uniform.test, call. = FALSE)
  )
}

.inflect_criterion_label <- function(criterion) {
  switch(
    criterion,
    combined = "dip and IQR QC pass rate",
    dip = "dip-test QC pass rate",
    iqr = "IQR-spread QC pass rate",
    stop("Unknown QC criterion: ", criterion, call. = FALSE)
  )
}

.inflect_resolve_cores <- function(multicore, cores) {
  if (length(multicore) != 1L || !is.logical(multicore) || is.na(multicore)) {
    stop("`multicore` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!isTRUE(multicore)) {
    return(cores)
  }
  if (is.null(cores)) {
    detected <- parallel::detectCores()
    if (length(detected) != 1L ||
        is.na(detected) ||
        !is.finite(detected) ||
        detected < 3L) {
      stop(
        "Could not resolve at least two workers; supply `cores >= 2` explicitly.",
        call. = FALSE
      )
    }
    cores <- detected - 1
  }
  if (length(cores) != 1L ||
      !is.numeric(cores) ||
      is.na(cores) ||
      !is.finite(cores) ||
      cores < 2L ||
      cores > .Machine$integer.max ||
      cores != floor(cores)) {
    stop(
      "`cores` must be a single integer >= 2 when `multicore = TRUE`.",
      call. = FALSE
    )
  }
  as.integer(cores)
}

.inflect_parallel_plan <- function(multicore,
                                   cores,
                                   n_tasks,
                                   os_type = .Platform$OS.type) {
  cores <- .inflect_resolve_cores(multicore, cores)
  if (length(n_tasks) != 1L ||
      !is.numeric(n_tasks) ||
      is.na(n_tasks) ||
      !is.finite(n_tasks) ||
      n_tasks < 0 ||
      n_tasks != floor(n_tasks)) {
    stop("`n_tasks` must be a non-negative integer.", call. = FALSE)
  }
  use_parallel <- isTRUE(multicore) &&
    n_tasks > 1L &&
    identical(os_type, "unix")
  effective_workers <- if (use_parallel) {
    min(as.integer(cores), as.integer(n_tasks))
  } else {
    1L
  }
  backend <- if (use_parallel) {
    "parallel::mclapply"
  } else if (isTRUE(multicore) && !identical(os_type, "unix")) {
    "serial_platform_fallback"
  } else {
    "serial"
  }
  list(
    use_parallel = use_parallel,
    requested_multicore = isTRUE(multicore),
    requested_cores = if (is.null(cores)) NA_integer_ else as.integer(cores),
    effective_workers = effective_workers,
    backend = backend,
    os_type = os_type
  )
}

.inflect_normalize_metaclustering <- function(metaclustering,
                                              n_nodes,
                                              label = "`metaclustering`",
                                              expected_k = NULL) {
  if (!is.integer(metaclustering)) {
    stop(label, " must be an integer vector", call. = FALSE)
  }
  if (length(metaclustering) != n_nodes) {
    stop(label, " must contain one entry per SOM node", call. = FALSE)
  }
  if (anyNA(metaclustering)) {
    stop(label, " must not contain missing values", call. = FALSE)
  }
  if (any(metaclustering < 1L)) {
    stop(label, " labels must be positive integers", call. = FALSE)
  }

  labels <- sort(unique(metaclustering))
  if (!is.null(expected_k) &&
      length(labels) != as.integer(expected_k)) {
    stop(
      label,
      " must contain exactly ",
      as.integer(expected_k),
      " unique cluster labels; found ",
      length(labels),
      ".",
      call. = FALSE
    )
  }
  if (!identical(labels, seq_len(max(labels)))) {
    metaclustering <- match(metaclustering, labels)
  }
  as.integer(metaclustering)
}

## Apply the requested zero rule literally. The default retains the complete
## transformed distribution. Compatibility mode excludes every non-positive
## value and never manufactures padding observations.
.inflect_marker_expression <- function(values, zeroes.in) {
  if (isFALSE(zeroes.in)) {
    return(values[values > 0])
  }
  values
}

## Build the p-value closure. With simulate.p.value / B in `dots` we cannot use the
## table and fall back to diptest::dip.test()'s Monte-Carlo branch; otherwise we use
## the fast statistic-plus-table path.
.inflect_with_seed <- function(seed, code) {
  seed <- .inflect_validate_seed(seed)
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
  set.seed(seed)
  force(code)
}

.inflect_make_p_of <- function(dots = list()) {
  if (isTRUE(dots$simulate.p.value)) {
    force(dots)
    function(me, seed = 1L) {
      .inflect_with_seed(
        seed,
        do.call(diptest::dip.test, c(list(me), dots))$p.value
      )
    }
  } else {
    function(me, seed = NULL) {
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

.inflect_seed_from_key <- function(seed, ..., stream = "") {
  seed <- .inflect_validate_seed(seed)
  modulus <- as.double(.Machine$integer.max)
  key <- enc2utf8(paste(c(..., stream), collapse = "\u001f"))
  bytes <- as.integer(charToRaw(key))
  hash <- as.double(seed) %% modulus
  if (length(bytes) > 0L) {
    for (byte in bytes) {
      ## The intermediate stays far below 2^53, so the modular arithmetic is
      ## exact and cannot overflow R's integer representation.
      hash <- (hash * 131 + byte + 1) %% modulus
    }
  }
  as.integer(hash)
}

.inflect_seed_for <- function(seed, index) {
  .inflect_seed_from_key(seed, "index", format(index, scientific = FALSE))
}

## Sample each SOM node once before scoring. The retained row indices are shared
## by every subtree, marker, test mode, and worker.
.inflect_sample_node_events <- function(node_events,
                                        max.events.per.node = NULL,
                                        seed = 1L) {
  max.events.per.node <- .inflect_validate_max_events_per_node(max.events.per.node)
  seed <- .inflect_validate_seed(seed)
  original_counts <- lengths(node_events)
  sampled <- node_events

  if (!is.null(max.events.per.node)) {
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

    for (node in seq_along(sampled)) {
      rows <- sampled[[node]]
      if (length(rows) > max.events.per.node) {
        set.seed(.inflect_seed_from_key(
          seed,
          "node",
          node,
          stream = "event_sample"
        ))
        keep <- sort(sample.int(length(rows), max.events.per.node))
        sampled[[node]] <- rows[keep]
      }
    }
  }

  retained_counts <- lengths(sampled)
  list(
    node_events = sampled,
    provenance = list(
      mode = if (is.null(max.events.per.node)) "all_events" else "capped_per_node",
      max_events_per_node = if (is.null(max.events.per.node)) {
        NA_integer_
      } else {
        max.events.per.node
      },
      original_events = sum(original_counts),
      retained_events = sum(retained_counts),
      per_node = data.frame(
        node = seq_along(node_events),
        original_events = as.integer(original_counts),
        retained_events = as.integer(retained_counts)
      )
    )
  )
}

.inflect_validate_worker_results <- function(rows_list,
                                             keys,
                                             marker_names,
                                             parallel = FALSE) {
  failures <- character(0)
  expected <- length(keys)
  if (!is.list(rows_list) || length(rows_list) != expected) {
    actual <- if (is.list(rows_list)) length(rows_list) else 0L
    failures <- c(
      failures,
      paste0("result count ", actual, " (expected ", expected, ")")
    )
    if (actual < expected) {
      for (idx in seq.int(actual + 1L, expected)) {
        failures <- c(
          failures,
          paste0("task ", idx, " (subtree ", keys[[idx]], "): missing result")
        )
      }
    } else if (actual > expected) {
      failures <- c(
        failures,
        paste0(actual - expected, " unexpected extra result(s)")
      )
    }
  }

  checked <- min(if (is.list(rows_list)) length(rows_list) else 0L, expected)
  detail_schema <- c(
    dip_pass = "logical",
    iqr_pass = "logical",
    combined_pass = "logical",
    criterion_pass = "logical",
    dip_p_value = "double",
    iqr = "double",
    event_count = "integer",
    test_event_count = "integer",
    excluded_nonpositive = "integer",
    excluded_nonfinite = "integer",
    failure_reason = "character"
  )
  for (idx in seq_len(checked)) {
    result <- rows_list[[idx]]
    reason <- NULL
    if (is.null(result)) {
      reason <- "NULL result"
    } else if (inherits(result, "try-error")) {
      reason <- paste0("try-error: ", trimws(as.character(result)[[1L]]))
    } else if (is.logical(result)) {
      ## Retain validation support for the deprecated internal Boolean-row
      ## schema so failures from older extension code remain intelligible.
      if (length(result) != length(marker_names)) {
        reason <- paste0(
          "length ",
          length(result),
          " (expected ",
          length(marker_names),
          ")"
        )
      } else if (!identical(names(result), marker_names)) {
        reason <- "marker names do not match the scoring schedule"
      }
    } else if (!is.list(result)) {
      reason <- paste0("type ", typeof(result), " (expected a QC detail list)")
    } else {
      missing_fields <- setdiff(names(detail_schema), names(result))
      if (length(missing_fields) > 0L) {
        reason <- paste0(
          "missing field(s): ",
          paste(missing_fields, collapse = ", ")
        )
      } else {
        field_failures <- character(0)
        for (field in names(detail_schema)) {
          value <- result[[field]]
          expected_type <- detail_schema[[field]]
          if (typeof(value) != expected_type) {
            field_failures <- c(
              field_failures,
              paste0(field, " has type ", typeof(value), " (expected ", expected_type, ")")
            )
          } else if (length(value) != length(marker_names)) {
            field_failures <- c(
              field_failures,
              paste0(field, " has length ", length(value), " (expected ", length(marker_names), ")")
            )
          } else if (!identical(names(value), marker_names)) {
            field_failures <- c(
              field_failures,
              paste0(field, " marker names do not match")
            )
          }
        }
        if (length(field_failures) > 0L) {
          reason <- paste(field_failures, collapse = ", ")
        }
      }
    }
    if (!is.null(reason)) {
      failures <- c(
        failures,
        paste0("task ", idx, " (subtree ", keys[[idx]], "): ", reason)
      )
    }
  }

  if (length(failures) > 0L) {
    prefix <- if (isTRUE(parallel)) {
      "Parallel QC worker failure"
    } else {
      "QC scoring failure"
    }
    stop(prefix, ": ", paste(failures, collapse = "; "), call. = FALSE)
  }
  rows_list
}

.inflect_named <- function(value, markers, mode) {
  out <- rep(value, length(markers))
  storage.mode(out) <- mode
  stats::setNames(out, markers)
}

.inflect_empty_qc_row <- function(markers) {
  list(
    dip_pass = .inflect_named(NA, markers, "logical"),
    iqr_pass = .inflect_named(NA, markers, "logical"),
    combined_pass = .inflect_named(NA, markers, "logical"),
    criterion_pass = .inflect_named(NA, markers, "logical"),
    dip_p_value = .inflect_named(NA_real_, markers, "double"),
    iqr = .inflect_named(NA_real_, markers, "double"),
    event_count = .inflect_named(0L, markers, "integer"),
    test_event_count = .inflect_named(0L, markers, "integer"),
    excluded_nonpositive = .inflect_named(0L, markers, "integer"),
    excluded_nonfinite = .inflect_named(0L, markers, "integer"),
    failure_reason = .inflect_named(NA_character_, markers, "character")
  )
}

.inflect_append_failure <- function(current, reason) {
  if (is.na(current) || !nzchar(current)) {
    return(reason)
  }
  paste(current, reason, sep = "; ")
}

.inflect_call_p_of <- function(p_of, me, seed) {
  p_formals <- names(formals(p_of))
  if ("seed" %in% p_formals || "..." %in% p_formals) {
    p_of(me, seed = seed)
  } else {
    p_of(me)
  }
}

## Compute the complete criterion-level evidence for one subtree. Both tests are
## evaluated regardless of the selected aggregate criterion so callers can tell
## an IQR pass from a dip pass and can audit discordant cells.
.inflect_qc_row_indexed <- function(data,
                                    rows,
                                    zeroes.in,
                                    uniform.test,
                                    th.pvalue,
                                    th.IQR,
                                    p_of,
                                    subsample = NULL,
                                    seed = 1L,
                                    subtree_key = "",
                                    marker_indices = seq_len(ncol(data)),
                                    marker_names = colnames(data)[marker_indices]) {
  markers <- marker_names
  out <- .inflect_empty_qc_row(markers)
  criterion <- .inflect_criterion(uniform.test)
  if (length(rows) == 0L) {
    out$failure_reason[] <- "no_events"
    return(out)
  }

  for (j in seq_along(markers)) {
    marker <- markers[[j]]
    marker_index <- marker_indices[[j]]
    raw_values <- data[rows, marker_index]
    finite <- is.finite(raw_values)
    values <- raw_values[finite]
    excluded_nonfinite <- sum(!finite)
    excluded_nonpositive <- if (isFALSE(zeroes.in)) {
      sum(values <= 0)
    } else {
      0L
    }
    me <- .inflect_marker_expression(values, zeroes.in)

    out$event_count[[j]] <- as.integer(length(raw_values))
    out$test_event_count[[j]] <- as.integer(length(me))
    out$excluded_nonpositive[[j]] <- as.integer(excluded_nonpositive)
    out$excluded_nonfinite[[j]] <- as.integer(excluded_nonfinite)

    if (length(me) <= 1L) {
      out$failure_reason[[j]] <- if (length(me) == 0L) {
        "no_testable_values"
      } else {
        "only_one_testable_value"
      }
      next
    }

    dip_values <- if (!is.null(subsample)) {
      .inflect_maybe_subsample(
        me,
        subsample,
        .inflect_seed_from_key(
          seed,
          subtree_key,
          marker,
          stream = "dip_subsample"
        )
      )
    } else {
      me
    }
    dip_result <- tryCatch(
      {
        p <- .inflect_call_p_of(
          p_of,
          dip_values,
          .inflect_seed_from_key(
            seed,
            subtree_key,
            marker,
            stream = "dip_simulation"
          )
        )
        if (length(p) != 1L || !is.numeric(p) || !is.finite(p) ||
            p < 0 || p > 1) {
          stop("dip test returned an invalid p-value")
        }
        list(value = as.numeric(p), error = NULL)
      },
      error = function(e) list(value = NA_real_, error = conditionMessage(e))
    )
    out$dip_p_value[[j]] <- dip_result$value
    if (is.null(dip_result$error)) {
      out$dip_pass[[j]] <- dip_result$value >= th.pvalue
    } else {
      out$failure_reason[[j]] <- .inflect_append_failure(
        out$failure_reason[[j]],
        paste0("dip_error: ", dip_result$error)
      )
    }

    iqr_result <- tryCatch(
      {
        value <- .inflect_iqr(me)
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
          stop("IQR calculation returned an invalid value")
        }
        list(value = as.numeric(value), error = NULL)
      },
      error = function(e) list(value = NA_real_, error = conditionMessage(e))
    )
    out$iqr[[j]] <- iqr_result$value
    if (is.null(iqr_result$error)) {
      out$iqr_pass[[j]] <- iqr_result$value < th.IQR
    } else {
      out$failure_reason[[j]] <- .inflect_append_failure(
        out$failure_reason[[j]],
        paste0("iqr_error: ", iqr_result$error)
      )
    }

    out$combined_pass[[j]] <- out$dip_pass[[j]] & out$iqr_pass[[j]]
    out$criterion_pass[[j]] <- switch(
      criterion,
      combined = out$combined_pass[[j]],
      dip = out$dip_pass[[j]],
      iqr = out$iqr_pass[[j]]
    )
  }
  out
}

## Deprecated Boolean-row wrapper retained for exact comparison tests. New
## code should consume `.inflect_qc_row_indexed()`.
.inflect_accuracy_row <- function(expr,
                                  zeroes.in,
                                  uniform.test,
                                  th.pvalue,
                                  th.IQR,
                                  p_of,
                                  subsample = NULL,
                                  seed = 1L) {
  .inflect_accuracy_row_indexed(
    data = expr,
    rows = seq_len(nrow(expr)),
    zeroes.in = zeroes.in,
    uniform.test = uniform.test,
    th.pvalue = th.pvalue,
    th.IQR = th.IQR,
    p_of = p_of,
    subsample = subsample,
    seed = seed
  )
}

## Indexed equivalent of `.inflect_accuracy_row()`. Only one marker vector is
## materialised at a time, avoiding full events-by-markers subtree copies.
.inflect_accuracy_row_indexed <- function(data,
                                          rows,
                                          zeroes.in,
                                          uniform.test,
                                          th.pvalue,
                                          th.IQR,
                                          p_of,
                                          subsample = NULL,
                                          seed = 1L,
                                          marker_indices = seq_len(ncol(data)),
                                          marker_names = colnames(data)[marker_indices]) {
  .inflect_qc_row_indexed(
    data = data,
    rows = rows,
    zeroes.in = zeroes.in,
    uniform.test = uniform.test,
    th.pvalue = th.pvalue,
    th.IQR = th.IQR,
    p_of = p_of,
    subsample = subsample,
    seed = seed,
    marker_indices = marker_indices,
    marker_names = marker_names
  )$criterion_pass
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
    data = data,
    marker_indices = match(ordered.markers, view$prettyColnames),
    ordered.markers = ordered.markers,
    mapping = view$map$mapping[, 1]
  )
}

.inflect_validate_qc_schedule <- function(set.i, n_nodes) {
  if (!is.numeric(set.i) ||
      length(set.i) == 0L ||
      anyNA(set.i) ||
      any(!is.finite(set.i)) ||
      any(set.i != floor(set.i))) {
    stop(
      "`set.i` must contain finite integer cluster counts.",
      call. = FALSE
    )
  }
  if (anyDuplicated(set.i)) {
    stop("`set.i` must contain unique cluster counts.", call. = FALSE)
  }
  if (is.unsorted(set.i, strictly = TRUE)) {
    stop("`set.i` must be strictly increasing.", call. = FALSE)
  }
  if (any(set.i < 1L) || any(set.i > n_nodes)) {
    stop(
      "`set.i` values must lie between 1 and the number of SOM nodes.",
      call. = FALSE
    )
  }
  as.integer(set.i)
}

.inflect_validate_metaclustering_list <- function(metaclustering.list,
                                                  set.i,
                                                  n_nodes) {
  if (!is.list(metaclustering.list)) {
    stop("`metaclustering.list` must be a named list.", call. = FALSE)
  }
  requested_names <- as.character(set.i)
  actual_names <- names(metaclustering.list)
  if (is.null(actual_names) ||
      anyNA(actual_names) ||
      any(!nzchar(actual_names)) ||
      anyDuplicated(actual_names)) {
    stop(
      "`metaclustering.list` must have unique names matching every requested k.",
      call. = FALSE
    )
  }
  missing_names <- setdiff(requested_names, actual_names)
  extra_names <- setdiff(actual_names, requested_names)
  if (length(missing_names) > 0L || length(extra_names) > 0L) {
    details <- c(
      if (length(missing_names)) {
        paste0("missing: ", paste(missing_names, collapse = ", "))
      },
      if (length(extra_names)) {
        paste0("unexpected: ", paste(extra_names, collapse = ", "))
      }
    )
    stop(
      "`metaclustering.list` must match `set.i` exactly (",
      paste(details, collapse = "; "),
      ").",
      call. = FALSE
    )
  }

  normalized <- vector("list", length(set.i))
  names(normalized) <- requested_names
  for (idx in seq_along(set.i)) {
    k <- set.i[[idx]]
    normalized[[idx]] <- .inflect_normalize_metaclustering(
      metaclustering = metaclustering.list[[as.character(k)]],
      n_nodes = n_nodes,
      label = paste0("`metaclustering.list[['", k, "']]`"),
      expected_k = k
    )
  }
  normalized
}

.inflect_zero_handling <- function(prep, zeroes.in, warn = TRUE) {
  rows <- lapply(seq_along(prep$ordered.markers), function(j) {
    values <- prep$data[, prep$marker_indices[[j]]]
    finite <- is.finite(values)
    negative <- sum(values[finite] < 0)
    zero <- sum(values[finite] == 0)
    data.frame(
      marker = prep$ordered.markers[[j]],
      total_events = length(values),
      finite_events = sum(finite),
      negative_values = negative,
      zero_values = zero,
      nonpositive_values = negative + zero,
      excluded_nonpositive = if (isFALSE(zeroes.in)) {
        negative + zero
      } else {
        0
      },
      excluded_nonfinite = sum(!finite),
      stringsAsFactors = FALSE
    )
  })
  per_marker <- do.call(rbind, rows)
  per_marker$excluded_fraction <- if (isFALSE(zeroes.in)) {
    per_marker$excluded_nonpositive / pmax(1, per_marker$finite_events)
  } else {
    0
  }

  negative_markers <- per_marker$marker[per_marker$negative_values > 0]
  if (isFALSE(zeroes.in) &&
      isTRUE(warn) &&
      length(negative_markers) > 0L) {
    exclusion_range <- range(
      100 * per_marker$excluded_fraction[per_marker$negative_values > 0]
    )
    warning(
      "`zeroes.in = FALSE` excludes every non-positive value. Negative ",
      "transformed values were found in ",
      length(negative_markers),
      " marker(s); their per-marker exclusion fractions range from ",
      format(round(exclusion_range[[1]], 2), trim = TRUE),
      "% to ",
      format(round(exclusion_range[[2]], 2), trim = TRUE),
      "%. Use `zeroes.in = TRUE` to score complete transformed distributions.",
      call. = FALSE
    )
  }
  list(
    zeroes_in = isTRUE(zeroes.in),
    rule = if (isTRUE(zeroes.in)) {
      "retain finite negative, zero, and positive transformed values"
    } else {
      "exclude all non-positive transformed values"
    },
    per_marker = per_marker
  )
}

.inflect_qc_rows_to_detail <- function(rows,
                                       marker_names,
                                       cluster_names) {
  fields <- names(.inflect_empty_qc_row(marker_names))
  detail <- lapply(fields, function(field) {
    matrix <- do.call(rbind, lapply(rows, `[[`, field))
    dimnames(matrix) <- list(cluster_names, marker_names)
    matrix
  })
  names(detail) <- fields
  detail
}

.inflect_qc_summary <- function(details, selected_criterion) {
  fields <- c(
    dip = "dip_pass",
    iqr = "iqr_pass",
    combined = "combined_pass"
  )
  rows <- list()
  row_index <- 0L
  for (k in names(details)) {
    for (criterion in names(fields)) {
      row_index <- row_index + 1L
      values <- details[[k]][[fields[[criterion]]]]
      total <- length(values)
      rows[[row_index]] <- data.frame(
        k = as.integer(k),
        criterion = criterion,
        selected = identical(criterion, selected_criterion),
        passed = sum(values %in% TRUE),
        failed = sum(values %in% FALSE),
        unresolved = sum(is.na(values)),
        total = total,
        qc_pass_rate = if (total == 0L) {
          NA_real_
        } else {
          sum(values %in% TRUE) * 100 / total
        },
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
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
    if (is.null(rows)) {
      rows <- integer(0)
    }
    row <- .inflect_accuracy_row_indexed(
      data = prep$data,
      rows = rows,
      zeroes.in = zeroes.in,
      uniform.test = uniform.test,
      th.pvalue = th.pvalue,
      th.IQR = th.IQR,
      p_of = p_of,
      subsample = subsample,
      seed = seed,
      marker_indices = prep$marker_indices,
      marker_names = prep$ordered.markers
    )
    accuracy.matrix[cl, ] <- row
    if (!is.null(cache)) {
      assign(key, row, envir = cache)
    }
  }
  accuracy.matrix
}
