.inflect_scalar <- function(x, name) {
  value <- x[[name]]
  if (is.null(value) || length(value) == 0) {
    return(NA)
  }
  value[[1]]
}

.inflect_collection_scalar <- function(collection, column, fun) {
  if (!is.data.frame(collection) ||
      !column %in% names(collection) ||
      length(collection[[column]]) == 0) {
    return(NA)
  }
  values <- collection[[column]]
  if (all(is.na(values))) {
    return(NA)
  }
  fun(values, na.rm = TRUE)
}

.inflect_accuracy_sets <- function(x) {
  accuracy.sets <- x$accuracy.sets
  if (is.null(accuracy.sets)) {
    accuracy.sets <- x$Accuracy.sets
  }
  accuracy.sets
}

.inflect_marker_count <- function(x) {
  accuracy.sets <- .inflect_accuracy_sets(x)
  if (is.null(accuracy.sets) || length(accuracy.sets) == 0) {
    return(NA_integer_)
  }

  first_accuracy <- accuracy.sets[[1]]
  if (is.null(first_accuracy)) {
    return(NA_integer_)
  }

  marker_count <- ncol(first_accuracy)
  if (is.null(marker_count)) {
    return(NA_integer_)
  }

  as.integer(marker_count)
}

.inflect_provenance_scalar <- function(provenance, name) {
  if (is.null(provenance) || !is.list(provenance)) {
    return(NA)
  }
  .inflect_scalar(provenance, name)
}

#' Print candidate values from a fastINFLECT result
#'
#' Displays the tested range, marker count, core run settings, and each row of
#' `x$selection`. The status beside each candidate shows whether its partition
#' was directly tested.
#'
#' @param x An `inflect.results` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.inflect.results <- function(x, ...) {
  result_summary <- summary(x)

  cat("fastINFLECT result\n")
  cat("  knee: ", result_summary$knee, "\n", sep = "")
  cat("  range: ", result_summary$range, "\n", sep = "")
  cat("  angle: ", result_summary$angle, "\n", sep = "")
  cat("  tested k values: ", result_summary$n_points, "\n", sep = "")
  cat(
    "  tested k range: ",
    result_summary$min_i,
    "-",
    result_summary$max_i,
    "\n",
    sep = ""
  )

  if (!is.na(result_summary$n_markers)) {
    cat("  markers: ", result_summary$n_markers, "\n", sep = "")
  }

  if (is.list(x$provenance) && length(x$provenance) > 0) {
    cat(
      "  provenance: uniform.test=",
      .inflect_provenance_scalar(x$provenance, "uniform.test"),
      "\n",
      sep = ""
    )
  }

  if (is.data.frame(x$selection) && nrow(x$selection) > 0) {
    cat("  k estimates and tested thresholds:\n")
    for (r in seq_len(nrow(x$selection))) {
      k <- x$selection$k[r]
      pass_rate <- x$selection$qc_pass_rate_at_k[r]
      status <- x$selection$k_status[r]
      cat(
        "    ", format(x$selection$method[r], width = 10), " k=",
        if (is.na(k)) "NA" else round(k),
        if (is.na(pass_rate)) {
          ""
        } else {
          paste0(" (", round(pass_rate, 1), "% QC pass)")
        },
        if (is.na(status)) "" else paste0(" [", status, "]"),
        "\n",
        sep = ""
      )
    }
  }

  invisible(x)
}

#' Summarize a fastINFLECT result
#'
#' Returns one row with the tested range, QC pass-rate range, run settings, and
#' the inflection, Kneedle, and threshold candidates.
#'
#' @param object An `inflect.results` object.
#' @param ... Unused.
#'
#' @return A one-row `data.frame` of scan settings and candidate values.
#' @export
summary.inflect.results <- function(object, ...) {
  collection <- object$scores
  if (is.null(collection)) {
    collection <- .inflect_score_frame(object$collection.U)
  }
  n_points <- if (is.data.frame(collection)) {
    nrow(collection)
  } else {
    NA_integer_
  }

  result_summary <- data.frame(
    knee = .inflect_scalar(object$lfunction, "knee"),
    range = .inflect_scalar(object$lfunction, "range"),
    angle = .inflect_scalar(object$lfunction, "angle"),
    n_points = n_points,
    min_i = .inflect_collection_scalar(collection, "k", min),
    max_i = .inflect_collection_scalar(collection, "k", max),
    min_qc_pass_rate = .inflect_collection_scalar(collection, "qc_pass_rate", min),
    max_qc_pass_rate = .inflect_collection_scalar(collection, "qc_pass_rate", max),
    n_markers = .inflect_marker_count(object),
    stringsAsFactors = FALSE
  )

  if (is.list(object$provenance) && length(object$provenance) > 0) {
    provenance <- object$provenance
    result_summary$uniform.test <- .inflect_provenance_scalar(provenance, "uniform.test")
    result_summary$th.pvalue <- .inflect_provenance_scalar(provenance, "th.pvalue")
    result_summary$th.IQR <- .inflect_provenance_scalar(provenance, "th.IQR")
    result_summary$package_version <- .inflect_provenance_scalar(provenance, "package_version")
    result_summary$criterion <- .inflect_provenance_scalar(provenance, "criterion")
  }

  if (is.data.frame(object$selection)) {
    sel <- object$selection
    pick <- function(method) {
      hit <- sel$k[sel$method == method]
      if (length(hit) == 0) NA_real_ else hit[[1]]
    }
    result_summary$k_inflection <- pick("inflection")
    result_summary$k_kneedle <- pick("kneedle")
    result_summary$k_threshold <- pick("threshold")
  }

  result_summary$min_unimodality <- result_summary$min_qc_pass_rate
  result_summary$max_unimodality <- result_summary$max_qc_pass_rate
  result_summary
}

#' Plot the fastINFLECT QC pass-rate curve
#'
#' Returns the diagnostic plot stored during [INFLECT()]. Use the curve to see
#' how the selected criterion changes across the literal `set.i` schedule.
#'
#' @param x An `inflect.results` object.
#' @param ... Unused.
#'
#' @return The stored diagnostic `ggplot` object.
#' @export
plot.inflect.results <- function(x, ...) {
  x$ggplot
}

#' Extract tested fastINFLECT QC scores
#'
#' Returns one row per tested `k`, with the aggregate `qc_pass_rate` and the
#' criterion selected by `uniform.test`.
#'
#' @param x An `inflect.results` object.
#' @param row.names `row.names` passed to `as.data.frame`.
#' @param optional `optional` passed to `as.data.frame`.
#' @param ... Additional arguments passed to `as.data.frame`.
#'
#' @return Canonical `x$scores` as a `data.frame`.
#' @export
as.data.frame.inflect.results <- function(x,
                                          row.names = NULL,
                                          optional = FALSE,
                                          ...) {
  scores <- x$scores
  if (is.null(scores)) {
    scores <- .inflect_score_frame(x$collection.U)
  }
  as.data.frame(
    scores,
    row.names = row.names,
    optional = optional,
    ...
  )
}
