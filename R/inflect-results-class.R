#' Create a fastINFLECT results object
#'
#' @description
#' `inflect.results` is an S3 list class containing the inflection point,
#' diagnostic curve, metaclustering results, and QC accuracy matrices returned by
#' [INFLECT()]. The class name remains `inflect.results` for compatibility with
#' the original INFLECT API.
#'
#' @param scores Canonical data frame containing `k`, `qc_pass_rate`, and
#'   `criterion`.
#' @param collection.U Deprecated score alias.
#' @param fittedcurve Data frame containing the coordinates of the fitted
#'   diagnostic curve.
#' @param lfunction List returned by [Lfunction()].
#' @param ggplot The diagnostic `ggplot2` object.
#' @param metaclustering.list List containing the iterative metaclustering
#'   results.
#' @param accuracy.sets List containing [FlowSOMQC()] results per
#'   metaclustering.
#' @param Accuracy.sets Backward-compatible alias for `accuracy.sets`.
#' @param qc.details Criterion-level matrices and numeric evidence by k.
#' @param dip_pass,iqr_pass,combined_pass,criterion_pass Named matrix lists.
#' @param provenance List containing normalized inputs and runtime metadata.
#'
#' @return An S3 object with class `inflect.results`.
#' @keywords internal
#' @noRd
new_inflect_results <- function(scores = NULL,
                                collection.U = NULL,
                                fittedcurve,
                                lfunction,
                                ggplot,
                                metaclustering.list = list(),
                                accuracy.sets = NULL,
                                Accuracy.sets = NULL,
                                qc.details = list(),
                                dip_pass = list(),
                                iqr_pass = list(),
                                combined_pass = list(),
                                criterion_pass = NULL,
                                selection = NULL,
                                provenance = list()) {
  normalize_scores <- get0(
    ".inflect_score_frame",
    mode = "function",
    inherits = TRUE
  )
  if (is.null(normalize_scores)) {
    normalize_scores <- function(x) {
      if (all(c("k", "qc_pass_rate") %in% names(x))) {
        return(x)
      }
      if (all(c("i", "Unimodality") %in% names(x))) {
        return(data.frame(
          k = x$i,
          qc_pass_rate = x$Unimodality,
          stringsAsFactors = FALSE
        ))
      }
      stop("Scores must contain `k` and `qc_pass_rate`.", call. = FALSE)
    }
  }
  if (is.null(scores)) {
    scores <- normalize_scores(collection.U)
  } else {
    scores <- normalize_scores(scores)
  }
  if (is.null(collection.U)) {
    collection.U <- data.frame(
      k = scores$k,
      qc_pass_rate = scores$qc_pass_rate,
      criterion = if ("criterion" %in% names(scores)) {
        scores$criterion
      } else {
        NA_character_
      },
      i = scores$k,
      Unimodality = scores$qc_pass_rate,
      stringsAsFactors = FALSE
    )
  }
  if (is.null(criterion_pass)) {
    criterion_pass <- accuracy.sets
  }
  if (is.null(accuracy.sets)) {
    accuracy.sets <- criterion_pass
  }
  if (is.null(accuracy.sets)) {
    accuracy.sets <- Accuracy.sets
  }
  if (is.null(accuracy.sets)) {
    accuracy.sets <- list()
  }
  if (!is.list(provenance)) {
    stop("`provenance` must be a list", call. = FALSE)
  }

  validate_inflect_results(structure(
    list(
      scores = scores,
      qc.scores = scores,
      collection.U = collection.U,
      fittedcurve = fittedcurve,
      lfunction = lfunction,
      ggplot = ggplot,
      metaclustering.list = metaclustering.list,
      Accuracy.sets = accuracy.sets,
      accuracy.sets = accuracy.sets,
      dip_pass = dip_pass,
      iqr_pass = iqr_pass,
      combined_pass = combined_pass,
      criterion_pass = accuracy.sets,
      qc.details = qc.details,
      selection = selection,
      provenance = provenance
    ),
    class = c("inflect.results", "list")
  ))
}

validate_inflect_results <- function(x) {
  if (!inherits(x, "inflect.results")) {
    stop("`x` must inherit from class 'inflect.results'", call. = FALSE)
  }

  required <- c(
    "collection.U",
    "scores",
    "qc.scores",
    "fittedcurve",
    "lfunction",
    "ggplot",
    "metaclustering.list",
    "Accuracy.sets",
    "accuracy.sets",
    "dip_pass",
    "iqr_pass",
    "combined_pass",
    "criterion_pass",
    "qc.details",
    "provenance"
  )
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0) {
    stop(
      "Missing inflect.results field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.data.frame(x$collection.U)) {
    stop("`collection.U` must be a data frame", call. = FALSE)
  }
  if (!is.data.frame(x$scores) ||
      !all(c("k", "qc_pass_rate") %in% names(x$scores))) {
    stop("`scores` must contain `k` and `qc_pass_rate`.", call. = FALSE)
  }
  if (!identical(x$scores, x$qc.scores)) {
    stop("`scores` and `qc.scores` must be identical.", call. = FALSE)
  }
  if (!is.data.frame(x$fittedcurve)) {
    stop("`fittedcurve` must be a data frame", call. = FALSE)
  }
  if (!is.list(x$lfunction)) {
    stop("`lfunction` must be a list", call. = FALSE)
  }
  if (!is.list(x$metaclustering.list)) {
    stop("`metaclustering.list` must be a list", call. = FALSE)
  }
  if (!is.list(x$accuracy.sets)) {
    stop("`accuracy.sets` must be a list", call. = FALSE)
  }
  for (field in c(
    "dip_pass",
    "iqr_pass",
    "combined_pass",
    "criterion_pass",
    "qc.details"
  )) {
    if (!is.list(x[[field]])) {
      stop("`", field, "` must be a list.", call. = FALSE)
    }
  }
  if (!is.list(x$provenance)) {
    stop("`provenance` must be a list", call. = FALSE)
  }
  if (!identical(x$Accuracy.sets, x$accuracy.sets)) {
    stop("`Accuracy.sets` and `accuracy.sets` must be identical", call. = FALSE)
  }

  x
}
