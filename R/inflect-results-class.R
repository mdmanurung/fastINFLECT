#' Create a fastINFLECT results object
#'
#' @description
#' `inflect.results` is an S3 list class containing the inflection point,
#' diagnostic curve, metaclustering results, and QC accuracy matrices returned by
#' [INFLECT()]. The class name remains `inflect.results` for compatibility with
#' the original INFLECT API.
#'
#' @param collection.U Data frame containing the Unimodality scores for each
#'   `i` in `set.i`.
#' @param fittedcurve Data frame containing the coordinates of the fitted
#'   diagnostic curve.
#' @param lfunction List returned by [Lfunction()].
#' @param ggplot The diagnostic `ggplot2` object.
#' @param metaclustering.list List containing the iterative metaclustering
#'   results.
#' @param accuracy.sets List containing [FlowSOMQC()] results per
#'   metaclustering.
#' @param Accuracy.sets Backward-compatible alias for `accuracy.sets`.
#' @param provenance List containing normalized inputs and runtime metadata.
#'
#' @return An S3 object with class `inflect.results`.
#' @keywords internal
#' @noRd
new_inflect_results <- function(collection.U,
                                fittedcurve,
                                lfunction,
                                ggplot,
                                metaclustering.list = list(),
                                accuracy.sets = NULL,
                                Accuracy.sets = NULL,
                                selection = NULL,
                                provenance = list()) {
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
      collection.U = collection.U,
      fittedcurve = fittedcurve,
      lfunction = lfunction,
      ggplot = ggplot,
      metaclustering.list = metaclustering.list,
      Accuracy.sets = accuracy.sets,
      accuracy.sets = accuracy.sets,
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
    "fittedcurve",
    "lfunction",
    "ggplot",
    "metaclustering.list",
    "Accuracy.sets",
    "accuracy.sets",
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
  if (!is.list(x$provenance)) {
    stop("`provenance` must be a list", call. = FALSE)
  }
  if (!identical(x$Accuracy.sets, x$accuracy.sets)) {
    stop("`Accuracy.sets` and `accuracy.sets` must be identical", call. = FALSE)
  }

  x
}
