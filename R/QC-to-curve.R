#' Plot a criterion-specific fastINFLECT QC curve
#'
#' @description
#' Fits the original four-parameter log-logistic diagnostic curve to
#' criterion-specific `qc_pass_rate` values and locates its geometric
#' inflection estimate. The estimate can fall outside the literal schedule; it
#' does not create a corresponding metaclustering partition.
#'
#' @param collection.U Result from \code{\link{iteration.QC}}, a canonical score
#'   data frame with `k` and `qc_pass_rate`, or the deprecated score aliases.
#' @param basedata `"Curve"` or `"Points"`.
#' @param ggtitle Optional plot title.
#'
#' @return A list containing canonical `scores`, `fittedcurve`, `lfunction`,
#'   and `ggplot`. Deprecated `collection.U` is retained as a migration alias.
#' @seealso \code{\link{INFLECT}}, \code{\link{iteration.QC}},
#'   \code{\link{Lfunction}}
#' @export
QC.to.curve <- function(collection.U,
                        basedata,
                        ggtitle = NULL) {
  if (!basedata %in% c("Curve", "Points")) {
    stop('basedata should match either "Curve" or "Points"', call. = FALSE)
  }
  scores <- .inflect_score_frame(collection.U)
  if (nrow(scores) < 5L) {
    stop("QC scores must contain at least five rows.", call. = FALSE)
  }
  criterion <- if ("criterion" %in% names(scores)) {
    unique(scores$criterion)
  } else if (is.list(collection.U) &&
             !is.null(collection.U$provenance$criterion)) {
    collection.U$provenance$criterion
  } else {
    "combined"
  }
  if (length(criterion) != 1L) {
    stop("QC scores must describe exactly one aggregate criterion.", call. = FALSE)
  }

  fit <- drc::drm(
    data = scores,
    formula = qc_pass_rate ~ k,
    fct = drc::LL.4()
  )
  curve_k <- seq_len(max(scores$k))
  fittedcurve <- data.frame(
    k = curve_k,
    qc_pass_rate = stats::predict(
      object = fit,
      newdata = data.frame(k = curve_k)
    )
  )

  if (basedata == "Curve") {
    lfunction.data <- fittedcurve[-seq_len(4), , drop = FALSE]
    result.full <- Lfunction(lfunction.data, cutoff = 1000)
  } else {
    result.full <- Lfunction(scores[, c("k", "qc_pass_rate")], cutoff = 1000)
  }

  setcolors <- RColorBrewer::brewer.pal(name = "Set1", n = 9)
  if (basedata == "Curve") {
    part1 <- fittedcurve[
      fittedcurve$k >= min(lfunction.data$k) &
        fittedcurve$k <= result.full$knee,
      ,
      drop = FALSE
    ]
    part2 <- fittedcurve[
      fittedcurve$k >= result.full$knee &
        fittedcurve$k <= result.full$range,
      ,
      drop = FALSE
    ]
  } else {
    part1 <- scores[scores$k <= result.full$knee, , drop = FALSE]
    part2 <- scores[
      scores$k >= result.full$knee & scores$k <= result.full$range,
      ,
      drop = FALSE
    ]
  }
  touchline1 <- stats::lm(qc_pass_rate ~ k, part1)
  touchline2 <- stats::lm(qc_pass_rate ~ k, part2)
  label <- data.frame(
    k = part2$k[[1]],
    qc_pass_rate = part1$qc_pass_rate[[1]],
    value = part2$k[[1]]
  )
  criterion_label <- .inflect_criterion_label(criterion)
  figure <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = fittedcurve,
      ggplot2::aes(x = k, y = qc_pass_rate),
      color = setcolors[[1]]
    ) +
    ggplot2::ggtitle(ggtitle) +
    ggplot2::geom_line(
      data = data.frame(k = part1$k, qc_pass_rate = stats::fitted(touchline1)),
      ggplot2::aes(x = k, y = qc_pass_rate),
      color = "black"
    ) +
    ggplot2::geom_line(
      data = data.frame(k = part2$k, qc_pass_rate = stats::fitted(touchline2)),
      ggplot2::aes(x = k, y = qc_pass_rate),
      linetype = 2,
      color = "black"
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = result.full$knee),
      linetype = 3,
      color = setcolors[[9]]
    ) +
    ggplot2::geom_point(
      data = scores,
      ggplot2::aes(x = k, y = qc_pass_rate),
      color = setcolors[[2]],
      alpha = 0.7
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("# of SOM metaclusters") +
    ggplot2::ylab(paste0(criterion_label, " (%)")) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(scores$k), 25)
    ) +
    ggplot2::geom_label(
      data = label,
      ggplot2::aes(x = k, y = qc_pass_rate, label = value)
    )

  ## Compatibility columns live only in the deprecated aliases.
  legacy_scores <- data.frame(
    k = scores$k,
    qc_pass_rate = scores$qc_pass_rate,
    criterion = criterion,
    i = scores$k,
    Unimodality = scores$qc_pass_rate,
    stringsAsFactors = FALSE
  )
  fittedcurve$x <- fittedcurve$k
  fittedcurve$y <- fittedcurve$qc_pass_rate
  list(
    scores = scores,
    fittedcurve = fittedcurve,
    lfunction = result.full,
    ggplot = figure,
    collection.U = legacy_scores
  )
}
