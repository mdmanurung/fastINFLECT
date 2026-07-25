#' @title Plot marker performance across metaclustering results
#'
#' @description Plots criterion-specific marker QC across metaclusterings. The
#' score is the percentage of clusters where a marker passed the aggregate
#' criterion selected in [INFLECT()]. It is a screening pass rate, not proof of
#' unimodality. The colour gradient denotes the tested cluster count.
#'
#' @param inflect.results A inflect.results object resulting from \code{\link{INFLECT}}.
#' @param ggtitle Optional. Character string to plot as title.
#' @param markers Which markers should be included in the plot? A vector of strings matching evaluated marker names. If \code{NULL}, all evaluated markers are displayed.
#'
#' @return A list with `marker.dataframe` and `plot`. The data frame contains
#' canonical `k`, `marker`, and `qc_pass_rate` columns. Deprecated `i`, `Marker`,
#' and `Performance` aliases are retained for compatibility.
#'
#' @examples
#' # Read in FlowSOM object from file. Downsampled clustering result of Levine32 dataset clustering.
#' # SOM-clustered to 375 clusters.
#' flowsom <- system.file("extdata", "Levine32sample.Rdata", package="fastINFLECT")
#' load(flowsom)
#' inflect.results <- INFLECT(
#'   FlowSOM.results = dataset,
#'   set.i = 5:12,
#'   multicore = FALSE
#' )
#'
#' # Display diagnostic graph
#' inflect.results$ggplot
#'
#' # Now check marker performance for all markers
#'marker.performance(inflect.results, ggtitle= "Levine32sample", markers=NULL)
#'
#'
#' @export

marker.performance <- function(inflect.results, ggtitle = NULL, markers = NULL) {
  if (!inherits(inflect.results, "inflect.results")) {
    stop("`inflect.results` must inherit from class 'inflect.results'", call. = FALSE)
  }

  accuracy.sets <- inflect.results$accuracy.sets
  if (is.null(accuracy.sets)) {
    accuracy.sets <- inflect.results$Accuracy.sets
  }
  if (is.null(accuracy.sets) || length(accuracy.sets) == 0) {
    stop("`inflect.results` must contain accuracy matrices", call. = FALSE)
  }

  if (!is.null(markers)){

    # Verify whether markers match colnames in accuracy matrices
    if (!is.character(markers)) {
      stop( "Error in marker.performance: markers should be a vector of strings" )
    } else  {
      if(!all(markers %in% colnames(accuracy.sets[[1]]))) {
        stop( "Error: markers should match colnames of matrices in inflect.results$Accuracy.sets" )
      }

    }

  } else {
    markers <- colnames(accuracy.sets[[1]])
  }

  success.rate <- lapply(accuracy.sets, function(x) {
    x <- x[, markers, drop = FALSE]
    colSums(x, na.rm = TRUE) * 100 / nrow(x)
  })
  success.rate <- as.data.frame(do.call(rbind, success.rate))
  success.rate$k <- as.numeric(rownames(success.rate))
  success.rate <- reshape2::melt(
    success.rate,
    variable.name = "marker",
    value.name = "qc_pass_rate",
    id.vars = "k"
  )
  success.rate$i <- success.rate$k
  success.rate$Marker <- success.rate$marker
  success.rate$Performance <- success.rate$qc_pass_rate


  # Dataframe is prepared, continue on to figure generation
  criterion_label <- if (!is.null(inflect.results$provenance$criterion)) {
    .inflect_criterion_label(inflect.results$provenance$criterion)
  } else {
    "criterion-specific QC pass rate"
  }
  figure <- ggplot2::ggplot(
    success.rate,
    ggplot2::aes(y = qc_pass_rate, x = marker)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0, color = "grey40") +
    ggplot2::geom_jitter(
      ggplot2::aes(color = k),
      shape = 1,
      alpha = 0.5
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_gradient(low = "#66C2A5", high = "#E41A1C") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::ylab(paste0(criterion_label, " within metaclusterings (%)"))
  if(!is.null(ggtitle)){figure <- figure + ggplot2::ggtitle(label = ggtitle)}
  return(list("marker.dataframe"= success.rate, "plot" = figure))
}
