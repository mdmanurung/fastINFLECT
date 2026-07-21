#' @title Plot marker performance across metaclustering results
#'
#' @description Plots the marker performance across metaclusterings in a boxplot. Performance score is the percentage of clusters
#' where a marker passed the \link[diptest]{dip.test} and IQR test per metaclustering. The color gradient denotes \code{set.i}.
#'
#' @param inflect.results A inflect.results object resulting from \code{\link{INFLECT}}.
#' @param ggtitle Optional. Character string to plot as title.
#' @param markers Which markers should be included in the plot? A vector of strings matching evaluated marker names. If \code{NULL}, all evaluated markers are displayed.
#'
#' @return \code{list} with 2 items. First is the melted dataframe with the marker performance percentage per marker per metaclustering. Second is a ggplot object, a boxplot with marker
#' performance on the y-axis, markers on the x-axis and the color scale denoting the amount of metaclusters evaluated.
#'
#' @examples
#' # Read in FlowSOM object from file. Downsampled clustering result of Levine32 dataset clustering.
#' # SOM-clustered to 375 clusters.
#' flowsom <- system.file("extdata", "Levine32sample.Rdata", package="fastINFLECT")
#' load(flowsom)
#' inflect.results<- INFLECT(FlowSOM.results= dataset, set.i= 5:12, multicore=FALSE, zeroes.in=FALSE)
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

  success.rate<- lapply(accuracy.sets, function(x){
    x<- x[, markers, drop = FALSE]
    apply(x, 2, sum, na.rm=TRUE)*100 / apply(x, 2, length)
  })
  success.rate<- as.data.frame(do.call(rbind, success.rate))
  success.rate$i <- as.numeric(rownames(success.rate))
  success.rate <- reshape2::melt(success.rate,variable.name="Marker", value.name= "Performance", id.vars="i")


  # Dataframe is prepared, continue on to figure generation
  figure <- ggplot2::ggplot(success.rate, ggplot2::aes(y = Performance, x = Marker, color = i)) +
    ggplot2::geom_boxplot(outlier.alpha = 0) +
    ggplot2::geom_jitter(shape = 1, alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_gradient(low = "#66C2A5", high = "#E41A1C") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::ylab("Marker performance (% passed within metaclusterings)")
  if(!is.null(ggtitle)){figure <- figure + ggplot2::ggtitle(label = ggtitle)}
  return(list("marker.dataframe"= success.rate, "plot" = figure))
}
