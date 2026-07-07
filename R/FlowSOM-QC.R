#' @title Cluster quality control using diptest and IQR check
#'
#' @description Unimodality score is calculated for cluster results. Per marker per cluster \link[diptest]{dip.test} is applied and inter-quartile range is assessed.
#'
#' @param FlowSOM.results A FlowSOM object with completed SOM clustering, either after full FlowSOM function or after BuildSOM
#' @param metaclustering Vector with metacluster codes for all SOM-clusters.
#' @param zeroes.in Should be values at and below \code{0} be included. Recommended default for mass cytometry data is \code{FALSE}
#' @param only.clustering.markers If \code{TRUE} only evaluates markers specified in FlowSOM.results$map$colsUsed
#' @param acquired_markers Vector of column names with marker data to be evaluated by INFLECT. Ignored if \code{only.clustering.markers == TRUE}
#' @param uniform.test What tests are performed per marker per cluster. Options are "both", "spread" , or "unimodality" as a string.
#' @param th.pvalue Threshold for rejecting Unimodality dip.test result. Default is \code{0.05}. For more information see \link[diptest]{dip.test}
#' @param th.IQR Threshold for rejecting marker distribution based on inter-quartile range. Default is arc-sinh transformed value of \code{2}.
#' @param verbose \code{logical} , default is \code{FALSE}
#' @param ... Additional arguments to pass to \code{\link[diptest]{dip.test}}.
#'
#' @return A \code{matrix} with evaluated markers in columns and clusters in rows. Each position in the matrix is \code{logical} indicating a pass or a fail.
#' @seealso \code{\link{INFLECT}} , \code{\link{iteration.QC}}
#'
#' @export
FlowSOMQC <- function(FlowSOM.results,
                      metaclustering,
                      zeroes.in = FALSE,
                      only.clustering.markers = TRUE,
                      acquired_markers = NULL,
                      uniform.test = c("both", "spread", "unimodality"),
                      th.pvalue = 0.05,
                      th.IQR = 2,
                      verbose = TRUE,
                      ...)
{
  uniform.test <- match.arg(uniform.test)

  if (is.null(FlowSOM.results)) {
    stop("Error in FlowSOM.QC: The 'FlowSOM.results' parameter can not be NULL")
  } else if (!inherits(FlowSOM.results, "FlowSOM")) {
    stop("Error in FlowSOM.QC: The 'FlowSOM.results' parameter required a 'FlowSOM' object")
  }
  if (is.null(metaclustering)) {
    stop("Error in FlowSOM.QC: The 'metaclustering' parameter can not be NULL")
  } else if (!is.integer(metaclustering)) {
    stop("Error in FlowSOM.QC: The 'metaclustering' parameter required a 'integer' of metaclustering results")
  }

  data <- FlowSOM.results$data
  if (isTRUE(FlowSOM.results$scale)) {
    for (j in seq_len(ncol(data))) {
      data[, j] <- data[, j] * FlowSOM.results$scaled.scale[j] + FlowSOM.results$scaled.center[j]
    }
  }
  colnames(data) <- FlowSOM.results$prettyColnames
  data <- cbind(data, "cluster" = metaclustering[FlowSOM.results$map$mapping[, 1]])
  clusters <- seq_len(max(metaclustering))

  clustering.markers <- FlowSOM.results$prettyColnames[FlowSOM.results$map$colsUsed]
  if (only.clustering.markers) {
    markers <- clustering.markers
  } else {
    if (!is.null(acquired_markers) && all(acquired_markers %in% FlowSOM.results$prettyColnames)) {
      markers <- acquired_markers
    } else {
      stop("Error in acquired_markers: The 'acquired_markers' vector must match names in 'FlowSOM.result$prettyColnames' ")
    }
  }

  ordered.markers <- c(
    gtools::mixedsort(intersect(markers, clustering.markers)),
    gtools::mixedsort(setdiff(markers, clustering.markers))
  )
  accuracy.matrix <- matrix(nrow = length(clusters), ncol = length(markers),
                            dimnames = list(clusters, ordered.markers))
  cluster.rows <- split(seq_len(nrow(data)), data[, "cluster"])

  count <- 0
  for (cluster in clusters) {
    if (verbose) {
      count <- count + 1
      message(paste0("Cluster: ", count, " on ", length(clusters)))
    }

    rows <- cluster.rows[[as.character(cluster)]]
    if (is.null(rows)) {
      expressions <- data[integer(0), ordered.markers, drop = FALSE]
    } else {
      expressions <- data[rows, ordered.markers, drop = FALSE]
    }

    for (marker in ordered.markers) {
      if (nrow(expressions) > 1) {
        values <- expressions[, marker]
        if (zeroes.in == FALSE) {
          marker.expression <- values[values > 0]
          if (length(marker.expression) < 5) {
            marker.expression <- c(rep(0, 5 - length(marker.expression)), marker.expression)
          }
        } else {
          marker.expression <- values
        }

        uniform <- TRUE
        if (uniform.test == "unimodality" || uniform.test == "both") {
          p.value <- diptest::dip.test(marker.expression,
                                       ...)$p.value
          uniform <- uniform && p.value >= th.pvalue
        }
        if (uniform.test == "spread" || uniform.test == "both") {
          marker.quantile <- stats::quantile(marker.expression)
          marker.iqr <- marker.quantile[4] - marker.quantile[2]
          uniform <- uniform && marker.iqr < th.IQR

        }

        accuracy.matrix[cluster, marker] <- uniform
      }
      else {
        accuracy.matrix[cluster, marker] <- NA
      }
    }

  }

  if (verbose) {
    message("[END] - generating Uniform Phenotypes QC")
  }
  invisible(accuracy.matrix)
}
