#' @title Cluster quality control using diptest and IQR check
#'
#' @description Unimodality score is calculated for cluster results. Per marker
#' per cluster \link[diptest]{dip.test} is applied and inter-quartile range is
#' assessed. This function preserves the original INFLECT marker-level QC
#' criterion so users can inspect or reproduce the base statistic directly; the
#' faster package-level sweep reuses the same criterion through memoised helper
#' code.
#'
#' @param FlowSOM.results A supported SOM object with completed SOM clustering. Supports \pkg{FlowSOM} objects and \pkg{kohonen} objects returned by \code{\link[kohonen]{som}} or \code{\link[kohonen]{xyf}}.
#' @param metaclustering Vector with metacluster codes for all SOM-clusters.
#' @param zeroes.in Should be values at and below \code{0} be included. Recommended default for mass cytometry data is \code{FALSE}
#' @param only.clustering.markers If \code{TRUE} only evaluates markers specified as clustering markers. For \pkg{kohonen} objects this is the first data layer.
#' @param acquired_markers Vector of column names with marker data to be evaluated by fastINFLECT. Ignored if \code{only.clustering.markers == TRUE}
#' @param uniform.test What tests are performed per marker per cluster. Options are "both", "spread" , or "unimodality" as a string.
#' @param th.pvalue Threshold for rejecting Unimodality dip.test result. Default is \code{0.05}. For more information see \link[diptest]{dip.test}
#' @param th.IQR Threshold for rejecting marker distribution based on inter-quartile range. Default is arc-sinh transformed value of \code{2}.
#' @param verbose \code{logical} , default is \code{TRUE}
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
  }
  if (!inherits(FlowSOM.results, "inflect_som_view")) {
    FlowSOM.results <- as_inflect_som(FlowSOM.results)
  }
  if (is.null(metaclustering)) {
    stop("Error in FlowSOM.QC: The 'metaclustering' parameter can not be NULL")
  } else if (!is.integer(metaclustering)) {
    stop("Error in FlowSOM.QC: The 'metaclustering' parameter required a 'integer' of metaclustering results")
  }

  prep <- .inflect_prepare_qc(
    view = FlowSOM.results,
    only.clustering.markers = only.clustering.markers,
    acquired_markers = acquired_markers
  )
  p_of <- .inflect_make_p_of(list(...))
  accuracy.matrix <- .inflect_accuracy_matrix(
    prep = prep,
    metaclustering = metaclustering,
    zeroes.in = zeroes.in,
    uniform.test = uniform.test,
    th.pvalue = th.pvalue,
    th.IQR = th.IQR,
    p_of = p_of,
    cache = NULL,
    verbose = verbose
  )

  if (verbose) {
    message("[END] - generating Uniform Phenotypes QC")
  }
  invisible(accuracy.matrix)
}
