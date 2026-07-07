#' @title Run iterative metaclustering on SOM-clustered FlowSOM object
#'
#' @description Metaclustering runs are determined for all cluster numbers in set.i . Hierarchical clustering is computed once and cut for all requested cluster numbers.
#'
#' @param FlowSOM.results A FlowSOM object with completed SOM clustering, either after full FlowSOM function or after BuildSOM
#' @param set.i Vector containing either the desired iterations to be tested
#' @param multicore Retained for backward compatibility. Hierarchical clustering is now computed once and is not parallelized.
#' @param cores Retained for backward compatibility.
#'
#' @return metaclustering.list A \code{list} of \code{arrays} with metacluster-codes for the SOM-clusters within the FlowSOM.results object.
#' @seealso \code{\link{INFLECT}}
#'
#' @export
iteration.metacluster <- function(FlowSOM.results, set.i, multicore = TRUE, cores = NULL) {
  codes <- FlowSOM.results$map$codes
  fit <- stats::hclust(stats::dist(codes, method = "minkowski"), method = "ward.D2")
  cutree.result <- stats::cutree(fit, k = set.i)

  if (is.null(dim(cutree.result))) {
    metaclustering.list <- list(cutree.result)
  } else {
    metaclustering.list <- lapply(seq_along(set.i), function(i) cutree.result[, i])
  }

  stats::setNames(metaclustering.list, as.character(set.i))
}
