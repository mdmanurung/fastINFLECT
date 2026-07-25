#' @title Run iterative metaclustering on SOM-clustered objects
#'
#' @description Computes metaclusterings for all cluster numbers in `set.i`.
#' The Ward.D2 hierarchy is built once from the SOM codebook and then cut at
#' each requested k, rather than recomputing a separate clustering per k.
#'
#' @param FlowSOM.results A supported SOM object with completed SOM clustering. Supports \pkg{FlowSOM} objects and \pkg{kohonen} objects returned by \code{\link[kohonen]{som}} or \code{\link[kohonen]{xyf}}.
#' @param set.i Literal, unique, strictly increasing integer cluster counts
#'   within the SOM-node range.
#' @param multicore Retained for backward compatibility. Hierarchical clustering is now computed once and is not parallelized.
#' @param cores Retained for backward compatibility.
#'
#' @return metaclustering.list A \code{list} of \code{arrays} with metacluster-codes for the SOM nodes within the input object.
#' @seealso \code{\link{INFLECT}}
#'
#' @export
iteration.metacluster <- function(FlowSOM.results, set.i, multicore = FALSE, cores = NULL) {
  FlowSOM.results <- as_inflect_som(FlowSOM.results)
  set.i <- .inflect_validate_qc_schedule(
    set.i,
    as.integer(FlowSOM.results$map$nNodes)
  )
  invisible(.inflect_resolve_cores(multicore, cores))
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
