#' @title Hierarchical (Ward) metaclustering of SOM codes
#'
#' @description Performs hierarchical metaclustering of SOM cluster codes using Ward's D2 linkage, adapted from FlowSOM.
#'
#' @param data Matrix with the median values for each clustering parameter for all SOM clusters.
#' @param nClus Amount of metaclusters to be obtained.
#'
#' @return  A \code{scalar} of \code{length(SOMclusters)} with metacluster identity.
#' @seealso \code{\link{INFLECT}}
#'
#' @export

metaClusteringhclust <- function(data, nClus){
  d <- stats::dist(data, method = "minkowski")
  fit <- stats::hclust(d, method="ward.D2")
  stats::cutree(fit, k=nClus)
}
