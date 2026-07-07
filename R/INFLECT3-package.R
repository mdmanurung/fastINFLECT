#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import FlowSOM
#' @importFrom graphics abline
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics text
#' @importFrom LearnGeom LinesAngles
#' @importFrom diptest dip.test
#' @importFrom doParallel registerDoParallel
#' @importFrom drc drm LL.4
#' @importFrom foreach foreach %dopar%
#' @importFrom gtools mixedsort
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom stats cutree density dist fitted hclust lm predict quantile resid setNames
## usethis namespace: end
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "i",
    "Marker",
    "Performance",
    "Ui",
    "Unimodality",
    "x",
    "y",
    "z"
  ))
}
