#' @description
#' fastINFLECT is a fast reimplementation of the original INFLECT package for
#' selecting a FlowSOM or kohonen metaclustering endpoint from marker
#' unimodality. The public `INFLECT()` workflow and the original diagnostic
#' interpretation are retained; the implementation replaces repeated per-k QC
#' work with memoised SOM-node subtree scoring, faster dip-test p-value lookup,
#' and small Rcpp accelerators.
#'
#' The original INFLECT method established the package's core idea: choose k by
#' scanning metaclusterings, measuring whether marker expression is unimodal
#' within each cluster, and locating the inflection point where additional
#' clusters stop improving that score. fastINFLECT acknowledges that work and
#' focuses this repository on making the same method practical for dense sweeps
#' and benchmarked comparisons.
#'
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
#' @importFrom diptest dip dip.test
#' @importFrom drc drm LL.4
#' @importFrom gtools mixedsort
#' @importFrom parallel detectCores mclapply
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Rcpp evalCpp
#' @importFrom reshape2 melt
#' @importFrom stats approx cutree density dist fitted hclust lm predict quantile resid setNames
#' @importFrom utils data
#' @useDynLib fastINFLECT, .registration = TRUE
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
