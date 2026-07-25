#' @description
#' fastINFLECT reimplements the INFLECT method for selecting a FlowSOM or
#' kohonen metaclustering endpoint from criterion-specific marker QC. The public
#' `INFLECT()` workflow is retained; the
#' implementation replaces repeated per-k QC work with memoised SOM-node
#' subtree scoring, faster dip-test p-value lookup, and Rcpp accelerators.
#'
#' The core INFLECT idea is to scan metaclusterings, calculate a marker-level QC
#' pass rate, and locate the inflection where additional clusters stop
#' improving that aggregate. Version 2.0 retains separate dip-test and
#' IQR-spread evidence, distinguishes fitted estimates from materialised
#' partitions, and avoids interpreting a pass rate as proof of unimodality.
#'
#' The original INFLECT implementation was developed by Jan Verhoeff in the
#' lab of JJ. Garcia-Vallejo and is available at
#' \url{https://github.com/jnverhoeff/GarciaVallejoLab}.
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
    "k",
    "Marker",
    "marker",
    "Performance",
    "Ui",
    "Unimodality",
    "qc_pass_rate",
    "value",
    "x",
    "y",
    "z"
  ))
}
