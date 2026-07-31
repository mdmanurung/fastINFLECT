#' @description
#' fastINFLECT evaluates a literal schedule of metacluster counts for a FlowSOM
#' or kohonen self-organising map. Start with [INFLECT()], then use `print()`,
#' `plot()`, and `as.data.frame()` to read the candidate values and tested QC
#' pass rates.
#'
#' The result retains the tested partitions, separate dip and IQR decisions,
#' detailed cluster-marker measurements, and run provenance. Candidate values
#' narrow the search for `k`; they do not establish biological validity or
#' prove unimodality.
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
