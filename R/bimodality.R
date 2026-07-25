#' @title Sarle's bimodality coefficient
#'
#' @description Computes Sarle's bimodality coefficient (BC) for a numeric sample,
#' a fast moment-based screen for bimodality:
#' \deqn{BC = (g^2 + 1) / (k + 3 (n-1)^2 / ((n-2)(n-3)))}
#' where \eqn{g} is the sample skewness and \eqn{k} the sample excess kurtosis.
#' BC lies in \code{(0, 1]}; the benchmark \code{5/9} (about 0.555) is the value for a
#' uniform distribution, and larger values indicate a more bimodal shape.
#'
#' BC is a moment-based diagnostic that can complement the dip-test evidence
#' retained by fastINFLECT. It is deliberately \emph{not} the
#' default QC statistic: because it is driven by skewness and kurtosis it can
#' flag heavy-tailed or strongly skewed single-mode markers (common in
#' cytometry) as bimodal. Use it only to rank marker-cluster pairs for
#' follow-up with explicit distributional tests and plots.
#'
#' @param x A numeric vector (missing values are dropped). At least four finite
#'   values are required.
#' @param na.rm \code{logical}; drop missing values before computing. Default
#'   \code{TRUE}.
#'
#' @return A single numeric bimodality coefficient, or \code{NA_real_} when there
#'   are too few points or the sample has zero variance.
#' @seealso \code{\link{FlowSOMQC}}, \code{\link{computemode}}
#' @examples
#' bimodality.coefficient(rnorm(500))                       # ~0.35, unimodal
#' bimodality.coefficient(c(rnorm(250), rnorm(250, 6)))     # high, bimodal
#' @export
bimodality.coefficient <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (n < 4L) {
    return(NA_real_)
  }
  m <- mean(x)
  d <- x - m
  m2 <- mean(d^2)
  if (m2 <= 0) {
    return(NA_real_)
  }
  m3 <- mean(d^3)
  m4 <- mean(d^4)

  ## sample-corrected skewness (G1) and excess kurtosis (G2), matching the SAS /
  ## mousetrap definitions of the coefficient.
  g1 <- m3 / m2^1.5
  G1 <- g1 * sqrt(n * (n - 1)) / (n - 2)
  g2 <- m4 / m2^2 - 3
  G2 <- ((n + 1) * g2 + 6) * (n - 1) / ((n - 2) * (n - 3))

  (G1^2 + 1) / (G2 + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}
