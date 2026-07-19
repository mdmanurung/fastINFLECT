#' @title Locate the knee of a diminishing-returns curve (Kneedle)
#'
#' @description Nonparametric knee/elbow detector following Satopää et al. (2011).
#' For a concave, increasing curve (as the fastINFLECT unimodality sweep typically is)
#' the knee is the x at which the normalised curve is farthest above the diagonal
#' joining its endpoints. Unlike the four-parameter log-logistic fit used by
#' \code{\link{QC.to.curve}}, it assumes no functional form, which makes it a robust
#' cross-check on dense sweeps.
#'
#' @param x Numeric vector of x coordinates (e.g. number of metaclusters).
#' @param y Numeric vector of y coordinates (e.g. unimodality score).
#' @param concave \code{logical}; \code{TRUE} (default) for a concave curve
#'   (diminishing returns). Set \code{FALSE} for a convex curve.
#' @param increasing \code{logical}; \code{TRUE} (default) if \code{y} rises with
#'   \code{x}.
#'
#' @return The \code{x} value at the detected knee, or \code{NA} if it cannot be
#'   determined (fewer than three finite points or a degenerate range).
#' @seealso \code{\link{inflect_threshold_k}}, \code{\link{QC.to.curve}}
#' @references Satopää, V., Albrecht, J., Irwin, D., & Raghavan, B. (2011).
#'   Finding a "kneedle" in a haystack. \emph{31st ICDCS Workshops}, 166-171.
#' @export
inflect_kneedle <- function(x, y, concave = TRUE, increasing = TRUE) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L) {
    return(NA_real_)
  }
  o <- order(x)
  x <- x[o]
  y <- y[o]

  rx <- diff(range(x))
  ry <- diff(range(y))
  if (rx == 0 || ry == 0) {
    return(NA_real_)
  }
  xn <- (x - min(x)) / rx
  yn <- (y - min(y)) / ry

  ## Signed distance from the chord joining the endpoints: yn = xn for an
  ## increasing curve, yn = 1 - xn for a decreasing one. A concave curve bulges
  ## above the chord (the knee is the maximum deviation); a convex curve bulges
  ## below it (the minimum). This is exact for all four shape combinations and
  ## reduces to the usual `max(yn - xn)` for the concave-increasing default.
  d <- if (increasing) yn - xn else yn - (1 - xn)
  idx <- if (concave) which.max(d) else which.min(d)
  x[idx]
}

#' @title Smallest k that reaches a target unimodality
#'
#' @description Returns the smallest number of metaclusters whose unimodality score
#' meets or exceeds \code{target}. This fastINFLECT addition operationalises the goal of admitting no
#' cluster with residual bimodal marker expression while avoiding over-clustering:
#' it is the first k at which (nearly) every (cluster, marker) pair is unimodal.
#'
#' @param collection.U Data frame with columns \code{i} and \code{Unimodality}
#'   (the percentage of unimodal (cluster, marker) pairs), as returned in
#'   \code{inflect.results$collection.U}.
#' @param target Desired unimodality. Values in \code{(0, 1]} are read as a
#'   fraction and values in \code{(1, 100]} as a percentage. Default \code{0.95}.
#'
#' @return The smallest \code{i} meeting the target, or \code{NA_integer_} if no
#'   tested k reaches it.
#' @seealso \code{\link{inflect_kneedle}}, \code{\link{INFLECT}}
#' @export
inflect_threshold_k <- function(collection.U, target = 0.95) {
  if (!is.data.frame(collection.U) ||
      !all(c("i", "Unimodality") %in% names(collection.U))) {
    stop("`collection.U` must be a data frame with `i` and `Unimodality` columns", call. = FALSE)
  }
  if (length(target) != 1L || !is.finite(target) || target <= 0) {
    stop("`target` must be a single positive number", call. = FALSE)
  }
  target_pct <- if (target <= 1) target * 100 else target
  ok <- collection.U$i[collection.U$Unimodality >= target_pct]
  if (length(ok) == 0L) {
    return(NA_integer_)
  }
  as.integer(min(ok))
}

## Assemble the recommended-k table reported alongside the LL.4 inflection point.
.inflect_selection <- function(collection.U, lfunction = NULL, target = 0.95) {
  df <- collection.U
  uni_at <- function(k) {
    if (is.na(k)) return(NA_real_)
    hit <- df$Unimodality[match(k, df$i)]
    if (length(hit) == 0L) NA_real_ else hit
  }
  knee <- if (!is.null(lfunction) && length(lfunction$knee) > 0) {
    lfunction$knee[[1]]
  } else {
    NA_real_
  }
  kneedle <- inflect_kneedle(df$i, df$Unimodality)
  thr <- inflect_threshold_k(df, target)

  target_pct <- if (target <= 1) target * 100 else target
  data.frame(
    method = c("inflection", "kneedle", "threshold"),
    k = c(knee, kneedle, thr),
    unimodality_at_k = c(uni_at(knee), uni_at(kneedle), uni_at(thr)),
    target = c(NA_real_, NA_real_, target_pct),
    stringsAsFactors = FALSE
  )
}
