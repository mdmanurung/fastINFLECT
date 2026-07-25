#' @title Locate the knee of a diminishing-returns curve (Kneedle)
#'
#' @description Nonparametric knee/elbow detector following Satopää et al. (2011).
#' For a concave, increasing curve (as a fastINFLECT QC sweep typically is)
#' the knee is the x at which the normalised curve is farthest above the diagonal
#' joining its endpoints. Unlike the four-parameter log-logistic fit used by
#' \code{\link{QC.to.curve}}, it assumes no functional form, which makes it a robust
#' cross-check on dense sweeps.
#'
#' @param x Numeric vector of x coordinates (e.g. number of metaclusters).
#' @param y Numeric vector of y coordinates (for example, a QC pass rate).
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

#' @title Smallest k that reaches a target QC pass rate
#'
#' @description Returns the smallest tested number of metaclusters whose
#' criterion-specific QC pass rate meets or exceeds `target`. This aggregate is
#' a screening metric; it does not establish that clusters are truly unimodal.
#'
#' @param collection.U Canonical score data frame with `k` and `qc_pass_rate`,
#'   or the deprecated `i` and `Unimodality` aliases.
#' @param target Desired pass rate. Values in \code{(0, 1]} are read as a
#'   fraction and values in \code{(1, 100]} as a percentage. Default \code{0.95}.
#'
#' @return The smallest \code{i} meeting the target, or \code{NA_integer_} if no
#'   tested k reaches it.
#' @seealso \code{\link{inflect_kneedle}}, \code{\link{INFLECT}}
#' @export
inflect_threshold_k <- function(collection.U, target = 0.95) {
  scores <- .inflect_score_frame(collection.U)
  target_pct <- .inflect_target_percent(target)
  ok <- scores$k[scores$qc_pass_rate >= target_pct]
  if (length(ok) == 0L) {
    return(NA_integer_)
  }
  as.integer(min(ok))
}

.inflect_target_percent <- function(target) {
  if (length(target) != 1L ||
      !is.numeric(target) ||
      is.na(target) ||
      !is.finite(target) ||
      target <= 0) {
    stop("`target` must be a single finite number in (0, 100].", call. = FALSE)
  }
  target_pct <- if (target <= 1) target * 100 else target
  if (target_pct > 100) {
    stop("`target` must be a single finite number in (0, 100].", call. = FALSE)
  }
  as.numeric(target_pct)
}

.inflect_score_frame <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (is.data.frame(x$scores)) {
      x <- x$scores
    } else if (is.data.frame(x$U.set)) {
      x <- x$U.set
    }
  }
  if (!is.data.frame(x)) {
    stop("QC scores must be supplied as a data frame.", call. = FALSE)
  }
  if (all(c("k", "qc_pass_rate") %in% names(x))) {
    scores <- data.frame(
      k = x$k,
      qc_pass_rate = x$qc_pass_rate,
      stringsAsFactors = FALSE
    )
    if ("criterion" %in% names(x)) {
      scores$criterion <- x$criterion
    }
  } else if (all(c("i", "Unimodality") %in% names(x))) {
    scores <- data.frame(
      k = x$i,
      qc_pass_rate = x$Unimodality,
      stringsAsFactors = FALSE
    )
  } else {
    stop(
      "QC scores must contain `k` and `qc_pass_rate` columns.",
      call. = FALSE
    )
  }
  if (length(scores$k) == 0L ||
      anyNA(scores$k) ||
      anyNA(scores$qc_pass_rate) ||
      any(!is.finite(scores$k)) ||
      any(!is.finite(scores$qc_pass_rate)) ||
      any(scores$qc_pass_rate < 0 | scores$qc_pass_rate > 100)) {
    stop("QC scores must contain finite pass rates in [0, 100].", call. = FALSE)
  }
  scores
}

## Assemble the recommended-k table reported alongside the LL.4 inflection point.
.inflect_selection <- function(collection.U,
                               lfunction = NULL,
                               target = 0.95,
                               fittedcurve = NULL) {
  df <- .inflect_score_frame(collection.U)
  score_at <- function(k) {
    if (is.na(k)) {
      return(list(score = NA_real_, directly_tested = FALSE, source = NA_character_))
    }
    tested_idx <- match(k, df$k)
    if (!is.na(tested_idx)) {
      return(list(
        score = df$qc_pass_rate[[tested_idx]],
        directly_tested = TRUE,
        source = "tested"
      ))
    }

    fitted_score <- NA_real_
    if (is.data.frame(fittedcurve)) {
      if (all(c("k", "qc_pass_rate") %in% names(fittedcurve))) {
        fitted <- fittedcurve[, c("k", "qc_pass_rate"), drop = FALSE]
      } else if (all(c("x", "y") %in% names(fittedcurve))) {
        fitted <- fittedcurve[, c("x", "y"), drop = FALSE]
        names(fitted) <- c("k", "qc_pass_rate")
      } else {
        fitted <- NULL
      }
    } else {
      fitted <- NULL
    }
    if (!is.null(fitted)) {
      keep <- is.finite(fitted$k) & is.finite(fitted$qc_pass_rate)
      fitted <- fitted[keep, , drop = FALSE]
      exact_idx <- match(k, fitted$k)
      if (!is.na(exact_idx)) {
        fitted_score <- fitted$qc_pass_rate[[exact_idx]]
      } else if (nrow(fitted) >= 2L) {
        fitted_score <- stats::approx(
          x = fitted$k,
          y = fitted$qc_pass_rate,
          xout = k,
          rule = 2,
          ties = mean
        )[["y"]]
      }
    }
    list(
      score = fitted_score,
      directly_tested = FALSE,
      source = if (is.na(fitted_score)) NA_character_ else "fitted"
    )
  }
  knee <- if (!is.null(lfunction) && length(lfunction$knee) > 0) {
    lfunction$knee[[1]]
  } else {
    NA_real_
  }
  kneedle <- inflect_kneedle(df$k, df$qc_pass_rate)
  thr <- inflect_threshold_k(df, target)

  target_pct <- .inflect_target_percent(target)
  ks <- c(knee, kneedle, thr)
  score_info <- lapply(ks, score_at)
  directly_tested <- vapply(score_info, `[[`, logical(1), "directly_tested")
  score_source <- vapply(score_info, `[[`, character(1), "source")
  qc_pass_rate <- vapply(score_info, `[[`, numeric(1), "score")
  data.frame(
    method = c("inflection", "kneedle", "threshold"),
    k = ks,
    qc_pass_rate_at_k = qc_pass_rate,
    directly_tested = directly_tested,
    partition_available = directly_tested,
    k_status = ifelse(
      is.na(ks),
      "not_available",
      ifelse(directly_tested, "tested_partition", "fitted_estimate_no_partition")
    ),
    score_source = score_source,
    target = c(NA_real_, NA_real_, target_pct),
    unimodality_at_k = qc_pass_rate,
    stringsAsFactors = FALSE
  )
}
