#' @title Determine kneepoint by minimizing fitting errors
#'
#' @description For each candidate split of a point collection, fits two lines and computes the weighted combined RMSE. Returns the split index with the smallest combined error. This is the split criterion used by \code{\link{Lfunction}}.
#'
#' @param dataframe Data frame with cluster counts and corresponding QC pass
#'   rates, typically from \code{\link{iteration.QC}}.
#' @seealso \code{\link{INFLECT}} , \code{\link{QC.to.curve}},\code{\link{Lfunction}}
#'
#' @return Integer denoting the row of dataframe whose coordinate provides the kneepoint
#' @export
#'
leastError <- function(dataframe) {
  if (nrow(dataframe) < 5) {
    stop("`dataframe` must contain at least five rows", call. = FALSE)
  }

  eval <- list()
  last_split <- nrow(dataframe) - 1L
  first_split <- min(5L, last_split)

  for (row in first_split:last_split) {
    part1 <- dataframe[1:row, ]
    resid1 <- stats::resid(stats::lm(y ~ x, part1))
    rmse1 <- sqrt(mean(resid1 ^ 2))

    part2 <- dataframe[row:nrow(dataframe), ]
    resid2 <- stats::resid(stats::lm(y ~ x, part2))
    rmse2 <- sqrt(mean(resid2 ^ 2))

    combine <-
      (rmse1 * (row - 1) / (nrow(dataframe) - 1)) + rmse2 * (nrow(dataframe) -
                                                               row) / (nrow(dataframe) - 1)

    eval[[row]] <- combine
    names(eval[[row]]) <- paste0("xpoint", row)
  }
  eval[sapply(eval, is.null)] <- NA
  least <- which.min(unlist(eval))
  least

}
