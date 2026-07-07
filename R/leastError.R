#' @title Determine kneepoint based on minimizing fitting errors
#'
#' @description Evaluates the errors for two fitted lines on both parts of a collection of points. Function calculates the errors for a split on each element and selects the one with least amount of error. Internal function of \code{\link{Lfunction}}
#'
#' @param dataframe Dataframe with set.i and corresponding unimodality scores, resulting from \code{\link{iteration.QC}}
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

  for (row in 5:(nrow(dataframe))) {
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
