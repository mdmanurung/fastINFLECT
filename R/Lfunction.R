#' @title Lfunction for determining Inflection Point
#'
#' @description Wrapper around \code{\link{leastError}} with an option to refine the Inflection Point by limiting the range of points evaluated for fitting line errors. Refinement is disabled by default.
#'
#' @param totaldataframe Data frame with cluster counts and corresponding QC
#'   pass rates, either observed points or a fitted curve.
#' @param cutoff Integer. Initial Inflection Point calculated on the entire curve is multiplied by \code{cutoff} to determine the new range of the curve that is used as input for \code{\link{leastError}}
#' @param plot Logical. If \code{TRUE}, draw the legacy diagnostic base plot.
#' @seealso \code{\link{INFLECT}} , \code{\link{QC.to.curve}},\code{\link{leastError}}
#'
#' @return \code{list} with 3 items: Inflection Point, the endpoint of the second touchline, and the angle between the two touchlines
#'
#' @export

Lfunction <- function(totaldataframe, cutoff = 1000, plot = FALSE) {
  if (ncol(totaldataframe) < 2) {
    stop("`totaldataframe` must contain at least two columns", call. = FALSE)
  }
  totaldataframe <- data.frame(
    x = totaldataframe[[1]],
    y = totaldataframe[[2]]
  )

  kneepoint <- leastError(totaldataframe)

  part1 <- totaldataframe[1:kneepoint, ]
  part2 <- totaldataframe[kneepoint:nrow(totaldataframe), ]
  test1 <- stats::lm(y ~ x, part1)
  test2 <- stats::lm(y ~ x, part2)
  angle1 <-
    LinesAngles(rev(test1$coefficients), rev(test2$coefficients))

  newrange <- max(5, as.integer(kneepoint * cutoff))
  oldrange <- newrange
  if (isTRUE(plot)) {
    angleplot(part1, part2, test1, test2, totaldataframe[kneepoint, 1], angle1)
  }
  if (newrange > nrow(totaldataframe)) {
    return(list(
      "knee" = totaldataframe[kneepoint, 1],
      "range" = totaldataframe[nrow(totaldataframe), 1],
      "angle" = angle1
    ))
  } else {
    repeat {
      oldangle <- angle1
      oldkneepoint <- kneepoint
      kneepoint <- leastError(totaldataframe[1:newrange, ])
      part1 <- totaldataframe[1:kneepoint, ]
      part2 <- totaldataframe[kneepoint:newrange, ]
      test1 <- stats::lm(y ~ x, part1)
      test2 <- stats::lm(y ~ x, part2)

      angle1 <- LinesAngles(rev(test1$coefficients), rev(test2$coefficients))
      if (isTRUE(plot)) {
        angleplot(part1, part2, test1, test2, main = totaldataframe[kneepoint, 1], angle1)
      }
      if (kneepoint >= oldkneepoint) {
        oldrange <- newrange
        break
      }
      oldrange <- newrange

      newrange <- min(nrow(totaldataframe), max(20, as.integer(kneepoint * cutoff)))

    }
    return(list(
      "knee" = totaldataframe[oldkneepoint, 1],
      "range" = totaldataframe[oldrange, 1],
      "angle" = oldangle
    ))
  }

}
