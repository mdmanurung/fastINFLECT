#' @title Plot diagnostic fastINFLECT curve and find inflection point
#'
#' @description Plots the `collection.U` score against `set.i` and uses
#' \code{\link{Lfunction}} to locate the inflection point for the optimal k.
#' The score curve may come from \code{\link{iteration.QC}} or any compatible
#' source; the curve-fitting procedure follows the original INFLECT method.
#'
#' @param collection.U List returned by \code{\link{iteration.QC}}. The first element must be a dataframe with set.i and corresponding unimodality scores.
#' @param basedata Data to be used to calculate inflection point, given as a string. Options are \code{Curve} and \code{Points}
#' @param ggtitle Optional title for resulting diagnostic graph. Default \code{NULL}
#' @return A \code{list} with 4 items. First is a data.frame with the unimodality scores for each metaclustering. Second is a dataframe with the fitted curve. Third is the result of \code{\link{Lfunction}}, lastly the diagnostic plot created from the other three items.
#'
#' @seealso \code{\link{INFLECT}} , \code{\link{iteration.QC}},\code{\link{Lfunction}}
#'
#' @export
QC.to.curve <-
  function(collection.U,
           basedata,
           ggtitle = NULL) {
    if (!basedata %in% c("Curve", "Points")) {
      stop('basedata should match either "Curve" or "Points"', call. = FALSE)
    }
    df.points <- collection.U[[1]]
    if (!is.data.frame(df.points) ||
        !all(c("i", "Unimodality") %in% colnames(df.points))) {
      stop("`collection.U` must contain a data frame with `i` and `Unimodality` columns", call. = FALSE)
    }
    if (nrow(df.points) < 5) {
      stop("`collection.U` must contain at least five rows", call. = FALSE)
    }

    drc <- drc::drm(data = df.points,
                    formula = Unimodality ~ i,
                    fct = drc::LL.4())
    curve.i <- seq_len(max(df.points$i))
    df.curve <-
      data.frame(x = curve.i,
                 y = stats::predict(object = drc, newdata = data.frame(i = curve.i)))

    if (basedata == "Curve") {
      lfunction.data <- df.curve[-seq_len(4), ]
      result.full <- Lfunction(lfunction.data, cutoff = 1000)

    }

    if (basedata == "Points") {
      result.full <- Lfunction(df.points, cutoff = 1000)
    }

    setcolors <- RColorBrewer::brewer.pal(name = "Set1", n = 9)

    if (basedata == "Curve") {
      part1 <- df.curve[df.curve$x >= min(lfunction.data$x) &
                          df.curve$x <= result.full$knee, , drop = FALSE]
      part2 <- df.curve[df.curve$x >= result.full$knee &
                          df.curve$x <= result.full$range, , drop = FALSE]
      touchline1 <- stats::lm(y ~ x, part1)
      touchline2 <- stats::lm(y ~ x, part2)
      label <-
        data.frame(x = part2[1, 1], y = part1[1, 2], z = part2[1, 1])
      title <- ggtitle

      figure <-
        ggplot2::ggplot() +
        ggplot2::geom_line(data = df.curve, ggplot2::aes(x = x, y = y), color = setcolors[1]) +
        ggplot2::ggtitle(title) +
        ggplot2::geom_line(data = data.frame(x = part1[, 1], y = stats::fitted(touchline1)), ggplot2::aes(x = x, y = y), color = "black") +
        ggplot2::geom_line(data = data.frame(x = part2[, 1], y = stats::fitted(touchline2)), ggplot2::aes(x = x, y = y), linetype = 2, color = "black") +
        ggplot2::geom_vline(ggplot2::aes(xintercept = result.full$knee), linetype = 3, color = setcolors[9]) +
        ggplot2::geom_point(data = df.points, ggplot2::aes(x = i, y = Unimodality), color = setcolors[2], alpha = 0.7) +
        ggplot2::theme_classic() +
        ggplot2::xlab("# of SOM metaclusters") +
        ggplot2::ylab("% unimodal distributions across clusters") +
        ggplot2::scale_x_continuous(breaks = seq(0, max(df.points$i), 25)) +
        ggplot2::geom_label(data = label, ggplot2::aes(x = x, y = y, label = z))

    }
    if (basedata == "Points") {
      part1 <- df.points[df.points$i <= result.full$knee, , drop = FALSE]
      part2 <- df.points[df.points$i >= result.full$knee &
                           df.points$i <= result.full$range, , drop = FALSE]
      touchline1 <- stats::lm(Unimodality ~ i, part1)
      touchline2 <- stats::lm(Unimodality ~ i, part2)
      label <-
        data.frame(x = part2[1, 1], y = part1[1, 2], z = part2[1, 1])
      title <- ggtitle

      figure <-
        ggplot2::ggplot() +
        ggplot2::geom_line(data = df.curve, ggplot2::aes(x = x, y = y), color = setcolors[1]) +
        ggplot2::ggtitle(title) +
        ggplot2::geom_line(data = data.frame(x = part1[, 1], y = stats::fitted(touchline1)), ggplot2::aes(x = x, y = y), color = "black") +
        ggplot2::geom_line(data = data.frame(x = part2[, 1], y = stats::fitted(touchline2)), ggplot2::aes(x = x, y = y), linetype = 2, color = "black") +
        ggplot2::geom_vline(ggplot2::aes(xintercept = result.full$knee), linetype = 3, color = setcolors[9]) +
        ggplot2::geom_point(data = df.points, ggplot2::aes(x = i, y = Unimodality), color = setcolors[2], alpha = 0.7) +
        ggplot2::theme_classic() +
        ggplot2::xlab("# of SOM metaclusters") +
        ggplot2::ylab("% unimodal distributions across clusters") +
        ggplot2::scale_x_continuous(breaks = seq(0, max(df.points$i), 25)) +
        ggplot2::geom_label(data = label, ggplot2::aes(x = x, y = y, label = z))
    }

    return(
      list(
        collection.U = df.points,
        fittedcurve = df.curve ,
        lfunction = result.full,
        ggplot = figure
      )
    )

  }
