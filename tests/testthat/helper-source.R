find_package_root <- function(start = getwd()) {
  start <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(start, "DESCRIPTION")) &&
        dir.exists(file.path(start, "R"))) {
      return(start)
    }
    parent <- dirname(start)
    if (identical(parent, start)) {
      stop("Could not find package root")
    }
    start <- parent
  }
}

pkg_root <- find_package_root()

source_pkg_file <- function(file, envir = globalenv()) {
  sys.source(file.path(pkg_root, "R", file), envir = envir)
}

make_fake_flowsom <- function(values, cols_used = seq_len(ncol(values))) {
  values <- as.matrix(values)
  mapping <- matrix(rep(1L, nrow(values)), ncol = 1)
  object <- list(
    data = values,
    scale = FALSE,
    prettyColnames = colnames(values),
    map = list(
      colsUsed = cols_used,
      mapping = mapping,
      nNodes = 1L,
      codes = values[1, , drop = FALSE]
    )
  )
  class(object) <- "FlowSOM"
  object
}
