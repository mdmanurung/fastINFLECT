as_inflect_som <- function(object) {
  if (inherits(object, "inflect_som_view")) {
    return(object)
  }
  if (inherits(object, "FlowSOM")) {
    return(as_inflect_flowsom(object))
  }
  if (inherits(object, "kohonen")) {
    return(as_inflect_kohonen(object))
  }

  stop(
    "`FlowSOM.results` must inherit from 'FlowSOM' or 'kohonen'",
    call. = FALSE
  )
}

as_inflect_flowsom <- function(object) {
  if (is.null(object$map) || is.null(object$map$codes)) {
    stop("FlowSOM object must contain `map$codes`", call. = FALSE)
  }
  if (is.null(object$data)) {
    stop("FlowSOM object must contain `data`", call. = FALSE)
  }
  if (is.null(object$map$mapping)) {
    stop("FlowSOM object must contain `map$mapping`", call. = FALSE)
  }

  data <- inflect_matrix(object$data, "`FlowSOM.results$data`")
  if (!is.numeric(data)) {
    stop("FlowSOM `data` must be numeric", call. = FALSE)
  }
  storage.mode(data) <- "double"

  pretty_colnames <- object$prettyColnames
  if (is.null(pretty_colnames)) {
    pretty_colnames <- colnames(data)
  }
  if (is.null(pretty_colnames)) {
    pretty_colnames <- paste0("marker", seq_len(ncol(data)))
  }
  if (length(pretty_colnames) != ncol(data)) {
    stop("FlowSOM marker names must match the number of data columns", call. = FALSE)
  }

  codes <- inflect_matrix(object$map$codes, "`FlowSOM.results$map$codes`")
  if (!is.numeric(codes)) {
    stop("FlowSOM `map$codes` must be numeric", call. = FALSE)
  }
  storage.mode(codes) <- "double"

  n_nodes <- object$map$nNodes
  if (is.null(n_nodes)) {
    n_nodes <- nrow(codes)
  }
  n_nodes <- as.integer(n_nodes)
  if (length(n_nodes) != 1L || is.na(n_nodes) || n_nodes < 1L) {
    stop("FlowSOM `map$nNodes` must be a positive integer", call. = FALSE)
  }
  if (n_nodes != nrow(codes)) {
    stop("FlowSOM `map$nNodes` must match the number of code rows", call. = FALSE)
  }

  mapping <- object$map$mapping
  if (is.null(dim(mapping))) {
    mapping <- matrix(mapping, ncol = 1)
  }
  mapping <- as.integer(mapping[, 1])
  if (length(mapping) != nrow(data)) {
    stop("FlowSOM `map$mapping` length must match the number of data rows", call. = FALSE)
  }
  if (anyNA(mapping)) {
    stop("FlowSOM `map$mapping` must not contain missing values", call. = FALSE)
  }
  if (any(mapping < 1L) || any(mapping > n_nodes)) {
    stop("FlowSOM `map$mapping` contains node ids outside `seq_len(map$nNodes)`", call. = FALSE)
  }

  cols_used <- object$map$colsUsed
  if (is.null(cols_used)) {
    cols_used <- seq_len(ncol(data))
  }
  cols_used <- as.integer(cols_used)
  if (length(cols_used) == 0L || anyNA(cols_used) ||
      any(cols_used < 1L) || any(cols_used > ncol(data))) {
    stop("FlowSOM `map$colsUsed` must contain valid data column indices", call. = FALSE)
  }

  new_inflect_som_view(
    data = data,
    pretty_colnames = pretty_colnames,
    cols_used = cols_used,
    mapping = mapping,
    n_nodes = n_nodes,
    codes = codes,
    scale = isTRUE(object$scale),
    scaled_scale = object$scaled.scale,
    scaled_center = object$scaled.center,
    source = list(
      type = "FlowSOM",
      class = class(object),
      data_layer = NA_character_,
      code_layers = NA_character_
    )
  )
}

as_inflect_kohonen <- function(object) {
  data_layers <- inflect_as_layer_list(object$data, "`FlowSOM.results$data`")
  code_layers <- inflect_as_layer_list(object$codes, "`FlowSOM.results$codes`")
  data_layer_names <- inflect_layer_names(data_layers)
  code_layer_names <- inflect_layer_names(code_layers)

  data_index <- 1L
  data <- inflect_matrix(
    data_layers[[data_index]],
    paste0("kohonen data layer `", data_layer_names[[data_index]], "`")
  )
  if (!is.numeric(data)) {
    stop("Kohonen data layer used for QC must be numeric", call. = FALSE)
  }
  storage.mode(data) <- "double"

  pretty_colnames <- colnames(data)
  if (is.null(pretty_colnames)) {
    pretty_colnames <- paste0("marker", seq_len(ncol(data)))
  }

  code_indices <- inflect_kohonen_code_indices(object, code_layers)
  codes <- inflect_bind_kohonen_codes(
    code_layers = code_layers,
    code_indices = code_indices,
    layer_names = code_layer_names,
    distance_weights = object$distance.weights,
    user_weights = object$user.weights
  )

  mapping <- object$unit.classif
  if (is.null(mapping)) {
    stop("Kohonen object must contain `unit.classif`", call. = FALSE)
  }
  if (length(mapping) != nrow(data)) {
    stop(
      "Kohonen `unit.classif` length must match the selected data layer rows",
      call. = FALSE
    )
  }
  mapping <- as.integer(mapping)
  if (anyNA(mapping)) {
    stop("Kohonen `unit.classif` must not contain missing values", call. = FALSE)
  }
  if (any(mapping < 1) || any(mapping > nrow(codes))) {
    stop("Kohonen `unit.classif` contains unit ids outside the code matrix", call. = FALSE)
  }

  new_inflect_som_view(
    data = data,
    pretty_colnames = pretty_colnames,
    cols_used = seq_len(ncol(data)),
    mapping = mapping,
    n_nodes = nrow(codes),
    codes = codes,
    scale = FALSE,
    scaled_scale = NULL,
    scaled_center = NULL,
    source = list(
      type = "kohonen",
      class = class(object),
      data_layer = data_layer_names[[data_index]],
      code_layers = code_layer_names[code_indices]
    )
  )
}

new_inflect_som_view <- function(data,
                                 pretty_colnames,
                                 cols_used,
                                 mapping,
                                 n_nodes,
                                 codes,
                                 scale,
                                 scaled_scale,
                                 scaled_center,
                                 source) {
  structure(
    list(
      data = data,
      scale = scale,
      scaled.scale = scaled_scale,
      scaled.center = scaled_center,
      prettyColnames = pretty_colnames,
      map = list(
        colsUsed = cols_used,
        mapping = matrix(as.integer(mapping), ncol = 1),
        nNodes = as.integer(n_nodes),
        codes = codes
      ),
      inflect_source = source
    ),
    class = c("inflect_som_view", "list")
  )
}

inflect_matrix <- function(x, label) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop(label, " must be a matrix or data frame", call. = FALSE)
  }
  x
}

inflect_as_layer_list <- function(x, label) {
  if (is.null(x)) {
    stop(label, " must not be NULL", call. = FALSE)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    return(list(X = x))
  }
  if (!is.list(x) || length(x) == 0) {
    stop(label, " must be a matrix, data frame, or non-empty list", call. = FALSE)
  }
  x
}

inflect_layer_names <- function(layers) {
  layer_names <- names(layers)
  if (is.null(layer_names)) {
    layer_names <- rep("", length(layers))
  }
  missing <- is.na(layer_names) | layer_names == ""
  layer_names[missing] <- paste0("layer", which(missing))
  layer_names
}

inflect_kohonen_code_indices <- function(object, code_layers) {
  code_indices <- object$whatmap
  if (is.null(code_indices)) {
    code_indices <- seq_along(code_layers)
  }
  code_indices <- as.integer(code_indices)
  code_indices <- code_indices[!is.na(code_indices)]
  code_indices <- code_indices[code_indices >= 1 & code_indices <= length(code_layers)]
  if (length(code_indices) == 0) {
    stop("Kohonen object does not reference any usable code layers", call. = FALSE)
  }
  unique(code_indices)
}

inflect_bind_kohonen_codes <- function(code_layers,
                                       code_indices,
                                       layer_names,
                                       distance_weights,
                                       user_weights) {
  matrices <- lapply(code_indices, function(layer_index) {
    layer_name <- layer_names[[layer_index]]
    matrix <- inflect_matrix(
      code_layers[[layer_index]],
      paste0("kohonen code layer `", layer_name, "`")
    )
    if (!is.numeric(matrix)) {
      return(NULL)
    }
    storage.mode(matrix) <- "double"
    matrix <- matrix * sqrt(inflect_layer_weight(distance_weights, layer_index)) *
      sqrt(inflect_layer_weight(user_weights, layer_index))

    columns <- colnames(matrix)
    if (is.null(columns)) {
      columns <- paste0("V", seq_len(ncol(matrix)))
    }
    if (length(code_indices) > 1L) {
      columns <- paste0(layer_name, ".", columns)
    }
    colnames(matrix) <- columns
    matrix
  })
  matrices <- Filter(Negate(is.null), matrices)

  if (length(matrices) == 0) {
    stop("Kohonen object must contain at least one numeric code layer", call. = FALSE)
  }

  n_rows <- vapply(matrices, nrow, integer(1))
  if (length(unique(n_rows)) != 1L) {
    stop("Kohonen code layers must have the same number of rows", call. = FALSE)
  }

  codes <- do.call(cbind, matrices)
  colnames(codes) <- make.unique(colnames(codes))
  codes
}

inflect_layer_weight <- function(weights, layer_index) {
  if (is.null(weights) || length(weights) < layer_index || is.na(weights[[layer_index]])) {
    return(1)
  }
  weight <- as.numeric(weights[[layer_index]])
  if (!is.finite(weight) || weight < 0) {
    return(1)
  }
  weight
}
