normalize_set_i <- function(set.i, FlowSOM.results) {
  FlowSOM.results <- as_inflect_som(FlowSOM.results)
  if (is.null(set.i)) {
    stop("Error in iteration.metacluster: The 'set.i' parameter can not be NULL")
  }
  if (!is.numeric(set.i) || anyNA(set.i)) {
    stop("Error in iteration.metacluster: The 'set.i' parameter must be numeric")
  }

  n_nodes <- as.integer(FlowSOM.results$map$nNodes)
  if (length(n_nodes) != 1L || is.na(n_nodes) || n_nodes < 1L) {
    stop("Error in iteration.metacluster: SOM node count must be a positive integer", call. = FALSE)
  }
  if (n_nodes < 5L) {
    stop("INFLECT requires at least five SOM nodes to fit the diagnostic curve", call. = FALSE)
  }

  bounded_seq <- function(from, to, by = 1L) {
    from <- as.integer(from)
    to <- as.integer(to)
    by <- as.integer(by)
    if (from > to) {
      return(integer(0))
    }
    seq(from, to, by = by)
  }

  if (length(set.i) == 2) {
    max_default_k <- as.integer(floor(0.9 * n_nodes))
    if (length(bounded_seq(5L, max_default_k)) < 5L) {
      max_default_k <- n_nodes
    }

    if (set.i[1] > max_default_k) {
      set.i <- bounded_seq(max(1L, max_default_k - 4L), max_default_k)
    } else if (set.i[2] < max_default_k) {
      set.i <-
        c(
          bounded_seq(5L, set.i[1]),
          bounded_seq(set.i[1] + 5L, set.i[2], 5L),
          bounded_seq(set.i[2] + 10L, max_default_k, 10L)
        )
    } else {
      set.i <-
        c(
          bounded_seq(5L, set.i[1]),
          bounded_seq(set.i[1] + 5L, max_default_k, 5L)
        )
    }

    if (length(set.i) < 5L) {
      set.i <- sort(unique(c(
        set.i,
        bounded_seq(max(1L, n_nodes - 4L), n_nodes)
      )))
    }
  } else if (length(set.i) < 50) {
    warning("Warning: Number of points in set.i is low, beware of noisy diagnostic curves ")
  }

  set.i <- unique(as.integer(set.i))
  if (length(set.i) == 0 || any(set.i < 1)) {
    stop("Error in iteration.metacluster: The 'set.i' parameter must contain positive cluster numbers")
  }
  if (any(set.i > n_nodes)) {
    stop("Error in iteration.metacluster: The 'set.i' parameter must not exceed the number of SOM nodes", call. = FALSE)
  }
  if (length(set.i) < 5L) {
    warning("Warning: Number of points in set.i is low, beware of noisy diagnostic curves ")
  }

  set.i
}

resolve_inflect_markers <- function(FlowSOM.results,
                                    only.clustering.markers = TRUE,
                                    acquired_markers = NULL) {
  FlowSOM.results <- as_inflect_som(FlowSOM.results)
  clustering.markers <- FlowSOM.results$prettyColnames[FlowSOM.results$map$colsUsed]
  if (only.clustering.markers) {
    markers <- clustering.markers
  } else {
    if (!is.null(acquired_markers) && all(acquired_markers %in% FlowSOM.results$prettyColnames)) {
      markers <- acquired_markers
    } else {
      stop("Error in acquired_markers: The 'acquired_markers' vector must match names in 'FlowSOM.result$prettyColnames' ")
    }
  }

  c(
    gtools::mixedsort(intersect(markers, clustering.markers)),
    gtools::mixedsort(setdiff(markers, clustering.markers))
  )
}

build_inflect_provenance <- function(FlowSOM.results,
                                     set.i,
                                     uniform.test,
                                     th.pvalue,
                                     th.IQR,
                                     zeroes.in,
                                     basedata,
                                     only.clustering.markers,
                                     acquired_markers,
                                     markers,
                                     elapsed_seconds,
                                     diptest_args = list(),
                                     max.n.diptest = NULL,
                                     seed = 1L,
                                     target = 0.95,
                                     timestamp = Sys.time()) {
  FlowSOM.results <- as_inflect_som(FlowSOM.results)
  source <- FlowSOM.results$inflect_source
  list(
    set.i = set.i,
    uniform.test = uniform.test,
    th.pvalue = th.pvalue,
    th.IQR = th.IQR,
    zeroes.in = zeroes.in,
    basedata = basedata,
    max.n.diptest = if (is.null(max.n.diptest)) NA_integer_ else as.integer(max.n.diptest),
    seed = seed,
    target = target,
    only.clustering.markers = only.clustering.markers,
    acquired_markers = acquired_markers,
    markers = markers,
    n_events = nrow(FlowSOM.results$data),
    n_nodes = FlowSOM.results$map$nNodes,
    n_clustering_markers = length(FlowSOM.results$map$colsUsed),
    som_type = source$type,
    som_class = source$class,
    som_data_layer = source$data_layer,
    som_code_layers = source$code_layers,
    package_version = inflect_package_version(),
    r_version = R.version.string,
    timestamp = format(timestamp, "%Y-%m-%dT%H:%M:%OS%z"),
    elapsed_seconds = elapsed_seconds,
    diptest_args = diptest_args
  )
}

inflect_package_version <- function() {
  version <- tryCatch(
    as.character(utils::packageVersion("fastINFLECT")),
    error = function(e) NA_character_
  )
  if (!is.na(version)) {
    return(version)
  }

  desc_path <- system.file("DESCRIPTION", package = "fastINFLECT")
  if (nzchar(desc_path)) {
    desc_version <- tryCatch(
      unname(read.dcf(desc_path)[1, "Version"]),
      error = function(e) NA_character_
    )
    if (!is.na(desc_version)) {
      return(desc_version)
    }
  }

  NA_character_
}
