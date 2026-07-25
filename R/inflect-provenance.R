normalize_set_i <- function(set.i, FlowSOM.results) {
  FlowSOM.results <- as_inflect_som(FlowSOM.results)
  if (is.null(set.i)) {
    stop("`set.i` is required and cannot be NULL.", call. = FALSE)
  }
  n_nodes <- as.integer(FlowSOM.results$map$nNodes)
  if (length(n_nodes) != 1L || is.na(n_nodes) || n_nodes < 1L) {
    stop("Error in iteration.metacluster: SOM node count must be a positive integer", call. = FALSE)
  }
  if (n_nodes < 5L) {
    stop("INFLECT requires at least five SOM nodes to fit the diagnostic curve", call. = FALSE)
  }
  set.i <- .inflect_validate_qc_schedule(set.i, n_nodes)
  if (length(set.i) < 5L) {
    stop(
      "`set.i` must contain at least five cluster counts for curve fitting.",
      call. = FALSE
    )
  }

  set.i
}

#' Build an explicitly bounded adaptive k schedule
#'
#' @description
#' Constructs the legacy adaptive schedule explicitly: dense single-k steps at
#' low k, five-k steps at medium k, and ten-k steps at high k. `max_k` is an
#' explicit upper bound, and the returned schedule never exceeds it.
#'
#' @param n_nodes Number of SOM nodes.
#' @param max_k Largest metacluster count to include. Must not exceed `n_nodes`.
#' @param dense_until Last k in the dense, one-k part of the schedule.
#' @param medium_until Last k in the medium, five-k part of the schedule.
#'
#' @return A strictly increasing integer vector suitable for `INFLECT(set.i=)`.
#' @export
inflect_adaptive_set_i <- function(n_nodes,
                                   max_k,
                                   dense_until = 25L,
                                   medium_until = 100L) {
  values <- list(
    n_nodes = n_nodes,
    max_k = max_k,
    dense_until = dense_until,
    medium_until = medium_until
  )
  valid_scalar_integer <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
      x >= -.Machine$integer.max && x <= .Machine$integer.max &&
      x == floor(x)
  }
  if (!all(vapply(values, valid_scalar_integer, logical(1)))) {
    stop(
      "`n_nodes`, `max_k`, `dense_until`, and `medium_until` must be finite integer scalars.",
      call. = FALSE
    )
  }

  n_nodes <- as.integer(n_nodes)
  max_k <- as.integer(max_k)
  dense_until <- as.integer(dense_until)
  medium_until <- as.integer(medium_until)
  if (n_nodes < 5L) {
    stop("`n_nodes` must be at least 5.", call. = FALSE)
  }
  if (max_k < 5L || max_k > n_nodes) {
    stop("`max_k` must be between 5 and `n_nodes`.", call. = FALSE)
  }
  if (dense_until < 5L || medium_until < dense_until) {
    stop(
      "`dense_until` must be at least 5 and no greater than `medium_until`.",
      call. = FALSE
    )
  }

  bounded_seq <- function(from, to, by) {
    from <- as.double(from)
    to <- as.double(to)
    by <- as.double(by)
    if (from > to) {
      return(integer(0))
    }
    length_out <- floor((to - from) / by) + 1
    if (!is.finite(length_out) || length_out > 1e7) {
      stop(
        "The requested adaptive schedule is too large to materialize safely.",
        call. = FALSE
      )
    }
    as.integer(from + seq.int(0, length_out - 1) * by)
  }
  dense_end <- min(as.double(max_k), as.double(dense_until))
  medium_end <- min(as.double(max_k), as.double(medium_until))
  medium_start <- as.double(dense_until) + 5
  high_start <- as.double(medium_until) + 10
  schedule <- c(
    bounded_seq(5, dense_end, 1),
    bounded_seq(medium_start, medium_end, 5),
    bounded_seq(high_start, as.double(max_k), 10)
  )
  schedule <- sort(unique(c(schedule, max_k)))
  if (length(schedule) < 5L) {
    stop(
      "The adaptive schedule must contain at least five values; increase `max_k`.",
      call. = FALSE
    )
  }
  as.integer(schedule)
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
                                     max.events.per.node = NULL,
                                     seed = 1L,
                                     target = 0.95,
                                     qc_provenance = list(),
                                     timestamp = Sys.time()) {
  FlowSOM.results <- as_inflect_som(FlowSOM.results)
  source <- FlowSOM.results$inflect_source
  if (!is.list(qc_provenance)) {
    stop("`qc_provenance` must be a list.", call. = FALSE)
  }
  event_sampling <- qc_provenance$event_sampling
  if (is.null(event_sampling)) {
    node_counts <- as.integer(tabulate(
      FlowSOM.results$map$mapping[, 1],
      nbins = FlowSOM.results$map$nNodes
    ))
    event_sampling <- list(
      mode = "all_events",
      max_events_per_node = NA_integer_,
      original_events = nrow(FlowSOM.results$data),
      retained_events = nrow(FlowSOM.results$data),
      per_node = data.frame(
        node = seq_len(FlowSOM.results$map$nNodes),
        original_events = node_counts,
        retained_events = node_counts
      )
    )
  }
  fit_provenance <- list(
    set.i = set.i,
    basedata = basedata,
    target = target,
    only.clustering.markers = only.clustering.markers,
    acquired_markers = acquired_markers,
    markers = markers,
    n_nodes = FlowSOM.results$map$nNodes,
    n_clustering_markers = length(FlowSOM.results$map$colsUsed),
    som_type = source$type,
    som_class = source$class,
    som_data_layer = source$data_layer,
    som_code_layers = source$code_layers,
    som_user_weights = source$user_weights,
    som_distance_weights = source$distance_weights,
    som_effective_layer_weights = source$effective_layer_weights,
    package_version = inflect_package_version(),
    r_version = R.version.string,
    timestamp = format(timestamp, "%Y-%m-%dT%H:%M:%OS%z"),
    elapsed_seconds = elapsed_seconds
  )
  provenance <- utils::modifyList(qc_provenance, fit_provenance)

  ## Compatibility scalars retained for existing provenance readers. Their
  ## values come from the QC boundary rather than being recomputed here.
  provenance$uniform.test <- uniform.test
  provenance$th.pvalue <- th.pvalue
  provenance$th.IQR <- th.IQR
  provenance$zeroes.in <- zeroes.in
  provenance$max.n.diptest <- if (is.null(max.n.diptest)) {
    NA_integer_
  } else {
    as.integer(max.n.diptest)
  }
  provenance$max.events.per.node <- if (is.null(max.events.per.node)) {
    NA_integer_
  } else {
    as.integer(max.events.per.node)
  }
  provenance$seed <- seed
  provenance$n_events <- event_sampling$original_events
  provenance$n_events_retained <- event_sampling$retained_events
  provenance$event_sampling <- event_sampling
  provenance$diptest_args <- diptest_args
  provenance
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
