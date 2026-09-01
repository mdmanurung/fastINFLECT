#' Score marker-level QC criteria across metaclustering schedules
#'
#' @description
#' Scores every requested metaclustering using separate Hartigan dip-test and
#' IQR-spread criteria. Each distinct SOM-node subtree is evaluated once and the
#' per-k matrices are assembled from that cache. The aggregate
#' `qc_pass_rate` is explicitly tied to the criterion selected by
#' `uniform.test`; it is not evidence that a distribution is truly unimodal.
#' All selected marker values must be finite; negative and zero values are
#' retained unchanged.
#'
#' @param FlowSOM.results A supported SOM object with completed SOM clustering.
#'   Supports \pkg{FlowSOM} objects and \pkg{kohonen} objects returned by
#'   \code{\link[kohonen]{som}} or \code{\link[kohonen]{xyf}}.
#' @param metaclustering.list Named list containing exactly one node-label
#'   vector for every value in `set.i`.
#' @param set.i Literal, unique, strictly increasing integer cluster counts.
#' @param workers Number of QC workers. `1L` is serial; values above one use
#'   fork-based \code{\link[parallel]{mclapply}} where supported.
#' @param markers Marker names to score. `NULL` uses the SOM clustering markers.
#' @param uniform.test Aggregate criterion: `"both"` selects the combined dip
#'   and IQR pass, `"spread"` selects IQR only, and `"unimodality"` selects the
#'   dip test only. Both component tests are always retained in `qc.details`.
#' @param th.pvalue Dip-test pass threshold. A cell passes when
#'   `p_value >= th.pvalue`.
#' @param th.IQR IQR-spread pass threshold. A cell passes when `IQR < th.IQR`.
#' @param max.n.diptest Optional positive dip-test sample cap of at least four.
#'   Sampling is deterministic per subtree and marker.
#' @param max.events.per.node Optional positive integer. Each SOM node is
#'   sampled once before subtree assembly; the retained rows feed every marker
#'   and both tests.
#' @param seed Non-negative base seed.
#' @param progress Show stage messages and a serial progress bar. Defaults to
#'   `interactive()`.
#'
#' @return A list containing canonical `scores`, separate named lists of
#'   `dip_pass`, `iqr_pass`, `combined_pass`, and selected `criterion_pass`
#'   matrices, full `qc.details`, and `provenance`. Deprecated `U.set` and
#'   `Accuracy.matrixes` aliases are retained for compatibility.
#' @seealso \code{\link{INFLECT}}, \code{\link{iteration.metacluster}},
#'   \code{\link{FlowSOMQC}}
#' @export
iteration.QC <- function(FlowSOM.results,
                         metaclustering.list,
                         set.i,
                         workers = 1L,
                         markers = NULL,
                         uniform.test = c("both", "spread", "unimodality"),
                         th.pvalue = 0.05,
                         th.IQR = 2,
                         max.n.diptest = NULL,
                         max.events.per.node = NULL,
                         seed = 1L,
                         progress = interactive()) {
  uniform.test <- match.arg(uniform.test)
  progress <- .inflect_validate_progress(progress)
  workers <- .inflect_validate_workers(workers)
  .inflect_validate_qc_arguments(th.pvalue, th.IQR)
  max.n.diptest <- .inflect_validate_max_n_diptest(max.n.diptest)
  max.events.per.node <- .inflect_validate_max_events_per_node(max.events.per.node)
  seed <- .inflect_validate_seed(seed)
  .inflect_warn_sampling(
    max.n.diptest = max.n.diptest,
    max.events.per.node = max.events.per.node
  )
  view <- as_inflect_som(FlowSOM.results)
  n_nodes <- as.integer(view$map$nNodes)
  set.i <- .inflect_validate_qc_schedule(set.i, n_nodes)
  metaclustering.list <- .inflect_validate_metaclustering_list(
    metaclustering.list,
    set.i,
    n_nodes
  )

  prep <- .inflect_prepare_qc(
    view = view,
    markers = markers
  )
  p_of <- .inflect_make_p_of()

  all_node_events <- split(
    seq_len(nrow(prep$data)),
    factor(prep$mapping, levels = seq_len(n_nodes))
  )
  sampled_nodes <- .inflect_sample_node_events(
    node_events = all_node_events,
    max.events.per.node = max.events.per.node,
    seed = seed
  )
  node_events <- sampled_nodes$node_events
  event_sampling <- sampled_nodes$provenance
  rm(all_node_events, sampled_nodes)

  keys_per_k <- vector("list", length(set.i))
  names(keys_per_k) <- as.character(set.i)
  unique_members <- list()
  for (k in set.i) {
    mc <- metaclustering.list[[as.character(k)]]
    keys <- character(k)
    for (cluster in seq_len(k)) {
      nodes <- which(mc == cluster)
      key <- paste0(nodes, collapse = ",")
      keys[[cluster]] <- key
      if (is.null(unique_members[[key]])) {
        unique_members[[key]] <- nodes
      }
    }
    keys_per_k[[as.character(k)]] <- keys
  }

  subtree_keys <- names(unique_members)
  parallel_plan <- .inflect_parallel_plan(
    workers = workers,
    n_tasks = length(subtree_keys)
  )
  if (identical(parallel_plan$backend, "serial_platform_fallback")) {
    warning(
      "Parallel QC is unavailable on this platform; using one worker.",
      call. = FALSE
    )
  }
  .inflect_progress_message(
    progress,
    "[fastINFLECT] Scoring ",
    length(subtree_keys),
    " unique subtrees with ",
    parallel_plan$effective_workers,
    " worker(s)"
  )
  eval_one <- function(idx) {
    subtree_key <- subtree_keys[[idx]]
    nodes <- unique_members[[subtree_key]]
    rows <- if (length(nodes) > 0L) {
      unlist(node_events[nodes], use.names = FALSE)
    } else {
      integer(0)
    }
    .inflect_qc_row_indexed(
      data = prep$data,
      rows = rows,
      uniform.test = uniform.test,
      th.pvalue = th.pvalue,
      th.IQR = th.IQR,
      p_of = p_of,
      subsample = max.n.diptest,
      seed = seed,
      subtree_key = subtree_key,
      marker_indices = prep$marker_indices,
      marker_names = prep$ordered.markers
    )
  }

  if (parallel_plan$use_parallel) {
    rows_list <- withCallingHandlers(
      parallel::mclapply(
        seq_along(subtree_keys),
        eval_one,
        mc.cores = parallel_plan$effective_workers
      ),
      warning = function(w) {
        if (grepl("scheduled core.*encountered error", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  } else {
    progress_bar <- if (progress) {
      utils::txtProgressBar(min = 0, max = length(subtree_keys), style = 3)
    } else {
      NULL
    }
    on.exit(if (!is.null(progress_bar)) close(progress_bar), add = TRUE)
    rows_list <- vector("list", length(subtree_keys))
    for (idx in seq_along(subtree_keys)) {
      rows_list[[idx]] <- eval_one(idx)
      if (!is.null(progress_bar)) {
        utils::setTxtProgressBar(progress_bar, idx)
      }
    }
    if (!is.null(progress_bar)) {
      close(progress_bar)
      progress_bar <- NULL
    }
  }
  rows_list <- .inflect_validate_worker_results(
    rows_list = rows_list,
    keys = subtree_keys,
    marker_names = prep$ordered.markers,
    parallel = parallel_plan$use_parallel
  )
  names(rows_list) <- subtree_keys
  rm(node_events, unique_members)

  qc.details <- vector("list", length(set.i))
  names(qc.details) <- as.character(set.i)
  for (idx in seq_along(set.i)) {
    k <- set.i[[idx]]
    keys <- keys_per_k[[as.character(k)]]
    qc.details[[idx]] <- .inflect_qc_rows_to_detail(
      rows = rows_list[keys],
      marker_names = prep$ordered.markers,
      cluster_names = as.character(seq_len(k))
    )
  }
  rm(rows_list, keys_per_k)

  dip_pass <- lapply(qc.details, `[[`, "dip_pass")
  iqr_pass <- lapply(qc.details, `[[`, "iqr_pass")
  combined_pass <- lapply(qc.details, `[[`, "combined_pass")
  criterion_pass <- lapply(qc.details, `[[`, "criterion_pass")
  criterion <- .inflect_criterion(uniform.test)
  criterion_summary <- .inflect_qc_summary(qc.details, criterion)
  selected_summary <- criterion_summary[criterion_summary$selected, , drop = FALSE]
  scores <- data.frame(
    k = as.integer(selected_summary$k),
    qc_pass_rate = selected_summary$qc_pass_rate,
    criterion = criterion,
    stringsAsFactors = FALSE
  )

  ## Deprecated aliases are deliberately isolated from the canonical contract.
  U.set <- data.frame(
    i = scores$k,
    Unimodality = scores$qc_pass_rate
  )
  source <- view$inflect_source
  provenance <- list(
    criterion = criterion,
    criterion_label = .inflect_criterion_label(criterion),
    uniform.test = uniform.test,
    thresholds = list(
      dip_p_value = th.pvalue,
      iqr = th.IQR
    ),
    th.pvalue = th.pvalue,
    th.IQR = th.IQR,
    markers = prep$ordered.markers,
    marker_indices = prep$marker_indices,
    marker_selection = if (is.null(markers)) "clustering_markers" else "explicit",
    value_handling = list(
      rule = "require finite QC marker data and retain all values unchanged",
      validation = "passed"
    ),
    sampling_order = c(
      "sample each SOM node once",
      "assemble subtrees in ascending SOM-node order",
      "select one marker at a time",
      "require finite selected marker values",
      "retain negative, zero, and positive values unchanged",
      "apply marker-specific dip cap if requested"
    ),
    seeds = list(
      base = seed,
      derivation = "overflow-safe hash of subtree, marker, and RNG stream"
    ),
    seed = seed,
    max.n.diptest = if (is.null(max.n.diptest)) NA_integer_ else max.n.diptest,
    max.events.per.node = if (is.null(max.events.per.node)) {
      NA_integer_
    } else {
      max.events.per.node
    },
    event_sampling = event_sampling,
    n_events = event_sampling$original_events,
    n_events_retained = event_sampling$retained_events,
    parallel = parallel_plan,
    effective_workers = parallel_plan$effective_workers,
    model = list(
      som_type = source$type,
      som_class = source$class,
      data_layer = source$data_layer,
      code_layers = source$code_layers,
      user_weights = source$user_weights,
      distance_weights = source$distance_weights,
      effective_layer_weights = source$effective_layer_weights
    ),
    scoring_mode = if (is.null(max.events.per.node)) {
      "indexed_all_events"
    } else {
      "indexed_node_capped"
    },
    criterion_summary = criterion_summary
  )

  .inflect_progress_message(progress, "[fastINFLECT] QC scoring complete")

  list(
    scores = scores,
    dip_pass = dip_pass,
    iqr_pass = iqr_pass,
    combined_pass = combined_pass,
    criterion_pass = criterion_pass,
    qc.details = qc.details,
    provenance = provenance,
    U.set = U.set,
    Accuracy.matrixes = criterion_pass
  )
}
