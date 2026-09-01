#' Score one SOM metaclustering with dip and IQR criteria
#'
#' @description
#' Computes separate Hartigan dip-test and IQR-spread evidence for every
#' cluster-marker pair. The returned matrix contains the pass decision selected
#' by `uniform.test`; its `qc.details` attribute contains `dip_pass`,
#' `iqr_pass`, `combined_pass`, p-values, IQRs, counts, and failure
#' reasons. All selected marker values must be finite; negative and zero values
#' are retained unchanged. A pass is criterion-specific and does not prove true
#' unimodality.
#'
#' @param FlowSOM.results A supported \pkg{FlowSOM} or \pkg{kohonen} SOM object.
#' @param metaclustering Integer vector with one metacluster label per SOM node.
#' @param markers Marker names to score. `NULL` uses the SOM clustering markers.
#' @param uniform.test Aggregate criterion: `"both"` (dip and IQR), `"spread"`
#'   (IQR), or `"unimodality"` (dip).
#' @param th.pvalue Dip-test pass threshold.
#' @param th.IQR IQR pass threshold.
#' @param max.n.diptest Optional dip-test sample cap of at least four.
#' @param seed Non-negative seed. Simulated dip p-values and optional
#'   subsampling use deterministic subtree-marker streams and preserve the
#'   caller's RNG state.
#' @param progress Show scoring progress. Defaults to `interactive()`.
#'
#' @return Invisibly, the selected-criterion logical matrix. Attributes
#'   `qc.details` and `provenance` retain the separated evidence.
#' @seealso \code{\link{INFLECT}}, \code{\link{iteration.QC}}
#' @export
FlowSOMQC <- function(FlowSOM.results,
                      metaclustering,
                      markers = NULL,
                      uniform.test = c("both", "spread", "unimodality"),
                      th.pvalue = 0.05,
                      th.IQR = 2,
                      max.n.diptest = NULL,
                      seed = 1L,
                      progress = interactive()) {
  uniform.test <- match.arg(uniform.test)
  progress <- .inflect_validate_progress(progress)
  .inflect_validate_qc_arguments(th.pvalue, th.IQR)
  max.n.diptest <- .inflect_validate_max_n_diptest(max.n.diptest)
  seed <- .inflect_validate_seed(seed)
  .inflect_warn_sampling(max.n.diptest = max.n.diptest)

  if (is.null(FlowSOM.results)) {
    stop("`FlowSOM.results` cannot be NULL.", call. = FALSE)
  }
  FlowSOM.results <- as_inflect_som(FlowSOM.results)
  if (is.null(metaclustering)) {
    stop("`metaclustering` cannot be NULL.", call. = FALSE)
  }
  metaclustering <- .inflect_normalize_metaclustering(
    metaclustering = metaclustering,
    n_nodes = FlowSOM.results$map$nNodes,
    label = "`metaclustering`"
  )

  prep <- .inflect_prepare_qc(
    view = FlowSOM.results,
    markers = markers
  )
  p_of <- .inflect_make_p_of()
  event_cluster <- metaclustering[prep$mapping]
  cluster_rows <- split(
    seq_len(nrow(prep$data)),
    factor(event_cluster, levels = seq_len(max(metaclustering)))
  )

  .inflect_progress_message(
    progress,
    "[fastINFLECT] Scoring ",
    length(cluster_rows),
    " metaclusters"
  )
  progress_bar <- if (progress) {
    utils::txtProgressBar(min = 0, max = length(cluster_rows), style = 3)
  } else {
    NULL
  }
  on.exit(if (!is.null(progress_bar)) close(progress_bar), add = TRUE)
  rows <- vector("list", length(cluster_rows))
  for (cluster in seq_along(cluster_rows)) {
    nodes <- which(metaclustering == cluster)
    subtree_key <- paste0(nodes, collapse = ",")
    rows[[cluster]] <- .inflect_qc_row_indexed(
      data = prep$data,
      rows = cluster_rows[[cluster]],
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
    if (!is.null(progress_bar)) {
      utils::setTxtProgressBar(progress_bar, cluster)
    }
  }
  if (!is.null(progress_bar)) {
    close(progress_bar)
    progress_bar <- NULL
  }
  details <- .inflect_qc_rows_to_detail(
    rows,
    prep$ordered.markers,
    as.character(seq_along(cluster_rows))
  )
  criterion <- .inflect_criterion(uniform.test)
  result <- details$criterion_pass
  attr(result, "qc.details") <- details
  attr(result, "provenance") <- list(
    criterion = criterion,
    criterion_label = .inflect_criterion_label(criterion),
    uniform.test = uniform.test,
    thresholds = list(dip_p_value = th.pvalue, iqr = th.IQR),
    markers = prep$ordered.markers,
    marker_selection = if (is.null(markers)) "clustering_markers" else "explicit",
    value_handling = list(
      rule = "require finite QC marker data and retain all values unchanged",
      validation = "passed"
    ),
    seed = seed,
    max.n.diptest = if (is.null(max.n.diptest)) NA_integer_ else max.n.diptest
  )

  .inflect_progress_message(progress, "[fastINFLECT] QC scoring complete")
  invisible(result)
}
