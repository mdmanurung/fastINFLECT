#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not resolve the synthesis script path.", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
output_dir <- file.path(
  repo_root,
  "inst",
  "benchmarks",
  "real-model-selection"
)
paths <- list(
  selection = file.path(output_dir, "selection-evidence.rds"),
  primary = file.path(
    repo_root,
    "inst",
    "benchmarks",
    "real-model-modality",
    "modality-evidence.rds"
  ),
  selected = file.path(
    repo_root,
    "inst",
    "benchmarks",
    "selected-k-modality",
    "modality-evidence.rds"
  ),
  warning = file.path(
    output_dir,
    "curve-fit-warning-diagnostic.rds"
  ),
  reshard = file.path(
    repo_root,
    "inst",
    "benchmarks",
    "real-model-modality",
    "refinement-reshard-provenance.tsv"
  ),
  comparison = file.path(output_dir, "audited-inflect-k-comparison.csv"),
  comparators = file.path(output_dir, "comparator-solution-quality.csv"),
  synthesis = file.path(output_dir, "bmv-k-selection-synthesis.rds"),
  report = file.path(output_dir, "README.md")
)

input_names <- c("selection", "primary", "selected", "warning", "reshard")
missing_inputs <- input_names[!file.exists(unlist(paths[input_names]))]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing final evidence: ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

atomic_save_rds <- function(object, path) {
  temporary <- tempfile("synthesis-", tmpdir = dirname(path))
  saveRDS(object, temporary, compress = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not publish ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  temporary <- tempfile("synthesis-", tmpdir = dirname(path))
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    stop("Could not publish ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_write_lines <- function(lines, path) {
  temporary <- tempfile("synthesis-", tmpdir = dirname(path))
  writeLines(lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not publish ", path, call. = FALSE)
  }
  invisible(path)
}

selection <- readRDS(paths$selection)
primary <- readRDS(paths$primary)
selected <- readRDS(paths$selected)
curve_warning <- readRDS(paths$warning)

stopifnot(
  selection$n_events == 39050953L,
  selection$n_nodes == 900L,
  identical(selection$schedule, 25L:100L),
  identical(
    selection$full_inflect$derived_metric_provenance$version,
    "exclude_zero_event_nodes_v2"
  ),
  all(is.finite(
    selection$full_inflect$partition_metrics$explained_event_weighted
  )),
  isTRUE(primary$completeness$every_partition_materialized),
  isTRUE(selected$completeness$every_partition_materialized),
  primary$completeness$expected_base_combinations ==
    primary$completeness$observed_base_combinations,
  primary$completeness$expected_refinement_combinations ==
    primary$completeness$observed_refinement_combinations,
  selected$completeness$expected_base_combinations ==
    selected$completeness$observed_base_combinations,
  selected$completeness$expected_refinement_combinations ==
    selected$completeness$observed_refinement_combinations,
  identical(
    curve_warning$classification,
    "transient_LL4_optimizer_exploration_warning"
  ),
  identical(
    curve_warning$selection_md5,
    unname(tools::md5sum(paths$selection))
  )
)
model_md5 <- c(
  selection = selection$model_md5,
  primary = primary$stage$model_md5,
  selected = selected$stage$model_md5
)
if (length(unique(unname(model_md5))) != 1L) {
  stop("Evidence bundles do not share one model MD5.", call. = FALSE)
}

all_classifications <- rbind(primary$classifications, selected$classifications)
all_solution_claims <- rbind(primary$solution_claims, selected$solution_claims)
all_status_fractions <- rbind(
  primary$summaries$status_fraction,
  selected$summaries$status_fraction
)
all_sample_manifest <- rbind(
  primary$stage$sample_manifest,
  selected$stage$sample_manifest
)

status_levels <- c(
  "detected_multimodality",
  "ambiguous",
  "no_detected_multimodality",
  "unresolved_discrete"
)
status_counts <- as.data.frame.matrix(xtabs(
  ~ solution_id + status,
  data = all_classifications
))
status_counts$solution_id <- rownames(status_counts)
rownames(status_counts) <- NULL
for (status in status_levels) {
  if (!status %in% names(status_counts)) {
    status_counts[[status]] <- 0L
  }
}
status_counts$total_cluster_marker_pairs <- rowSums(
  status_counts[status_levels]
)
status_counts$detected_fraction <- status_counts$detected_multimodality /
  status_counts$total_cluster_marker_pairs
status_counts$ambiguous_fraction <- status_counts$ambiguous /
  status_counts$total_cluster_marker_pairs
status_counts$unresolved_fraction <- status_counts$unresolved_discrete /
  status_counts$total_cluster_marker_pairs

occupancy_rows <- unique(all_sample_manifest[c(
  "solution_id",
  "cluster",
  "cluster_events"
)])
occupancy <- do.call(rbind, lapply(
  split(occupancy_rows, occupancy_rows$solution_id),
  function(x) {
    positive <- x$cluster_events > 0L
    data.frame(
      solution_id = x$solution_id[[1L]],
      nominal_k = nrow(x),
      event_populated_clusters = sum(positive),
      empty_event_clusters = sum(!positive),
      smallest_nonempty_cluster_events = min(x$cluster_events[positive]),
      largest_cluster_event_fraction = max(x$cluster_events) /
        selection$n_events,
      stringsAsFactors = FALSE
    )
  }
))
rownames(occupancy) <- NULL

criterion <- selection$full_inflect$criterion_summary
criterion_wide <- reshape(
  criterion[c("k", "criterion", "qc_pass_rate", "unresolved")],
  idvar = "k",
  timevar = "criterion",
  direction = "wide"
)
names(criterion_wide) <- sub(
  "^qc_pass_rate[.]",
  "",
  names(criterion_wide)
)
names(criterion_wide) <- sub(
  "^unresolved[.]",
  "unresolved_",
  names(criterion_wide)
)
names(criterion_wide)[names(criterion_wide) %in% c(
  "dip",
  "iqr",
  "combined"
)] <- paste0(
  names(criterion_wide)[names(criterion_wide) %in% c(
    "dip",
    "iqr",
    "combined"
  )],
  "_pass_rate"
)

metrics <- selection$full_inflect$partition_metrics
metrics <- metrics[
  metrics$code_space == "model_effective_codes",
  ,
  drop = FALSE
]
metric_occupancy <- metrics[c(
  "k",
  "event_populated_clusters",
  "empty_event_clusters"
)]
metrics$event_populated_clusters <- NULL
metrics$empty_event_clusters <- NULL
solution_ids <- unique(all_classifications$solution_id)
effective_ids <- grep(
  "^inflect_ward_k[0-9]+$",
  solution_ids,
  value = TRUE
)
effective <- data.frame(
  solution_id = effective_ids,
  k = as.integer(sub("^inflect_ward_k", "", effective_ids)),
  stringsAsFactors = FALSE
)
effective <- merge(effective, criterion_wide, by = "k", all.x = TRUE)
effective <- merge(effective, metrics, by = "k", all.x = TRUE)
effective <- merge(effective, occupancy, by = "solution_id", all.x = TRUE)
effective <- merge(effective, status_counts, by = "solution_id", all.x = TRUE)
effective <- merge(
  effective,
  all_solution_claims[c(
    "solution_id",
    "strict_solution_no_detected_multimodality",
    "statement"
  )],
  by = "solution_id",
  all.x = TRUE
)
effective <- effective[order(effective$k), ]
rownames(effective) <- NULL
matched_metric_occupancy <- metric_occupancy[
  match(effective$k, metric_occupancy$k),
  ,
  drop = FALSE
]
if (
  anyNA(matched_metric_occupancy$k) ||
    !identical(
      as.integer(effective$event_populated_clusters),
      as.integer(matched_metric_occupancy$event_populated_clusters)
    ) ||
    !identical(
      as.integer(effective$empty_event_clusters),
      as.integer(matched_metric_occupancy$empty_event_clusters)
    )
) {
  stop(
    "Independent event-occupancy summaries do not agree.",
    call. = FALSE
  )
}

sample_iqr <- selection$sample_iqr_summary
comparators <- merge(
  all_solution_claims,
  occupancy,
  by = "solution_id",
  all.x = TRUE,
  sort = FALSE
)
comparators <- merge(
  comparators,
  status_counts,
  by = "solution_id",
  all.x = TRUE,
  sort = FALSE
)
comparators <- merge(
  comparators,
  sample_iqr,
  by = "solution_id",
  all.x = TRUE,
  sort = FALSE
)
comparators <- comparators[order(comparators$family, comparators$k), ]
rownames(comparators) <- NULL

selections <- selection$full_inflect$criterion_selections
selected_k <- function(criterion_name, method_name) {
  value <- selections$k[
    selections$criterion == criterion_name &
      selections$method == method_name
  ]
  if (length(value) != 1L) {
    stop("Expected one criterion/method selection.", call. = FALSE)
  }
  as.integer(value)
}
combined_inflection <- selected_k("combined", "inflection")
combined_kneedle <- selected_k("combined", "kneedle")
dip_inflection <- selected_k("dip", "inflection")
dip_kneedle <- selected_k("dip", "kneedle")
spread_inflection <- selected_k("iqr", "inflection")
spread_kneedle <- selected_k("iqr", "kneedle")
operational_k <- combined_kneedle
operational_quality <- effective[effective$k == operational_k, , drop = FALSE]
if (nrow(operational_quality) != 1L) {
  stop("The operational combined Kneedle was not directly audited.", call. = FALSE)
}

consensus_plateau <- selection$consensus$plateau$k[[1L]]
consensus_at_operational <- selection$consensus$summary[
  selection$consensus$summary$k == operational_k,
  ,
  drop = FALSE
]
consensus_at_spread <- selection$consensus$summary[
  selection$consensus$summary$k == spread_kneedle,
  ,
  drop = FALSE
]

metric_recommendations <- selection$full_inflect$metric_recommendations
metric_k <- function(code_space, method) {
  value <- metric_recommendations$k[
    metric_recommendations$code_space == code_space &
      metric_recommendations$method == method
  ]
  if (length(value) != 1L) {
    stop("Expected one internal-metric recommendation.", call. = FALSE)
  }
  as.integer(value)
}
effective_metric_k <- c(
  silhouette = metric_k("model_effective_codes", "silhouette_max"),
  calinski_harabasz = metric_k(
    "model_effective_codes",
    "calinski_harabasz_max"
  ),
  davies_bouldin = metric_k(
    "model_effective_codes",
    "davies_bouldin_min"
  ),
  explained_unweighted = metric_k(
    "model_effective_codes",
    "explained_unweighted_kneedle"
  ),
  explained_event_weighted = metric_k(
    "model_effective_codes",
    "explained_event_weighted_kneedle"
  )
)
historical_metric_k <- c(
  silhouette = metric_k("historical_80_20_codes", "silhouette_max"),
  calinski_harabasz = metric_k(
    "historical_80_20_codes",
    "calinski_harabasz_max"
  ),
  davies_bouldin = metric_k(
    "historical_80_20_codes",
    "davies_bouldin_min"
  ),
  explained_unweighted = metric_k(
    "historical_80_20_codes",
    "explained_unweighted_kneedle"
  ),
  explained_event_weighted = metric_k(
    "historical_80_20_codes",
    "explained_event_weighted_kneedle"
  )
)

weighted_status <- all_status_fractions[
  all_status_fractions$solution_id ==
    paste0("inflect_ward_k", operational_k),
  ,
  drop = FALSE
]

synthesis <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
  decision = list(
    operational_k = operational_k,
    lower_resolution_inflection_k = combined_inflection,
    spread_sensitivity_upper_k = spread_kneedle,
    status = paste0(
      "Operational metaclustering recommendation, not a claim of biological ",
      "truth or global unimodality."
    )
  ),
  inflect_selections = list(
    combined = c(inflection = combined_inflection, kneedle = combined_kneedle),
    dip = c(inflection = dip_inflection, kneedle = dip_kneedle),
    iqr_spread = c(
      inflection = spread_inflection,
      kneedle = spread_kneedle
    )
  ),
  consensus = list(
    plateau = selection$consensus$plateau,
    operational_k = consensus_at_operational,
    spread_k = consensus_at_spread
  ),
  internal_metric_recommendations =
    metric_recommendations,
  audited_inflect_k = effective,
  comparator_solutions = comparators,
  operational_status_fractions = weighted_status,
  curve_fit_warning = curve_warning,
  model_md5 = model_md5,
  input_md5 = unname(tools::md5sum(c(
    paths$selection,
    paths$primary,
    paths$selected,
    paths$warning,
    paths$reshard
  ))),
  input_files = normalizePath(c(
    paths$selection,
    paths$primary,
    paths$selected,
    paths$warning,
    paths$reshard
  ), mustWork = TRUE),
  script_file = script_file,
  script_md5 = unname(tools::md5sum(script_file)),
  session_info = utils::sessionInfo()
)

format_value <- function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  format(round(x, digits), trim = TRUE, scientific = FALSE)
}
plateau_text <- if (is.na(consensus_plateau)) {
  "No five-k consensus plateau met both delta-area <= 0.01 and minimum seed ARI >= 0.95."
} else {
  paste0(
    "The first five-k consensus plateau starts at k=",
    consensus_plateau,
    "."
  )
}
quality <- operational_quality
global_modality_text <- if (any(
  effective$strict_solution_no_detected_multimodality
)) {
  paste0(
    "At least one directly audited k met the workflow's strict no-detected ",
    "criterion; this remains absence of detected univariate multimodality, ",
    "not proof of true unimodality."
  )
} else {
  "No directly audited k supports a global unimodality claim."
}
report <- c(
  "# BMV cluster-number synthesis",
  "",
  paste0("- Generated: ", synthesis$generated_at),
  paste0("- Full model: ", format(selection$n_events, big.mark = ","), " events, 900 SOM nodes, 27 markers."),
  paste0("- Model MD5: `", unname(model_md5[[1L]]), "`."),
  "",
  "## Recommendation",
  "",
  paste0(
    "- Use nominal **k=", operational_k,
    "** as the operational fitted-model partition. It is the combined and dip ",
    "Kneedle and was directly materialised and independently audited."
  ),
  paste0(
    "- Treat k=", combined_inflection,
    " as the lower-resolution parametric elbow and k=", spread_kneedle,
    " as the spread/variance sensitivity endpoint; the defensible sensitivity ",
    "range is therefore ", operational_k, "-", spread_kneedle, "."
  ),
  paste0(
    "- k=", operational_k, " contains ",
    quality$event_populated_clusters, " event-populated clusters and ",
    quality$empty_event_clusters, " empty nominal cluster; its largest cluster ",
    "contains ", format_value(100 * quality$largest_cluster_event_fraction, 2L),
    "% of events."
  ),
  paste0(
    "- Final modality audit at k=", operational_k, ": ",
    quality$detected_multimodality, " detected, ",
    quality$ambiguous, " ambiguous, ",
    quality$unresolved_discrete, " unresolved, and ",
    quality$no_detected_multimodality,
    " no-detected cluster-marker entries."
  ),
  paste0("- ", global_modality_text),
  "",
  "## Method comparison",
  "",
  paste0(
    "- Combined INFLECT: inflection k=", combined_inflection,
    "; Kneedle k=", combined_kneedle, "."
  ),
  paste0(
    "- Dip-only: inflection k=", dip_inflection,
    "; Kneedle k=", dip_kneedle, "."
  ),
  paste0(
    "- IQR spread: inflection k=", spread_inflection,
    "; Kneedle k=", spread_kneedle, "."
  ),
  paste0("- Consensus metaclustering: ", plateau_text),
  paste0(
    "- At k=", operational_k, ", consensus minimum seed ARI=",
    format_value(consensus_at_operational$min_seed_ari),
    " and median relative delta area=",
    format_value(consensus_at_operational$median_relative_delta_area, 4L),
    "."
  ),
  paste0(
    "- Fitted-model internal metrics: silhouette k=",
    effective_metric_k[["silhouette"]], ", Calinski-Harabasz k=",
    effective_metric_k[["calinski_harabasz"]], ", Davies-Bouldin k=",
    effective_metric_k[["davies_bouldin"]], ", unweighted dispersion k=",
    effective_metric_k[["explained_unweighted"]],
    ", event-weighted dispersion k=",
    effective_metric_k[["explained_event_weighted"]], "."
  ),
  paste0(
    "- Historical intended 80/20 code space: silhouette k=",
    historical_metric_k[["silhouette"]], ", Calinski-Harabasz k=",
    historical_metric_k[["calinski_harabasz"]], ", Davies-Bouldin k=",
    historical_metric_k[["davies_bouldin"]],
    ", unweighted dispersion k=",
    historical_metric_k[["explained_unweighted"]],
    ", event-weighted dispersion k=",
    historical_metric_k[["explained_event_weighted"]], "."
  ),
  paste0(
    "- The captured `NaNs produced` warning was reproducibly classified as ",
    "temporary LL.4 optimizer exploration: all final fitted values were finite ",
    "and the rerun inflections exactly matched 28/38/27."
  ),
  "- Silhouette, Calinski-Harabasz, Davies-Bouldin, and explained-dispersion recommendations are retained in `metric-recommendations.csv`; they diagnose different trade-offs and are not interchangeable votes.",
  "",
  "## Interpretation limits",
  "",
  "- IQR spread is a dispersion screen, not a modality test.",
  "- Dip plus ACR evidence addresses univariate marker marginals with BH correction within each 27-marker cluster/method/seed/sample-size family.",
  "- Consensus stability measures reproducibility conditional on k; it does not establish biological correctness or select k by itself.",
  "- The fitted model's effective layer weighting differs from the historical intended 80/20 code space, so the historical k=50/62 partitions remain sensitivity comparators rather than substitutes for the fitted-model recommendation.",
  "- The slow 5,000-event primary refinement tail was resumed from immutable per-task checkpoints with a wider logical sharding; hashes and job mappings are retained in `refinement-reshard-provenance.tsv`.",
  "",
  "See `audited-inflect-k-comparison.csv` and `comparator-solution-quality.csv` for exact counts and metrics."
)

atomic_write_csv(effective, paths$comparison)
atomic_write_csv(comparators, paths$comparators)
atomic_save_rds(synthesis, paths$synthesis)
atomic_write_lines(report, paths$report)
message("BMV k-selection synthesis complete: ", paths$report)
