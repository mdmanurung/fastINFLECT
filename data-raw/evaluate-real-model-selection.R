#!/usr/bin/env Rscript

## Resumable cluster-number and spread benchmark for the 39,050,953-event BMV
## model. This is deliberately separate from validate-real-model-modality.R:
## the latter tests marker modality, whereas this script compares k-selection
## diagnostics and retains criterion-specific fastINFLECT evidence.

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0L) {
    return(default)
  }
  sub(prefix, "", hit[[length(hit)]], fixed = TRUE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub(
  "^--file=",
  "",
  script_args[startsWith(script_args, "--file=")]
)
repo_root <- if (length(script_file) > 0L) {
  normalizePath(file.path(dirname(script_file[[1L]]), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

mode <- arg_value("mode", "run")
model_path <- Sys.getenv(
  "INFLECT_LARGE_MODEL",
  unset = paste0(
    "/exports/para-lipg-hpc/Xuran/results/bmv/aurora_exvivo/cluster/",
    "model/xyf_som_model_with_LD_filter.rds"
  )
)
selection_work_dir <- Sys.getenv(
  "INFLECT_SELECTION_WORKDIR",
  unset = file.path(repo_root, "data-raw", ".real-model-selection-work")
)
modality_work_dir <- Sys.getenv(
  "INFLECT_MODALITY_WORKDIR",
  unset = file.path(repo_root, "data-raw", ".real-model-modality-work")
)
evidence_dir <- Sys.getenv(
  "INFLECT_SELECTION_EVIDENCE_DIR",
  unset = file.path(repo_root, "inst", "benchmarks", "real-model-selection")
)

schedule <- 25L:100L
candidate_k <- c(50L, 55L, 60L, 62L, 65L)
seeds <- c(1L, 42L, 2026L)
consensus_reps <- 100L
iqr_threshold <- 2
consensus_plateau_delta <- 0.01
consensus_plateau_width <- 5L
derived_metrics_version <- "exclude_zero_event_nodes_v2"

dir.create(selection_work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(evidence_dir, recursive = TRUE, showWarnings = FALSE)
paths <- list(
  contract = file.path(selection_work_dir, "run-contract.rds"),
  full_inflect = file.path(selection_work_dir, "full-inflect-selection.rds"),
  consensus_dir = file.path(selection_work_dir, "consensus"),
  iqr_dir = file.path(selection_work_dir, "sample-iqr"),
  final = file.path(evidence_dir, "selection-evidence.rds")
)
dir.create(paths$consensus_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$iqr_dir, recursive = TRUE, showWarnings = FALSE)

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = "checkpoint-", tmpdir = dirname(path))
  saveRDS(object, temporary, compress = TRUE)
  if (!file.rename(temporary, path)) {
    copied <- file.copy(temporary, path, overwrite = TRUE)
    unlink(temporary)
    if (!isTRUE(copied)) {
      stop("Could not save checkpoint: ", path, call. = FALSE)
    }
  }
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = "table-", tmpdir = dirname(path))
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    copied <- file.copy(temporary, path, overwrite = TRUE)
    unlink(temporary)
    if (!isTRUE(copied)) {
      stop("Could not save CSV: ", path, call. = FALSE)
    }
  }
  invisible(path)
}

read_status_kb <- function(field = "VmHWM") {
  status_path <- "/proc/self/status"
  if (!file.exists(status_path)) {
    return(NA_real_)
  }
  status <- readLines(status_path, warn = FALSE)
  hit <- grep(paste0("^", field, ":"), status, value = TRUE)
  if (length(hit) == 0L) {
    return(NA_real_)
  }
  as.numeric(sub(".*?([0-9]+) kB.*", "\\1", hit[[1L]]))
}

package_version_or_na <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    return(NA_character_)
  }
  as.character(utils::packageVersion(package))
}

selection_source_files <- function() {
  files <- c(
    normalizePath(script_file[[1L]], mustWork = TRUE),
    file.path(repo_root, "data-raw", "run-real-model-consensus-seed.R"),
    file.path(repo_root, "data-raw", "resume-real-model-selection.sbatch"),
    list.files(
      file.path(repo_root, "R"),
      pattern = "[.]R$",
      full.names = TRUE
    ),
    file.path(repo_root, "DESCRIPTION")
  )
  normalizePath(files[file.exists(files)], mustWork = TRUE)
}

selection_contract <- function(stage) {
  source_files <- selection_source_files()
  list(
    schema_version = 3L,
    model_path = stage$model_path,
    model_md5 = stage$model_md5,
    model_size_bytes = stage$model_size_bytes,
    stage_generated_at = stage$generated_at,
    stage_source_md5 = stage$source_provenance$source_md5,
    source_files = normalizePath(source_files, mustWork = TRUE),
    source_md5 = unname(tools::md5sum(source_files)),
    schedule = schedule,
    candidate_k = candidate_k,
    seeds = seeds,
    consensus_reps = consensus_reps,
    consensus_p_item = 0.9,
    consensus_p_feature = 1,
    consensus_plateau_delta = consensus_plateau_delta,
    consensus_plateau_width = consensus_plateau_width,
    iqr_threshold = iqr_threshold,
    sample_iqr_events = 5000L,
    finite_marker_values_required = TRUE,
    max_n_diptest = NA_integer_,
    max_events_per_node = NA_integer_,
    derived_metrics_version = derived_metrics_version
  )
}

ensure_contract <- function(stage) {
  current <- selection_contract(stage)
  if (!file.exists(paths$contract)) {
    atomic_save_rds(current, paths$contract)
    return(current)
  }
  recorded <- readRDS(paths$contract)
  if (!identical(recorded, current)) {
    stop(
      "Selection inputs or source changed. Use a new ",
      "INFLECT_SELECTION_WORKDIR instead of reusing stale checkpoints.",
      call. = FALSE
    )
  }
  recorded
}

adjusted_rand_index <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2L) {
    return(NA_real_)
  }
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  total_pairs <- choose2(n)
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- (sum_rows + sum_cols) / 2
  denominator <- maximum - expected
  if (denominator == 0) {
    return(if (identical(as.integer(x), as.integer(y))) 1 else NA_real_)
  }
  (sum_cells - expected) / denominator
}

scale_codebook <- function(codes) {
  scaled <- scale(codes)
  scaled[is.na(scaled)] <- 0
  scaled
}

weighted_dispersion <- function(codes, labels, weights) {
  codes <- as.matrix(codes)
  labels <- as.vector(labels)
  weights <- as.numeric(weights)
  if (nrow(codes) != length(labels) ||
      nrow(codes) != length(weights) ||
      any(!is.finite(codes)) ||
      any(!is.finite(weights)) ||
      any(weights < 0) ||
      !any(weights > 0)) {
    stop("Invalid codes, labels, or nonnegative weights.", call. = FALSE)
  }
  ## Zero-event SOM nodes remain in the nominal partition and occupancy
  ## diagnostics, but they cannot define an event-weighted centroid. Excluding
  ## them here prevents a zero-only metacluster from producing 0 / 0 centroids.
  positive <- weights > 0
  codes <- codes[positive, , drop = FALSE]
  labels <- labels[positive]
  weights <- weights[positive]
  labels <- as.integer(factor(labels))
  weighted_codes <- codes * weights
  cluster_weight <- rowsum(weights, labels, reorder = TRUE)
  centroids <- rowsum(weighted_codes, labels, reorder = TRUE) /
    as.numeric(cluster_weight)
  centered <- codes - centroids[labels, , drop = FALSE]
  within <- sum(weights * rowSums(centered^2))
  overall <- colSums(weighted_codes) / sum(weights)
  total <- sum(weights * rowSums((codes - rep(
    overall,
    each = nrow(codes)
  ))^2))
  list(
    within = within,
    total = total,
    explained = if (total == 0) NA_real_ else 1 - within / total,
    centroids = centroids
  )
}

davies_bouldin <- function(codes, labels) {
  labels <- as.integer(factor(labels))
  counts <- tabulate(labels)
  centroids <- rowsum(codes, labels, reorder = TRUE) / counts
  scatter <- vapply(
    seq_along(counts),
    function(cluster) {
      rows <- labels == cluster
      mean(sqrt(rowSums(
        (codes[rows, , drop = FALSE] -
          rep(centroids[cluster, ], each = sum(rows)))^2
      )))
    },
    numeric(1)
  )
  centroid_distance <- as.matrix(stats::dist(centroids))
  ratio <- outer(scatter, scatter, "+") / centroid_distance
  diag(ratio) <- -Inf
  if (any(!is.finite(apply(ratio, 1L, max)))) {
    return(NA_real_)
  }
  mean(apply(ratio, 1L, max))
}

partition_metrics <- function(codes,
                              node_weights,
                              code_space,
                              ks = schedule) {
  distance <- stats::dist(codes)
  fit <- stats::hclust(distance, method = "ward.D2")
  unweighted_total <- weighted_dispersion(
    codes,
    rep(1L, nrow(codes)),
    rep(1, nrow(codes))
  )$total
  rows <- vector("list", length(ks))

  for (i in seq_along(ks)) {
    k <- ks[[i]]
    labels <- stats::cutree(fit, k = k)
    unweighted <- weighted_dispersion(
      codes,
      labels,
      rep(1, nrow(codes))
    )
    weighted <- weighted_dispersion(codes, labels, node_weights)
    between <- unweighted_total - unweighted$within
    ch <- (between / (k - 1)) /
      (unweighted$within / (nrow(codes) - k))
    silhouette <- cluster::silhouette(labels, distance)
    cluster_events <- as.numeric(rowsum(node_weights, labels, reorder = TRUE))

    rows[[i]] <- data.frame(
      code_space = code_space,
      k = k,
      silhouette = mean(silhouette[, "sil_width"]),
      calinski_harabasz = ch,
      davies_bouldin = davies_bouldin(codes, labels),
      explained_unweighted = unweighted$explained,
      explained_event_weighted = weighted$explained,
      smallest_cluster_nodes = min(tabulate(labels)),
      smallest_cluster_events = min(cluster_events),
      smallest_cluster_event_fraction = min(cluster_events) /
        sum(cluster_events),
      event_populated_clusters = sum(cluster_events > 0),
      empty_event_clusters = sum(cluster_events == 0),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

metric_recommendations <- function(metrics) {
  output <- list()
  index <- 0L
  for (space in unique(metrics$code_space)) {
    x <- metrics[metrics$code_space == space, , drop = FALSE]
    recommendations <- list(
      silhouette_max = x$k[[which.max(x$silhouette)]],
      calinski_harabasz_max = x$k[[which.max(x$calinski_harabasz)]],
      davies_bouldin_min = x$k[[which.min(x$davies_bouldin)]],
      explained_unweighted_kneedle = fastINFLECT::inflect_kneedle(
        x$k,
        x$explained_unweighted
      ),
      explained_event_weighted_kneedle = fastINFLECT::inflect_kneedle(
        x$k,
        x$explained_event_weighted
      )
    )
    for (method in names(recommendations)) {
      index <- index + 1L
      output[[index]] <- data.frame(
        code_space = space,
        method = method,
        k = as.integer(recommendations[[method]]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, output)
}

criterion_selections <- function(criterion_summary) {
  select_internal <- getFromNamespace(".inflect_selection", "fastINFLECT")
  rows <- list()
  row_index <- 0L

  for (criterion in c("dip", "iqr", "combined")) {
    scores <- criterion_summary[
      criterion_summary$criterion == criterion,
      c("k", "qc_pass_rate"),
      drop = FALSE
    ]
    scores$criterion <- criterion
    selected <- tryCatch(
      {
        diagnostic <- fastINFLECT::QC.to.curve(collection.U = scores)
        value <- select_internal(
          collection.U = scores,
          lfunction = diagnostic$lfunction,
          fittedcurve = diagnostic$fittedcurve,
          target = 0.95
        )
        value$error <- NA_character_
        value
      },
      error = function(e) {
        data.frame(
          method = c("inflection", "kneedle", "threshold"),
          k = c(
            NA_real_,
            fastINFLECT::inflect_kneedle(
              scores$k,
              scores$qc_pass_rate
            ),
            fastINFLECT::inflect_threshold_k(scores, target = 0.95)
          ),
          qc_pass_rate_at_k = NA_real_,
          directly_tested = FALSE,
          partition_available = FALSE,
          k_status = "selection_fit_failed",
          score_source = NA_character_,
          target = c(NA_real_, NA_real_, 95),
          unimodality_at_k = NA_real_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
    selected$criterion <- criterion
    row_index <- row_index + 1L
    rows[[row_index]] <- selected
  }
  do.call(rbind, rows)
}

candidate_details <- function(result) {
  fields <- c(
    "dip_p_value",
    "iqr",
    "dip_pass",
    "iqr_pass",
    "combined_pass",
    "criterion_pass",
    "event_count",
    "test_event_count",
    "failure_reason"
  )
  output <- list()
  output_index <- 0L

  for (k in candidate_k) {
    detail <- result$qc.details[[as.character(k)]]
    template <- detail$dip_pass
    row <- expand.grid(
      cluster = as.integer(rownames(template)),
      marker = colnames(template),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    row$k <- k
    for (field in fields) {
      if (!is.null(detail[[field]])) {
        row[[field]] <- as.vector(detail[[field]])
      }
    }
    output_index <- output_index + 1L
    output[[output_index]] <- row
  }
  do.call(rbind, output)
}

refresh_cached_partition_metrics <- function(checkpoint, checkpoint_path) {
  current_metrics <- checkpoint$partition_metrics
  already_current <- identical(
    checkpoint$derived_metric_provenance$version,
    derived_metrics_version
  ) &&
    is.data.frame(current_metrics) &&
    "explained_event_weighted" %in% names(current_metrics) &&
    all(is.finite(current_metrics$explained_event_weighted)) &&
    all(c(
      "event_populated_clusters",
      "empty_event_clusters"
    ) %in% names(current_metrics))
  if (already_current) {
    return(checkpoint)
  }
  required <- c(
    "inflect_codes",
    "historical_codes",
    "node_counts",
    "schedule"
  )
  if (!all(required %in% names(checkpoint))) {
    stop(
      "Cached full checkpoint lacks inputs required to repair derived metrics.",
      call. = FALSE
    )
  }
  input_md5 <- unname(tools::md5sum(checkpoint_path))
  previous_nonfinite <- if (is.data.frame(current_metrics) &&
      "explained_event_weighted" %in% names(current_metrics)) {
    sum(!is.finite(current_metrics$explained_event_weighted))
  } else {
    NA_integer_
  }
  checkpoint$partition_metrics <- rbind(
    partition_metrics(
      checkpoint$inflect_codes,
      checkpoint$node_counts,
      "model_effective_codes",
      ks = checkpoint$schedule
    ),
    partition_metrics(
      checkpoint$historical_codes,
      checkpoint$node_counts,
      "historical_80_20_codes",
      ks = checkpoint$schedule
    )
  )
  checkpoint$metric_recommendations <- metric_recommendations(
    checkpoint$partition_metrics
  )
  checkpoint$derived_metric_provenance <- list(
    version = derived_metrics_version,
    repaired = TRUE,
    repaired_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    input_checkpoint_md5 = input_md5,
    previous_nonfinite_event_weighted = previous_nonfinite,
    rule = paste0(
      "Zero-event SOM nodes remain in nominal partitions and occupancy ",
      "counts but are excluded from event-weighted centroid calculations."
    )
  )
  atomic_save_rds(checkpoint, checkpoint_path)
  checkpoint
}

run_full_inflect <- function() {
  if (file.exists(paths$full_inflect)) {
    checkpoint <- readRDS(paths$full_inflect)
    return(refresh_cached_partition_metrics(checkpoint, paths$full_inflect))
  }
  if (!file.exists(model_path)) {
    stop("Model does not exist: ", model_path, call. = FALSE)
  }
  pkgload::load_all(repo_root, quiet = TRUE)
  load_timing <- system.time(model <- readRDS(model_path))
  view <- getFromNamespace("as_inflect_som", "fastINFLECT")(model)
  stopifnot(
    nrow(view$data) == 39050953L,
    view$map$nNodes == 900L,
    length(view$map$colsUsed) == 27L
  )
  node_counts <- tabulate(
    view$map$mapping[, 1L],
    nbins = view$map$nNodes
  )
  inflect_codes <- as.matrix(view$map$codes)
  historical_codes <- cbind(
    scale_codebook(as.matrix(model$codes[[1L]])) * sqrt(0.8),
    scale_codebook(as.matrix(model$codes[[2L]])) * sqrt(0.2)
  )

  run_warnings <- character()
  run_timing <- system.time({
    result <- withCallingHandlers(
      fastINFLECT::INFLECT(
        FlowSOM.results = model,
        set.i = schedule,
        workers = 2L,
        uniform.test = "both",
        max.n.diptest = NULL,
        max.events.per.node = NULL,
        seed = 42L,
        progress = FALSE
      ),
      warning = function(w) {
        run_warnings <<- c(run_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })
  stopifnot(
    identical(names(result$metaclustering.list), as.character(schedule)),
    all(vapply(
      seq_along(schedule),
      function(i) {
        length(unique(result$metaclustering.list[[i]])) == schedule[[i]]
      },
      logical(1)
    ))
  )

  metrics <- rbind(
    partition_metrics(inflect_codes, node_counts, "model_effective_codes"),
    partition_metrics(historical_codes, node_counts, "historical_80_20_codes")
  )
  checkpoint <- list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    model_path = normalizePath(model_path, mustWork = TRUE),
    n_events = nrow(view$data),
    n_nodes = view$map$nNodes,
    markers = view$prettyColnames[view$map$colsUsed],
    schedule = schedule,
    candidate_k = candidate_k,
    criterion_summary = result$provenance$criterion_summary,
    criterion_selections = criterion_selections(
      result$provenance$criterion_summary
    ),
    candidate_details = candidate_details(result),
    package_selection = result$selection,
    partition_metrics = metrics,
    metric_recommendations = metric_recommendations(metrics),
    candidate_partitions = result$metaclustering.list[
      as.character(candidate_k)
    ],
    inflect_codes = inflect_codes,
    historical_codes = historical_codes,
    node_counts = node_counts,
    model_weights = result$provenance$model,
    runtime = list(
      model_load_seconds = unname(load_timing[["elapsed"]]),
      inflect_seconds = unname(run_timing[["elapsed"]]),
      peak_rss_kb = read_status_kb(),
      warnings = unique(run_warnings)
    ),
    derived_metric_provenance = list(
      version = derived_metrics_version,
      repaired = FALSE,
      generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
      rule = paste0(
        "Zero-event SOM nodes remain in nominal partitions and occupancy ",
        "counts but are excluded from event-weighted centroid calculations."
      )
    ),
    provenance = result$provenance,
    dependency_versions = c(
      R = as.character(getRversion()),
      fastINFLECT = package_version_or_na("fastINFLECT"),
      FlowSOM = package_version_or_na("FlowSOM"),
      ConsensusClusterPlus = package_version_or_na("ConsensusClusterPlus"),
      cluster = package_version_or_na("cluster"),
      diptest = package_version_or_na("diptest")
    )
  )
  rm(model, view, result)
  invisible(gc())
  atomic_save_rds(checkpoint, paths$full_inflect)
  checkpoint
}

consensus_cdf_area <- function(consensus_matrix) {
  values <- consensus_matrix[upper.tri(consensus_matrix)]
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(NA_real_)
  }
  ## Integral from 0 to 1 of the empirical CDF: E(1 - X).
  1 - mean(values)
}

run_consensus_seed <- function(codes, seed) {
  output <- file.path(
    paths$consensus_dir,
    sprintf("consensus-seed-%04d.rds", seed)
  )
  if (file.exists(output)) {
    return(readRDS(output))
  }
  started <- proc.time()[["elapsed"]]
  results <- ConsensusClusterPlus::ConsensusClusterPlus(
    t(codes),
    maxK = max(schedule),
    reps = consensus_reps,
    pItem = 0.9,
    pFeature = 1,
    title = selection_work_dir,
    plot = NULL,
    writeTable = FALSE,
    verbose = FALSE,
    clusterAlg = "hc",
    distance = "euclidean",
    seed = seed
  )
  all_k <- 2L:max(schedule)
  area <- vapply(
    all_k,
    function(k) consensus_cdf_area(results[[k]]$consensusMatrix),
    numeric(1)
  )
  delta <- c(area[[1L]], diff(area) / head(area, -1L))
  classes <- lapply(
    schedule,
    function(k) as.integer(results[[k]]$consensusClass)
  )
  names(classes) <- as.character(schedule)
  checkpoint <- list(
    seed = seed,
    cdf = data.frame(
      seed = seed,
      k = all_k,
      cdf_area = area,
      relative_delta_area = delta,
      stringsAsFactors = FALSE
    ),
    classes = classes,
    wall_seconds = proc.time()[["elapsed"]] - started,
    peak_rss_kb = read_status_kb()
  )
  rm(results)
  invisible(gc())
  atomic_save_rds(checkpoint, output)
  checkpoint
}

summarize_consensus <- function(seed_results) {
  cdf <- do.call(rbind, lapply(seed_results, `[[`, "cdf"))
  rows <- vector("list", length(schedule))
  pair_names <- utils::combn(as.character(seeds), 2L, simplify = FALSE)

  for (i in seq_along(schedule)) {
    k <- schedule[[i]]
    deltas <- cdf$relative_delta_area[cdf$k == k]
    areas <- cdf$cdf_area[cdf$k == k]
    aris <- vapply(
      pair_names,
      function(pair) {
        adjusted_rand_index(
          seed_results[[pair[[1L]]]]$classes[[as.character(k)]],
          seed_results[[pair[[2L]]]]$classes[[as.character(k)]]
        )
      },
      numeric(1)
    )
    rows[[i]] <- data.frame(
      k = k,
      mean_cdf_area = mean(areas),
      median_relative_delta_area = stats::median(deltas),
      min_seed_ari = min(aris),
      median_seed_ari = stats::median(aris),
      stringsAsFactors = FALSE
    )
  }
  summary <- do.call(rbind, rows)
  qualifies <- summary$median_relative_delta_area <= consensus_plateau_delta &
    summary$min_seed_ari >= 0.95
  plateau_k <- NA_integer_
  if (length(qualifies) >= consensus_plateau_width) {
    for (i in seq_len(length(qualifies) - consensus_plateau_width + 1L)) {
      window <- seq.int(i, i + consensus_plateau_width - 1L)
      if (all(qualifies[window])) {
        plateau_k <- summary$k[[i]]
        break
      }
    }
  }
  list(
    per_seed_cdf = cdf,
    summary = summary,
    plateau = data.frame(
      method = "first_5k_consensus_plateau",
      k = plateau_k,
      relative_delta_threshold = consensus_plateau_delta,
      minimum_seed_ari = 0.95,
      window_width = consensus_plateau_width,
      stringsAsFactors = FALSE
    )
  )
}

iqr_checkpoint_path <- function(solution_id, cluster, seed) {
  file.path(
    paths$iqr_dir,
    solution_id,
    sprintf("cluster-%03d-seed-%04d.rds", cluster, seed)
  )
}

run_sample_iqr <- function() {
  manifest_path <- file.path(modality_work_dir, "sample-manifest.rds")
  if (!file.exists(manifest_path)) {
    stop(
      "Modality stage has not produced sample-manifest.rds.",
      call. = FALSE
    )
  }
  manifest <- readRDS(manifest_path)
  for (i in seq_len(nrow(manifest))) {
    task <- manifest[i, , drop = FALSE]
    output <- iqr_checkpoint_path(
      task$solution_id,
      task$cluster,
      task$seed
    )
    if (file.exists(output)) {
      next
    }
    sample_object <- readRDS(task$sample_path)
    values <- sample_object$values
    rows <- lapply(
      seq_len(ncol(values)),
      function(j) {
        finite <- values[, j]
        finite <- finite[is.finite(finite)]
        iqr <- if (length(finite) <= 1L) {
          NA_real_
        } else {
          stats::IQR(finite, type = 7)
        }
        data.frame(
          solution_id = sample_object$solution_id,
          cluster = sample_object$cluster,
          cluster_events = sample_object$cluster_events,
          seed = sample_object$seed,
          marker = colnames(values)[[j]],
          sample_size = length(finite),
          iqr = iqr,
          iqr_pass = if (is.na(iqr)) NA else iqr < iqr_threshold,
          stringsAsFactors = FALSE
        )
      }
    )
    atomic_save_rds(do.call(rbind, rows), output)
  }
  checkpoints <- list.files(
    paths$iqr_dir,
    pattern = "[.]rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  do.call(rbind, lapply(checkpoints, readRDS))
}

sample_iqr_summary <- function(rows) {
  solutions <- unique(rows$solution_id)
  output <- vector("list", length(solutions))
  for (i in seq_along(solutions)) {
    solution <- rows[rows$solution_id == solutions[[i]], , drop = FALSE]
    valid <- !is.na(solution$iqr_pass)
    weights <- solution$cluster_events
    output[[i]] <- data.frame(
      solution_id = solutions[[i]],
      valid_cells = sum(valid),
      unresolved_cells = sum(!valid),
      unweighted_iqr_pass_rate = if (sum(valid) == 0L) {
        NA_real_
      } else {
        mean(solution$iqr_pass[valid]) * 100
      },
      event_weighted_iqr_pass_rate = if (sum(valid) == 0L) {
        NA_real_
      } else {
        stats::weighted.mean(
          solution$iqr_pass[valid],
          weights[valid]
        ) * 100
      },
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}

source_provenance <- function() {
  source_files <- selection_source_files()
  list(
    git_commit = system2(
      "git",
      c("-C", repo_root, "rev-parse", "HEAD"),
      stdout = TRUE
    )[[1L]],
    git_status = system2(
      "git",
      c("-C", repo_root, "status", "--short"),
      stdout = TRUE
    ),
    source_md5 = unname(tools::md5sum(source_files)),
    source_files = source_files
  )
}

run_analysis <- function() {
  started <- proc.time()[["elapsed"]]
  stage_path <- file.path(modality_work_dir, "stage.rds")
  if (!file.exists(stage_path)) {
    stop("Modality stage has not produced stage.rds.", call. = FALSE)
  }
  stage <- readRDS(stage_path)
  contract <- ensure_contract(stage)
  full <- run_full_inflect()
  seed_results <- lapply(
    seeds,
    function(seed) run_consensus_seed(full$inflect_codes, seed)
  )
  names(seed_results) <- as.character(seeds)
  consensus <- summarize_consensus(seed_results)
  sample_iqr <- run_sample_iqr()

  evidence <- list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    model_path = full$model_path,
    model_md5 = stage$model_md5,
    n_events = full$n_events,
    n_nodes = full$n_nodes,
    schedule = schedule,
    candidate_k = candidate_k,
    contract = contract,
    full_inflect = full[
      setdiff(names(full), c("inflect_codes", "historical_codes"))
    ],
    consensus = consensus,
    consensus_seed_runtime = do.call(
      rbind,
      lapply(seed_results, function(x) {
        data.frame(
          seed = x$seed,
          wall_seconds = x$wall_seconds,
          peak_rss_kb = x$peak_rss_kb,
          stringsAsFactors = FALSE
        )
      })
    ),
    consensus_seed_provenance = lapply(seed_results, function(x) {
      if (!is.null(x$producer)) {
        return(x$producer)
      }
      list(
        implementation = "serial main selection workflow",
        script_file = normalizePath(script_file[[1L]], mustWork = TRUE),
        script_md5 = unname(tools::md5sum(script_file[[1L]]))
      )
    }),
    sample_iqr = sample_iqr,
    sample_iqr_summary = sample_iqr_summary(sample_iqr),
    parameters = list(
      consensus_reps = consensus_reps,
      consensus_p_item = 0.9,
      consensus_p_feature = 1,
      consensus_plateau_delta = consensus_plateau_delta,
      consensus_plateau_width = consensus_plateau_width,
      iqr_threshold = iqr_threshold,
      sample_iqr_events = 5000L,
      sample_iqr_value_handling = "require finite and retain all transformed values"
    ),
    source_provenance = source_provenance(),
    wall_seconds = proc.time()[["elapsed"]] - started,
    peak_rss_kb = read_status_kb(),
    session_info = utils::sessionInfo()
  )
  atomic_save_rds(evidence, paths$final)
  atomic_write_csv(
    evidence$full_inflect$criterion_summary,
    file.path(evidence_dir, "criterion-curves.csv")
  )
  atomic_write_csv(
    evidence$full_inflect$criterion_selections,
    file.path(evidence_dir, "criterion-selections.csv")
  )
  atomic_write_csv(
    evidence$full_inflect$partition_metrics,
    file.path(evidence_dir, "partition-metrics.csv")
  )
  atomic_write_csv(
    evidence$full_inflect$metric_recommendations,
    file.path(evidence_dir, "metric-recommendations.csv")
  )
  atomic_write_csv(
    evidence$consensus$summary,
    file.path(evidence_dir, "consensus-stability.csv")
  )
  atomic_write_csv(
    evidence$sample_iqr_summary,
    file.path(evidence_dir, "sample-iqr-summary.csv")
  )
  message("Selection evidence complete: ", paths$final)
  invisible(evidence)
}

run_self_test <- function() {
  stopifnot(
    identical(adjusted_rand_index(1:5, 1:5), 1),
    isTRUE(all.equal(
      consensus_cdf_area(matrix(c(1, 0, 0, 1), 2L)),
      1
    ))
  )
  codes <- rbind(
    c(0, 0),
    c(0, 0.2),
    c(5, 5),
    c(5.2, 5)
  )
  metrics <- partition_metrics(codes, rep(1, 4L), "test", ks = 2L:3L)
  zero_weight <- weighted_dispersion(
    rbind(c(0, 0), c(2, 2), c(100, 100)),
    c(1L, 2L, 3L),
    c(1, 1, 0)
  )
  stopifnot(
    nrow(metrics) == 2L,
    isTRUE(all.equal(zero_weight$explained, 1)),
    all(is.finite(unlist(zero_weight[c("within", "total", "explained")]))),
    all(c(
      "silhouette",
      "calinski_harabasz",
      "davies_bouldin"
    ) %in% names(metrics))
  )
  message("Selection workflow self-test passed.")
}

switch(
  mode,
  run = run_analysis(),
  status = {
    print(list(
      full_inflect = file.exists(paths$full_inflect),
      consensus_seeds = list.files(
        paths$consensus_dir,
        pattern = "[.]rds$"
      ),
      iqr_checkpoints = length(list.files(
        paths$iqr_dir,
        pattern = "[.]rds$",
        recursive = TRUE
      )),
      final = file.exists(paths$final)
    ))
  },
  `self-test` = run_self_test(),
  stop("Unknown mode: ", mode, call. = FALSE)
)
