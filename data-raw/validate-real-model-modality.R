#!/usr/bin/env Rscript

## Checkpointed real-model modality audit for fastINFLECT 2.0.
##
## The workflow is intentionally staged because the source model contains
## 39,050,953 events and ACR mode tests use B = 1999:
##
##   Rscript data-raw/validate-real-model-modality.R --mode=stage
##   Rscript data-raw/validate-real-model-modality.R --mode=pilot
##   INFLECT_ARRAY_COUNT=64 SLURM_ARRAY_TASK_ID=0 \
##     Rscript data-raw/validate-real-model-modality.R --mode=base
##   Rscript data-raw/validate-real-model-modality.R --mode=classify
##   INFLECT_ARRAY_COUNT=64 SLURM_ARRAY_TASK_ID=0 \
##     Rscript data-raw/validate-real-model-modality.R --mode=refine
##   Rscript data-raw/validate-real-model-modality.R --mode=summarize
##
## Array indices are zero-based. Repeat base/refine for every index in
## 0:(INFLECT_ARRAY_COUNT - 1). Existing cluster checkpoints are reused.

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0L) {
    return(default)
  }
  sub(prefix, "", hit[[length(hit)]], fixed = TRUE)
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[startsWith(script_arg, "--file=")])
repo_root <- if (length(script_file) > 0L) {
  normalizePath(file.path(dirname(script_file[[1L]]), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

mode <- arg_value("mode", "status")
model_path <- Sys.getenv(
  "INFLECT_LARGE_MODEL",
  unset = paste0(
    "/exports/para-lipg-hpc/Xuran/results/bmv/aurora_exvivo/cluster/",
    "model/xyf_som_model_with_LD_filter.rds"
  )
)
work_dir <- Sys.getenv(
  "INFLECT_MODALITY_WORKDIR",
  unset = file.path(repo_root, "data-raw", ".real-model-modality-work")
)
evidence_dir <- Sys.getenv(
  "INFLECT_MODALITY_EVIDENCE_DIR",
  unset = file.path(repo_root, "inst", "benchmarks", "real-model-modality")
)
local_r_lib <- Sys.getenv(
  "INFLECT_MODALITY_R_LIB",
  unset = file.path(repo_root, "data-raw", ".real-model-lib")
)
if (dir.exists(local_r_lib)) {
  .libPaths(c(normalizePath(local_r_lib, mustWork = TRUE), .libPaths()))
}
shard_count <- as.integer(Sys.getenv("INFLECT_ARRAY_COUNT", unset = "1"))
shard_index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
acr_B <- 1999L
base_sample_size <- 2000L
refinement_sample_sizes <- c(1000L, 5000L)
max_staged_sample <- max(c(base_sample_size, refinement_sample_sizes))
seeds <- c(1L, 42L, 2026L)
q_threshold <- 0.05
borderline_q <- 0.10

if (is.na(shard_count) || shard_count < 1L) {
  stop("INFLECT_ARRAY_COUNT must be a positive integer.", call. = FALSE)
}
if (is.na(shard_index) || shard_index < 0L || shard_index >= shard_count) {
  stop(
    "SLURM_ARRAY_TASK_ID must be in [0, INFLECT_ARRAY_COUNT).",
    call. = FALSE
  )
}

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(evidence_dir, recursive = TRUE, showWarnings = FALSE)
dirs <- list(
  samples = file.path(work_dir, "samples"),
  base = file.path(work_dir, "base"),
  refine = file.path(work_dir, "refine"),
  logs = file.path(work_dir, "logs")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

paths <- list(
  stage = file.path(work_dir, "stage.rds"),
  sample_manifest = file.path(work_dir, "sample-manifest.rds"),
  pilot = file.path(work_dir, "pilot-100-acr.rds"),
  classifications = file.path(work_dir, "classifications.rds"),
  refinement_manifest = file.path(work_dir, "refinement-manifest.rds"),
  final_bundle = file.path(evidence_dir, "modality-evidence.rds")
)

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

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  max(x)
}

source_provenance <- function() {
  source_files <- c(
    normalizePath(script_file[[1L]], mustWork = TRUE),
    list.files(
      file.path(repo_root, "R"),
      pattern = "[.]R$",
      full.names = TRUE
    ),
    file.path(repo_root, "DESCRIPTION")
  )
  git_output <- function(...) {
    tryCatch(
      system2(
        "git",
        c("-C", repo_root, ...),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) NA_character_
    )
  }
  list(
    git_commit = git_output("rev-parse", "HEAD")[[1L]],
    git_status = git_output("status", "--short"),
    source_md5 = unname(tools::md5sum(source_files)),
    source_files = normalizePath(source_files, mustWork = TRUE)
  )
}

empty_modality_rows <- function() {
  data.frame(
    solution_id = character(),
    cluster = integer(),
    cluster_events = integer(),
    seed = integer(),
    marker = character(),
    requested_sample_size = integer(),
    sample_size = integer(),
    unique_values = integer(),
    unique_fraction = numeric(),
    max_tie_fraction = numeric(),
    boundary_fraction = numeric(),
    discrete_reason = character(),
    dip_p = numeric(),
    dip_statistic = numeric(),
    acr_p = numeric(),
    acr_statistic = numeric(),
    failure_reason = character(),
    dip_warnings = character(),
    acr_warnings = character(),
    dip_seconds = numeric(),
    acr_seconds = numeric(),
    dip_q = numeric(),
    acr_q = numeric(),
    bh_family = character(),
    stringsAsFactors = FALSE
  )
}

dependency_versions <- function() {
  packages <- c(
    "fastINFLECT",
    "diptest",
    "multimode",
    "FlowSOM",
    "kohonen",
    "ggplot2",
    "R"
  )
  stats::setNames(
    vapply(packages, function(package) {
      if (package == "R") {
        return(as.character(getRversion()))
      }
      if (package == "fastINFLECT") {
        return(unname(
          read.dcf(file.path(repo_root, "DESCRIPTION"))[1L, "Version"]
        ))
      }
      if (!requireNamespace(package, quietly = TRUE)) {
        return(NA_character_)
      }
      as.character(utils::packageVersion(package))
    }, character(1)),
    packages
  )
}

safe_seed <- function(seed, ...) {
  modulus <- as.double(.Machine$integer.max)
  key <- enc2utf8(paste(c(...), collapse = "\u001f"))
  bytes <- as.integer(charToRaw(key))
  value <- as.double(seed) %% modulus
  if (length(bytes) > 0L) {
    for (byte in bytes) {
      value <- (value * 131 + byte + 1) %% modulus
    }
  }
  as.integer(value)
}

with_preserved_seed <- function(seed, code) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

adjusted_rand_index <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2L) {
    return(NA_real_)
  }
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(n)
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- (sum_rows + sum_cols) / 2
  denominator <- maximum - expected
  if (denominator == 0) {
    return(if (identical(as.integer(x), as.integer(y))) 1 else NA_real_)
  }
  (sum_cells - expected) / denominator
}

validate_partition <- function(labels, expected_k, n_nodes, solution_id) {
  if (length(labels) != n_nodes || anyNA(labels)) {
    stop(
      solution_id,
      " does not label every SOM node exactly once.",
      call. = FALSE
    )
  }
  unique_labels <- sort(unique(as.integer(labels)))
  if (length(unique_labels) != expected_k) {
    stop(
      solution_id,
      " has ",
      length(unique_labels),
      " clusters; expected ",
      expected_k,
      ".",
      call. = FALSE
    )
  }
  as.integer(match(labels, unique_labels))
}

scale_codebook <- function(codes) {
  scaled <- scale(codes)
  scaled[is.na(scaled)] <- 0
  scaled
}

build_partitions <- function(model, view) {
  inflect_codes <- view$map$codes
  inflect_fit <- stats::hclust(
    stats::dist(inflect_codes),
    method = "ward.D2"
  )
  partitions <- list()
  solution_metadata <- list()

  for (k in c(55L, 60L, 62L, 65L)) {
    solution_id <- paste0("inflect_ward_k", k)
    partitions[[solution_id]] <- validate_partition(
      stats::cutree(inflect_fit, k = k),
      k,
      view$map$nNodes,
      solution_id
    )
    solution_metadata[[solution_id]] <- data.frame(
      solution_id = solution_id,
      family = "INFLECT Ward.D2",
      k = k,
      code_space = "model effective code matrix",
      seed = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  if (is.null(model$codes) || length(model$codes) < 2L) {
    stop(
      "Historical 80/20 reconstruction requires X and Y kohonen codebooks.",
      call. = FALSE
    )
  }
  x_codes <- scale_codebook(as.matrix(model$codes[[1L]]))
  y_codes <- scale_codebook(as.matrix(model$codes[[2L]]))
  historical_codes <- cbind(
    x_codes * sqrt(0.8),
    y_codes * sqrt(0.2)
  )
  historical_fit <- stats::hclust(
    stats::dist(historical_codes),
    method = "ward.D2"
  )
  for (k in c(50L, 62L)) {
    solution_id <- paste0("historical_80_20_ward_k", k)
    partitions[[solution_id]] <- validate_partition(
      stats::cutree(historical_fit, k = k),
      k,
      view$map$nNodes,
      solution_id
    )
    solution_metadata[[solution_id]] <- data.frame(
      solution_id = solution_id,
      family = "historical 80/20 Ward.D2",
      k = k,
      code_space = "standardized X/Y codebooks; 0.8/0.2 squared-distance weights",
      seed = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  consensus_ids <- character(0)
  for (seed in seeds) {
    solution_id <- paste0("flowsom_consensus_k62_seed", seed)
    consensus_ids <- c(consensus_ids, solution_id)
    consensus <- FlowSOM::metaClustering_consensus(
      inflect_codes,
      k = 62L,
      seed = seed
    )
    partitions[[solution_id]] <- validate_partition(
      consensus,
      62L,
      view$map$nNodes,
      solution_id
    )
    solution_metadata[[solution_id]] <- data.frame(
      solution_id = solution_id,
      family = "FlowSOM consensus",
      k = 62L,
      code_space = "model effective code matrix",
      seed = seed,
      stringsAsFactors = FALSE
    )
  }

  consensus_ari <- data.frame()
  pair_index <- 0L
  for (i in seq_len(length(consensus_ids) - 1L)) {
    for (j in seq.int(i + 1L, length(consensus_ids))) {
      pair_index <- pair_index + 1L
      consensus_ari <- rbind(
        consensus_ari,
        data.frame(
          solution_a = consensus_ids[[i]],
          solution_b = consensus_ids[[j]],
          ari = adjusted_rand_index(
            partitions[[consensus_ids[[i]]]],
            partitions[[consensus_ids[[j]]]]
          ),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  stable_consensus <- all(consensus_ari$ari >= 0.95)
  display_controls <- if (stable_consensus) {
    "flowsom_consensus_k62_seed42"
  } else {
    consensus_ids
  }

  list(
    partitions = partitions,
    solution_metadata = do.call(rbind, solution_metadata),
    inflect_codes = inflect_codes,
    historical_codes = historical_codes,
    consensus_ari = consensus_ari,
    consensus_stable = stable_consensus,
    display_controls = display_controls
  )
}

sample_path <- function(solution_id, cluster, seed) {
  file.path(
    dirs$samples,
    solution_id,
    sprintf("cluster-%03d-seed-%04d.rds", cluster, seed)
  )
}

base_path <- function(solution_id, cluster, seed) {
  file.path(
    dirs$base,
    solution_id,
    sprintf("cluster-%03d-seed-%04d.rds", cluster, seed)
  )
}

refine_path <- function(solution_id, cluster, seed, sample_size) {
  file.path(
    dirs$refine,
    solution_id,
    sprintf(
      "cluster-%03d-seed-%04d-n-%04d.rds",
      cluster,
      seed,
      sample_size
    )
  )
}

run_stage <- function() {
  if (!file.exists(model_path)) {
    stop("Model does not exist: ", model_path, call. = FALSE)
  }
  pkgload::load_all(repo_root, quiet = TRUE)
  stage_start <- proc.time()[["elapsed"]]
  message("Loading model: ", model_path)
  model <- readRDS(model_path)
  view <- as_inflect_som(model)
  if (nrow(view$data) != 39050953L) {
    stop(
      "Model event count changed: expected 39,050,953; found ",
      nrow(view$data),
      ".",
      call. = FALSE
    )
  }
  if (ncol(view$data) != 27L) {
    stop("Expected 27 marker columns; found ", ncol(view$data), ".", call. = FALSE)
  }

  partition_data <- build_partitions(model, view)
  node_counts <- tabulate(
    view$map$mapping[, 1L],
    nbins = view$map$nNodes
  )
  node_events <- split(
    seq_len(nrow(view$data)),
    factor(view$map$mapping[, 1L], levels = seq_len(view$map$nNodes))
  )

  manifest_rows <- list()
  manifest_index <- 0L
  for (solution_id in names(partition_data$partitions)) {
    partition <- partition_data$partitions[[solution_id]]
    for (cluster in seq_len(max(partition))) {
      cluster_nodes <- which(partition == cluster)
      cluster_events <- sum(node_counts[cluster_nodes])
      for (seed in seeds) {
        manifest_index <- manifest_index + 1L
        manifest_rows[[manifest_index]] <- data.frame(
          solution_id = solution_id,
          cluster = cluster,
          seed = seed,
          cluster_events = cluster_events,
          staged_events = min(cluster_events, max_staged_sample),
          sample_path = sample_path(solution_id, cluster, seed),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  sample_manifest <- do.call(rbind, manifest_rows)
  atomic_save_rds(sample_manifest, paths$sample_manifest)

  for (solution_id in names(partition_data$partitions)) {
    partition <- partition_data$partitions[[solution_id]]
    for (cluster in seq_len(max(partition))) {
      expected_paths <- vapply(
        seeds,
        function(seed) sample_path(solution_id, cluster, seed),
        character(1)
      )
      if (all(file.exists(expected_paths))) {
        next
      }
      cluster_nodes <- which(partition == cluster)
      cluster_rows <- unlist(node_events[cluster_nodes], use.names = FALSE)
      for (seed in seeds) {
        output <- sample_path(solution_id, cluster, seed)
        if (file.exists(output)) {
          next
        }
        take <- min(length(cluster_rows), max_staged_sample)
        sampled_rows <- if (take < length(cluster_rows)) {
          with_preserved_seed(
            safe_seed(0L, solution_id, cluster, seed, "event_sample"),
            sample(cluster_rows, size = take, replace = FALSE)
          )
        } else {
          cluster_rows
        }
        sampled_matrix <- view$data[
          sampled_rows,
          view$map$colsUsed,
          drop = FALSE
        ]
        colnames(sampled_matrix) <- view$prettyColnames[view$map$colsUsed]
        atomic_save_rds(
          list(
            solution_id = solution_id,
            cluster = cluster,
            seed = seed,
            cluster_events = length(cluster_rows),
            row_ids = sampled_rows,
            values = sampled_matrix
          ),
          output
        )
      }
      rm(cluster_rows)
    }
  }

  source <- view$inflect_source
  stage <- list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    model_path = normalizePath(model_path, mustWork = TRUE),
    model_md5 = unname(tools::md5sum(model_path)),
    model_size_bytes = file.info(model_path)$size,
    model_mtime = format(file.info(model_path)$mtime, "%Y-%m-%dT%H:%M:%OS%z"),
    n_events = nrow(view$data),
    n_nodes = view$map$nNodes,
    markers = view$prettyColnames[view$map$colsUsed],
    model_weights = list(
      user_weights = source$user_weights,
      distance_weights = source$distance_weights,
      effective_layer_weights = source$effective_layer_weights
    ),
    historical_weighting = list(
      x = 0.8,
      y = 0.2,
      source_notebook = paste0(
        "/exports/para-lipg-hpc/Xuran/scripts/bmv/1_cluster/",
        "01_xyf_SOM_exvivo.ipynb"
      ),
      reconstruction = paste0(
        "scale X and Y codebooks separately; replace scale-generated NA with 0; ",
        "multiply by sqrt(0.8) and sqrt(0.2); Ward.D2"
      )
    ),
    solution_metadata = partition_data$solution_metadata,
    partitions = partition_data$partitions,
    consensus_ari = partition_data$consensus_ari,
    consensus_stable = partition_data$consensus_stable,
    display_controls = partition_data$display_controls,
    sample_sizes = list(
      base = base_sample_size,
      refinement = refinement_sample_sizes,
      staged_max = max_staged_sample
    ),
    seeds = seeds,
    sample_manifest = sample_manifest,
    source_provenance = source_provenance(),
    dependency_versions = dependency_versions(),
    library_paths = .libPaths(),
    stage_wall_seconds = proc.time()[["elapsed"]] - stage_start,
    stage_peak_rss_kb = read_status_kb(),
    session_info = utils::sessionInfo()
  )
  rm(node_events, node_counts, partition_data, model, view)
  invisible(gc())
  atomic_save_rds(stage, paths$stage)
  message("Stage complete: ", paths$stage)
}

capture_condition <- function(code) {
  warnings <- character(0)
  value <- withCallingHandlers(
    tryCatch(
      list(value = force(code), error = NA_character_),
      error = function(e) {
        list(value = NULL, error = conditionMessage(e))
      }
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  value$warnings <- unique(warnings)
  value
}

discrete_diagnostics <- function(values) {
  finite <- values[is.finite(values)]
  n <- length(finite)
  if (n == 0L) {
    return(list(
      finite = finite,
      n = 0L,
      unique_n = 0L,
      unique_fraction = 0,
      max_tie_fraction = NA_real_,
      boundary_fraction = NA_real_,
      reason = "no_finite_values"
    ))
  }
  counts <- table(finite, useNA = "no")
  unique_n <- length(counts)
  max_tie_fraction <- max(counts) / n
  boundary_fraction <- max(
    mean(finite == min(finite)),
    mean(finite == max(finite))
  )
  reasons <- character(0)
  if (n < 4L) {
    reasons <- c(reasons, "fewer_than_4_values")
  }
  if (unique_n < 8L || unique_n / n < 0.01) {
    reasons <- c(reasons, "insufficient_unique_values")
  }
  if (max_tie_fraction >= 0.50) {
    reasons <- c(reasons, "excessive_ties")
  }
  if (boundary_fraction >= 0.20) {
    reasons <- c(reasons, "boundary_saturation")
  }
  list(
    finite = finite,
    n = n,
    unique_n = unique_n,
    unique_fraction = unique_n / n,
    max_tie_fraction = max_tie_fraction,
    boundary_fraction = boundary_fraction,
    reason = if (length(reasons) == 0L) {
      NA_character_
    } else {
      paste(unique(reasons), collapse = ";")
    }
  )
}

run_marker_tests <- function(sample_object,
                             sample_size,
                             markers = colnames(sample_object$values)) {
  if (!requireNamespace("multimode", quietly = TRUE)) {
    stop(
      "Package `multimode` is required for the ACR modality audit.",
      call. = FALSE
    )
  }
  available <- nrow(sample_object$values)
  use_n <- min(sample_size, available)
  values_matrix <- sample_object$values[seq_len(use_n), markers, drop = FALSE]
  output <- vector("list", length(markers))

  for (j in seq_along(markers)) {
    marker <- markers[[j]]
    diagnostics <- discrete_diagnostics(values_matrix[, marker])
    dip_result <- list(value = NULL, error = NA_character_, warnings = character(0))
    acr_result <- list(value = NULL, error = NA_character_, warnings = character(0))
    dip_seconds <- 0
    acr_seconds <- 0

    if (is.na(diagnostics$reason)) {
      dip_timing <- system.time({
        dip_result <- capture_condition(
          diptest::dip.test(diagnostics$finite)
        )
      })
      dip_seconds <- unname(dip_timing[["elapsed"]])
      acr_timing <- system.time({
        acr_result <- capture_condition(
          with_preserved_seed(
            safe_seed(
              0L,
              sample_object$solution_id,
              sample_object$cluster,
              sample_object$seed,
              marker,
              sample_size,
              "acr"
            ),
            multimode::modetest(
              diagnostics$finite,
              mod0 = 1,
              method = "ACR",
              B = acr_B
            )
          )
        )
      })
      acr_seconds <- unname(acr_timing[["elapsed"]])
    }

    dip_value <- dip_result$value
    acr_value <- acr_result$value
    failure_reasons <- c(
      if (!is.na(diagnostics$reason)) diagnostics$reason,
      if (!is.na(dip_result$error)) paste0("dip_error:", dip_result$error),
      if (!is.na(acr_result$error)) paste0("acr_error:", acr_result$error)
    )
    output[[j]] <- data.frame(
      solution_id = sample_object$solution_id,
      cluster = sample_object$cluster,
      cluster_events = sample_object$cluster_events,
      seed = sample_object$seed,
      marker = marker,
      requested_sample_size = sample_size,
      sample_size = diagnostics$n,
      unique_values = diagnostics$unique_n,
      unique_fraction = diagnostics$unique_fraction,
      max_tie_fraction = diagnostics$max_tie_fraction,
      boundary_fraction = diagnostics$boundary_fraction,
      discrete_reason = diagnostics$reason,
      dip_p = if (is.null(dip_value)) NA_real_ else dip_value$p.value,
      dip_statistic = if (is.null(dip_value)) {
        NA_real_
      } else {
        unname(dip_value$statistic)
      },
      acr_p = if (is.null(acr_value)) NA_real_ else acr_value$p.value,
      acr_statistic = if (is.null(acr_value)) {
        NA_real_
      } else {
        unname(acr_value$statistic)
      },
      failure_reason = if (length(failure_reasons) == 0L) {
        NA_character_
      } else {
        paste(failure_reasons, collapse = ";")
      },
      dip_warnings = paste(dip_result$warnings, collapse = " | "),
      acr_warnings = paste(acr_result$warnings, collapse = " | "),
      dip_seconds = dip_seconds,
      acr_seconds = acr_seconds,
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, output)
  result$dip_q <- stats::p.adjust(result$dip_p, method = "BH")
  result$acr_q <- stats::p.adjust(result$acr_p, method = "BH")
  result$bh_family <- paste0(
    ncol(sample_object$values),
    " markers within solution-cluster-method-seed-sample_size"
  )
  result
}

task_rows_for_shard <- function(manifest) {
  task_index <- seq_len(nrow(manifest)) - 1L
  manifest[task_index %% shard_count == shard_index, , drop = FALSE]
}

run_base <- function() {
  if (!file.exists(paths$sample_manifest)) {
    stop("Run --mode=stage first.", call. = FALSE)
  }
  manifest <- task_rows_for_shard(readRDS(paths$sample_manifest))
  message(
    "Base shard ",
    shard_index,
    "/",
    shard_count,
    ": ",
    nrow(manifest),
    " cluster-seed tasks"
  )
  for (i in seq_len(nrow(manifest))) {
    task <- manifest[i, , drop = FALSE]
    output <- base_path(task$solution_id, task$cluster, task$seed)
    if (file.exists(output)) {
      next
    }
    sample_object <- readRDS(task$sample_path)
    timing <- system.time({
      results <- run_marker_tests(sample_object, base_sample_size)
    })
    atomic_save_rds(
      list(
        task = task,
        results = results,
        wall_seconds = unname(timing[["elapsed"]]),
        peak_rss_kb = read_status_kb(),
        acr_B = acr_B
      ),
      output
    )
  }
}

run_pilot <- function() {
  if (!file.exists(paths$sample_manifest)) {
    stop("Run --mode=stage first.", call. = FALSE)
  }
  manifest <- readRDS(paths$sample_manifest)
  remaining <- 100L
  rows <- list()
  row_index <- 0L
  start <- proc.time()[["elapsed"]]
  for (i in seq_len(nrow(manifest))) {
    if (remaining == 0L) {
      break
    }
    sample_object <- readRDS(manifest$sample_path[[i]])
    selected_markers <- head(colnames(sample_object$values), remaining)
    result <- run_marker_tests(
      sample_object,
      base_sample_size,
      markers = selected_markers
    )
    row_index <- row_index + 1L
    rows[[row_index]] <- result
    remaining <- remaining - nrow(result)
  }
  pilot_rows <- do.call(rbind, rows)
  elapsed <- proc.time()[["elapsed"]] - start
  if (nrow(pilot_rows) != 100L) {
    stop("Pilot did not produce exactly 100 ACR tests.", call. = FALSE)
  }
  total_expected <- sum(
    readRDS(paths$stage)$solution_metadata$k
  ) * length(seeds) * length(readRDS(paths$stage)$markers)
  pilot <- list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    n_tests = nrow(pilot_rows),
    acr_B = acr_B,
    wall_seconds = elapsed,
    seconds_per_test = elapsed / nrow(pilot_rows),
    projected_base_test_seconds = elapsed / nrow(pilot_rows) * total_expected,
    peak_rss_kb = read_status_kb(),
    results = pilot_rows,
    session_info = utils::sessionInfo()
  )
  atomic_save_rds(pilot, paths$pilot)
  atomic_write_csv(
    pilot_rows,
    file.path(evidence_dir, "pilot-100-acr-tests.csv")
  )
  message(
    "Pilot: ",
    round(pilot$seconds_per_test, 3),
    " seconds per test; projected serial base wall time ",
    round(pilot$projected_base_test_seconds / 3600, 2),
    " hours."
  )
}

expected_checkpoint_paths <- function(manifest, kind = c("base", "refine")) {
  kind <- match.arg(kind)
  if (nrow(manifest) == 0L) {
    return(character())
  }
  if (kind == "base") {
    return(mapply(
      base_path,
      manifest$solution_id,
      manifest$cluster,
      manifest$seed,
      USE.NAMES = FALSE
    ))
  }
  mapply(
    refine_path,
    manifest$solution_id,
    manifest$cluster,
    manifest$seed,
    manifest$requested_sample_size,
    USE.NAMES = FALSE
  )
}

read_complete_checkpoints <- function(manifest, kind = c("base", "refine")) {
  kind <- match.arg(kind)
  checkpoint_paths <- expected_checkpoint_paths(manifest, kind)
  if (length(checkpoint_paths) == 0L) {
    return(list(
      rows = empty_modality_rows(),
      wall_seconds = 0,
      peak_rss_kb = NA_real_
    ))
  }
  missing <- checkpoint_paths[!file.exists(checkpoint_paths)]
  if (length(missing) > 0L) {
    stop(
      length(missing),
      " ",
      kind,
      " checkpoint(s) are missing; first missing: ",
      missing[[1L]],
      call. = FALSE
    )
  }
  checkpoints <- lapply(checkpoint_paths, readRDS)
  list(
    rows = do.call(rbind, lapply(checkpoints, `[[`, "results")),
    wall_seconds = sum(vapply(checkpoints, `[[`, numeric(1), "wall_seconds")),
    peak_rss_kb = safe_max(
      vapply(checkpoints, `[[`, numeric(1), "peak_rss_kb")
    )
  )
}

classify_pair <- function(rows) {
  rows <- rows[match(seeds, rows$seed), , drop = FALSE]
  if (nrow(rows) != length(seeds) || any(rows$seed != seeds)) {
    stop("Classification requires seeds 1, 42, and 2026.", call. = FALSE)
  }
  valid <- is.na(rows$failure_reason) &
    is.finite(rows$dip_q) &
    is.finite(rows$acr_q)
  dip_reject <- valid & rows$dip_q < q_threshold
  acr_reject <- valid & rows$acr_q < q_threshold
  both_reject <- dip_reject & acr_reject
  borderline <- valid & (
    (rows$dip_q >= q_threshold & rows$dip_q < borderline_q) |
      (rows$acr_q >= q_threshold & rows$acr_q < borderline_q)
  )

  if (!all(valid)) {
    status <- "unresolved_discrete"
    reason <- paste(unique(stats::na.omit(rows$failure_reason)), collapse = " | ")
    if (!nzchar(reason)) {
      reason <- "invalid_or_nonfinite_test_result"
    }
  } else if (sum(both_reject) >= 2L) {
    status <- "detected_multimodality"
    reason <- "both tests reject in at least two seeds"
  } else if (!any(dip_reject) && !any(acr_reject) && !any(borderline)) {
    status <- "no_detected_multimodality"
    reason <- "neither test rejects in any seed"
  } else {
    status <- "ambiguous"
    reasons <- character(0)
    if (any(xor(dip_reject, acr_reject))) {
      reasons <- c(reasons, "discordant_tests")
    }
    if (length(unique(dip_reject)) > 1L ||
        length(unique(acr_reject)) > 1L) {
      reasons <- c(reasons, "seed_instability")
    }
    if (any(borderline)) {
      reasons <- c(reasons, "borderline_q")
    }
    if (length(reasons) == 0L) {
      reasons <- "rejection_not_replicated_by_both_tests"
    }
    reason <- paste(unique(reasons), collapse = ";")
  }
  data.frame(
    solution_id = rows$solution_id[[1L]],
    cluster = rows$cluster[[1L]],
    cluster_events = rows$cluster_events[[1L]],
    marker = rows$marker[[1L]],
    status = status,
    reason = reason,
    valid_seeds = sum(valid),
    both_reject_seeds = sum(both_reject),
    dip_reject_seeds = sum(dip_reject),
    acr_reject_seeds = sum(acr_reject),
    min_dip_q = if (all(is.na(rows$dip_q))) NA_real_ else min(rows$dip_q, na.rm = TRUE),
    min_acr_q = if (all(is.na(rows$acr_q))) NA_real_ else min(rows$acr_q, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

size_strata <- function(values) {
  rank_value <- rank(values, ties.method = "average")
  pmin(4L, pmax(1L, ceiling(rank_value / length(values) * 4L)))
}

build_refinement_manifest <- function(classifications, sample_manifest) {
  classifications$size_stratum <- ave(
    classifications$cluster_events,
    classifications$solution_id,
    FUN = size_strata
  )
  classifications$refinement_reason <- NA_character_
  classifications$refinement_reason[
    classifications$status %in% c("detected_multimodality", "ambiguous")
  ] <- classifications$status[
    classifications$status %in% c("detected_multimodality", "ambiguous")
  ]

  pass_indices <- which(
    classifications$status == "no_detected_multimodality"
  )
  grouping <- interaction(
    classifications$solution_id[pass_indices],
    classifications$marker[pass_indices],
    classifications$size_stratum[pass_indices],
    drop = TRUE
  )
  for (indices in split(pass_indices, grouping)) {
    hashes <- vapply(indices, function(idx) {
      safe_seed(
        0L,
        classifications$solution_id[[idx]],
        classifications$cluster[[idx]],
        classifications$marker[[idx]],
        "pass_audit"
      )
    }, integer(1))
    take <- max(1L, ceiling(0.05 * length(indices)))
    selected <- indices[order(hashes, classifications$cluster[indices])][
      seq_len(take)
    ]
    classifications$refinement_reason[selected] <- "stratified_5_percent_pass_audit"
  }

  selected <- classifications[!is.na(classifications$refinement_reason), ]
  selected_clusters <- unique(selected[, c("solution_id", "cluster")])
  empty_tasks <- data.frame(
    solution_id = character(),
    cluster = integer(),
    seed = integer(),
    requested_sample_size = integer(),
    sample_path = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(selected_clusters) == 0L) {
    return(list(
      classifications = classifications,
      selected_pairs = selected,
      tasks = empty_tasks
    ))
  }
  task_rows <- list()
  index <- 0L
  for (i in seq_len(nrow(selected_clusters))) {
    solution_id <- selected_clusters$solution_id[[i]]
    cluster <- selected_clusters$cluster[[i]]
    for (seed in seeds) {
      for (sample_size in refinement_sample_sizes) {
        index <- index + 1L
        sample_row <- sample_manifest[
          sample_manifest$solution_id == solution_id &
            sample_manifest$cluster == cluster &
            sample_manifest$seed == seed,
          ,
          drop = FALSE
        ]
        if (nrow(sample_row) != 1L) {
          stop(
            "Expected one staged sample for ",
            solution_id,
            " cluster ",
            cluster,
            " seed ",
            seed,
            "; found ",
            nrow(sample_row),
            ".",
            call. = FALSE
          )
        }
        task_rows[[index]] <- data.frame(
          solution_id = solution_id,
          cluster = cluster,
          seed = seed,
          requested_sample_size = sample_size,
          sample_path = sample_row$sample_path[[1L]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(
    classifications = classifications,
    selected_pairs = selected,
    tasks = if (length(task_rows) == 0L) {
      empty_tasks
    } else {
      do.call(rbind, task_rows)
    }
  )
}

run_classify <- function() {
  if (!file.exists(paths$stage)) {
    stop("Run --mode=stage first.", call. = FALSE)
  }
  sample_manifest <- readRDS(paths$sample_manifest)
  base <- read_complete_checkpoints(sample_manifest, "base")
  expected_rows <- sum(
    readRDS(paths$stage)$solution_metadata$k
  ) * length(seeds) * length(readRDS(paths$stage)$markers)
  if (nrow(base$rows) != expected_rows) {
    stop(
      "Base evidence has ",
      nrow(base$rows),
      " rows; expected ",
      expected_rows,
      ".",
      call. = FALSE
    )
  }
  key <- interaction(
    base$rows$solution_id,
    base$rows$cluster,
    base$rows$marker,
    drop = TRUE,
    lex.order = TRUE
  )
  classifications <- do.call(
    rbind,
    lapply(split(base$rows, key), classify_pair)
  )
  rownames(classifications) <- NULL
  refinement <- build_refinement_manifest(classifications, sample_manifest)
  atomic_save_rds(
    list(
      base_rows = base$rows,
      classifications = refinement$classifications,
      selected_pairs = refinement$selected_pairs,
      base_wall_seconds = base$wall_seconds,
      base_peak_rss_kb = base$peak_rss_kb
    ),
    paths$classifications
  )
  atomic_save_rds(refinement$tasks, paths$refinement_manifest)
  atomic_write_csv(
    base$rows,
    file.path(evidence_dir, "base-tests.csv")
  )
  atomic_write_csv(
    refinement$classifications,
    file.path(evidence_dir, "classifications.csv")
  )
  atomic_write_csv(
    refinement$selected_pairs,
    file.path(evidence_dir, "refinement-selected-pairs.csv")
  )
  message(
    "Classification complete; refinement tasks: ",
    nrow(refinement$tasks)
  )
}

run_refine <- function() {
  if (!file.exists(paths$refinement_manifest)) {
    stop("Run --mode=classify first.", call. = FALSE)
  }
  manifest <- task_rows_for_shard(readRDS(paths$refinement_manifest))
  message(
    "Refinement shard ",
    shard_index,
    "/",
    shard_count,
    ": ",
    nrow(manifest),
    " cluster-seed-size tasks"
  )
  for (i in seq_len(nrow(manifest))) {
    task <- manifest[i, , drop = FALSE]
    output <- refine_path(
      task$solution_id,
      task$cluster,
      task$seed,
      task$requested_sample_size
    )
    if (file.exists(output)) {
      next
    }
    sample_object <- readRDS(task$sample_path)
    timing <- system.time({
      results <- run_marker_tests(
        sample_object,
        task$requested_sample_size
      )
    })
    atomic_save_rds(
      list(
        task = task,
        results = results,
        wall_seconds = unname(timing[["elapsed"]]),
        peak_rss_kb = read_status_kb(),
        acr_B = acr_B
      ),
      output
    )
  }
}

pair_key <- function(data) {
  paste(
    data$solution_id,
    data$cluster,
    data$marker,
    sep = "\u001f"
  )
}

empty_refinement_statuses <- function() {
  data.frame(
    solution_id = character(),
    cluster = integer(),
    cluster_events = integer(),
    marker = character(),
    requested_sample_size = integer(),
    status = character(),
    reason = character(),
    valid_seeds = integer(),
    both_reject_seeds = integer(),
    dip_reject_seeds = integer(),
    acr_reject_seeds = integer(),
    min_dip_q = numeric(),
    min_acr_q = numeric(),
    stringsAsFactors = FALSE
  )
}

integrate_refinement_status <- function(base_classifications,
                                        selected_pairs,
                                        refinement_rows) {
  final <- base_classifications
  final$base_status <- final$status
  final$base_reason <- final$reason
  final$refinement_status_n1000 <- NA_character_
  final$refinement_status_n5000 <- NA_character_
  final$refinement_stable <- NA
  final$final_status_reason <- "base_2000_event_evidence"

  if (nrow(selected_pairs) == 0L) {
    sensitivity <- final[FALSE, c(
      "solution_id",
      "cluster",
      "marker",
      "base_status",
      "refinement_status_n1000",
      "refinement_status_n5000",
      "status",
      "refinement_stable",
      "final_status_reason"
    )]
    names(sensitivity)[names(sensitivity) == "status"] <- "final_status"
    return(list(
      classifications = final,
      refinement_statuses = empty_refinement_statuses(),
      sensitivity = sensitivity
    ))
  }

  selected_keys <- pair_key(selected_pairs)
  target_rows <- refinement_rows[pair_key(refinement_rows) %in% selected_keys, ]
  expected_target_rows <- nrow(selected_pairs) *
    length(seeds) *
    length(refinement_sample_sizes)
  if (nrow(target_rows) != expected_target_rows) {
    stop(
      "Refinement evidence has ",
      nrow(target_rows),
      " targeted rows; expected ",
      expected_target_rows,
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(target_rows[c(
    "solution_id",
    "cluster",
    "marker",
    "seed",
    "requested_sample_size"
  )])) {
    stop("Refinement evidence contains duplicate targeted combinations.", call. = FALSE)
  }

  group_key <- interaction(
    target_rows$solution_id,
    target_rows$cluster,
    target_rows$marker,
    target_rows$requested_sample_size,
    drop = TRUE,
    lex.order = TRUE
  )
  statuses <- do.call(
    rbind,
    lapply(split(target_rows, group_key), function(rows) {
      classified <- classify_pair(rows)
      classified$requested_sample_size <- rows$requested_sample_size[[1L]]
      classified
    })
  )
  rownames(statuses) <- NULL
  statuses <- statuses[c(
    "solution_id",
    "cluster",
    "cluster_events",
    "marker",
    "requested_sample_size",
    "status",
    "reason",
    "valid_seeds",
    "both_reject_seeds",
    "dip_reject_seeds",
    "acr_reject_seeds",
    "min_dip_q",
    "min_acr_q"
  )]

  for (i in seq_len(nrow(selected_pairs))) {
    selected <- selected_pairs[i, , drop = FALSE]
    key <- pair_key(selected)
    final_index <- which(pair_key(final) == key)
    if (length(final_index) != 1L) {
      stop("Could not map a selected pair back to one base classification.", call. = FALSE)
    }
    pair_statuses <- statuses[pair_key(statuses) == key, ]
    size_index <- match(
      refinement_sample_sizes,
      pair_statuses$requested_sample_size
    )
    if (anyNA(size_index)) {
      stop("A selected pair is missing a required refinement sample size.", call. = FALSE)
    }
    pair_statuses <- pair_statuses[
      size_index,
      ,
      drop = FALSE
    ]
    if (nrow(pair_statuses) != length(refinement_sample_sizes) ||
        any(pair_statuses$requested_sample_size != refinement_sample_sizes)) {
      stop("A selected pair is missing a required refinement sample size.", call. = FALSE)
    }

    status_1000 <- pair_statuses$status[
      pair_statuses$requested_sample_size == 1000L
    ][[1L]]
    status_5000 <- pair_statuses$status[
      pair_statuses$requested_sample_size == 5000L
    ][[1L]]
    base_status <- final$base_status[[final_index]]
    all_statuses <- c(base_status, status_1000, status_5000)
    final$refinement_status_n1000[[final_index]] <- status_1000
    final$refinement_status_n5000[[final_index]] <- status_5000
    final$refinement_stable[[final_index]] <- length(unique(all_statuses)) == 1L

    if (any(all_statuses == "unresolved_discrete")) {
      final$status[[final_index]] <- "unresolved_discrete"
      final$reason[[final_index]] <- paste0(
        "required refinement unresolved; statuses n1000/n2000/n5000=",
        status_1000,
        "/",
        base_status,
        "/",
        status_5000
      )
      final$final_status_reason[[final_index]] <- "required_refinement_unresolved"
    } else if (base_status == "ambiguous") {
      final$status[[final_index]] <- "ambiguous"
      final$reason[[final_index]] <- paste0(
        final$base_reason[[final_index]],
        "; refinement statuses n1000/n5000=",
        status_1000,
        "/",
        status_5000
      )
      final$final_status_reason[[final_index]] <- "base_evidence_ambiguous"
    } else if (length(unique(all_statuses)) == 1L) {
      final$status[[final_index]] <- base_status
      final$final_status_reason[[final_index]] <-
        "stable_across_1000_2000_5000_events"
    } else {
      final$status[[final_index]] <- "ambiguous"
      final$reason[[final_index]] <- paste0(
        "sample-size sensitivity; statuses n1000/n2000/n5000=",
        status_1000,
        "/",
        base_status,
        "/",
        status_5000
      )
      final$final_status_reason[[final_index]] <- "sample_size_instability"
    }
  }

  selected_final <- final[pair_key(final) %in% selected_keys, ]
  sensitivity <- selected_final[c(
    "solution_id",
    "cluster",
    "marker",
    "base_status",
    "refinement_status_n1000",
    "refinement_status_n5000",
    "status",
    "refinement_stable",
    "final_status_reason"
  )]
  names(sensitivity)[names(sensitivity) == "status"] <- "final_status"
  list(
    classifications = final,
    refinement_statuses = statuses,
    sensitivity = sensitivity
  )
}

pairwise_solution_ari <- function(partitions) {
  ids <- names(partitions)
  rows <- list()
  index <- 0L
  for (i in seq_len(length(ids) - 1L)) {
    for (j in seq.int(i + 1L, length(ids))) {
      index <- index + 1L
      rows[[index]] <- data.frame(
        solution_a = ids[[i]],
        solution_b = ids[[j]],
        ari = adjusted_rand_index(partitions[[i]], partitions[[j]]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

status_summaries <- function(classifications) {
  status_levels <- c(
    "detected_multimodality",
    "ambiguous",
    "no_detected_multimodality",
    "unresolved_discrete"
  )
  rows <- list()
  index <- 0L
  for (solution_id in unique(classifications$solution_id)) {
    solution <- classifications[
      classifications$solution_id == solution_id,
      ,
      drop = FALSE
    ]
    total_weight <- sum(solution$cluster_events)
    for (status in status_levels) {
      index <- index + 1L
      selected <- solution$status == status
      rows[[index]] <- data.frame(
        solution_id = solution_id,
        status = status,
        unweighted_fraction = sum(selected) / nrow(solution),
        event_weighted_fraction = sum(
          solution$cluster_events[selected]
        ) / total_weight,
        count = sum(selected),
        stringsAsFactors = FALSE
      )
    }
  }
  status_fraction <- do.call(rbind, rows)

  cluster_rows <- list()
  index <- 0L
  cluster_key <- interaction(
    classifications$solution_id,
    classifications$cluster,
    drop = TRUE
  )
  for (cluster in split(classifications, cluster_key)) {
    index <- index + 1L
    cluster_rows[[index]] <- data.frame(
      solution_id = cluster$solution_id[[1L]],
      cluster = cluster$cluster[[1L]],
      cluster_events = cluster$cluster_events[[1L]],
      no_detected_or_unresolved = !any(cluster$status %in% c(
        "detected_multimodality",
        "unresolved_discrete"
      )),
      strict_no_flagged_marker = all(
        cluster$status == "no_detected_multimodality"
      ),
      stringsAsFactors = FALSE
    )
  }
  cluster_summary <- do.call(rbind, cluster_rows)
  solution_cluster <- do.call(
    rbind,
    lapply(
      split(cluster_summary, cluster_summary$solution_id),
      function(solution) {
        data.frame(
          solution_id = solution$solution_id[[1L]],
          fraction_clusters_no_detected_or_unresolved = mean(
            solution$no_detected_or_unresolved
          ),
          strict_fraction_clusters_no_flagged_marker = mean(
            solution$strict_no_flagged_marker
          ),
          strict_solution_no_detected_multimodality = all(
            solution$strict_no_flagged_marker
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  list(
    status_fraction = status_fraction,
    cluster_summary = cluster_summary,
    solution_cluster = solution_cluster
  )
}

status_palette <- c(
  detected_multimodality = "#D55E00",
  ambiguous = "#E69F00",
  no_detected_multimodality = "#D9D9D9",
  unresolved_discrete = "#0072B2"
)

status_symbol <- c(
  detected_multimodality = "D",
  ambiguous = "A",
  no_detected_multimodality = "",
  unresolved_discrete = "U"
)

plot_status_heatmaps <- function(classifications, path) {
  grDevices::pdf(path, width = 12, height = 10, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (solution_id in unique(classifications$solution_id)) {
    data <- classifications[
      classifications$solution_id == solution_id,
      ,
      drop = FALSE
    ]
    data$cluster_label <- factor(
      data$cluster,
      levels = rev(sort(unique(data$cluster)))
    )
    data$marker <- factor(data$marker, levels = unique(data$marker))
    data$symbol <- unname(status_symbol[data$status])
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = marker, y = cluster_label, fill = status)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.15) +
      ggplot2::geom_text(
        ggplot2::aes(label = symbol),
        size = 1.7,
        color = "black"
      ) +
      ggplot2::scale_fill_manual(
        values = status_palette,
        breaks = names(status_palette),
        drop = FALSE
      ) +
      ggplot2::labs(
        title = solution_id,
        subtitle = "D = detected, A = ambiguous, U = unresolved; blank = no detected multimodality",
        x = "Marker",
        y = "Cluster",
        fill = "Status"
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom"
      )
    print(plot)
  }
}

distribution_plot <- function(pair, sample_manifest, base_rows, title_prefix) {
  sample_rows <- sample_manifest[
    sample_manifest$solution_id == pair$solution_id &
      sample_manifest$cluster == pair$cluster,
    ,
    drop = FALSE
  ]
  data_rows <- lapply(seq_len(nrow(sample_rows)), function(i) {
    sample_object <- readRDS(sample_rows$sample_path[[i]])
    values <- sample_object$values[
      seq_len(min(base_sample_size, nrow(sample_object$values))),
      pair$marker
    ]
    data.frame(
      value = values,
      seed = factor(sample_object$seed, levels = seeds)
    )
  })
  plot_data <- do.call(rbind, data_rows)
  evidence <- base_rows[
    base_rows$solution_id == pair$solution_id &
      base_rows$cluster == pair$cluster &
      base_rows$marker == pair$marker,
    ,
    drop = FALSE
  ]
  subtitle <- paste0(
    "dip q: ",
    paste(format(signif(evidence$dip_q, 3)), collapse = ", "),
    "; ACR q: ",
    paste(format(signif(evidence$acr_q, 3)), collapse = ", ")
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(
      bins = 50,
      fill = "#56B4E9",
      color = "white",
      linewidth = 0.15
    ) +
    ggplot2::facet_wrap(ggplot2::vars(seed), ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = paste0(
        title_prefix,
        ": ",
        pair$solution_id,
        " / cluster ",
        pair$cluster,
        " / ",
        pair$marker
      ),
      subtitle = subtitle,
      x = "Complete transformed marker value",
      y = "Sampled events"
    ) +
    ggplot2::theme_minimal(base_size = 9)
}

plot_pair_pages <- function(pairs,
                            sample_manifest,
                            base_rows,
                            path,
                            title_prefix) {
  grDevices::pdf(path, width = 8, height = 9, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (nrow(pairs) == 0L) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste("No", title_prefix, "pairs"))
    return(invisible(path))
  }
  for (i in seq_len(nrow(pairs))) {
    print(distribution_plot(
      pairs[i, , drop = FALSE],
      sample_manifest,
      base_rows,
      title_prefix
    ))
  }
  invisible(path)
}

run_summarize <- function() {
  stage <- readRDS(paths$stage)
  classification_data <- readRDS(paths$classifications)
  refinement_manifest <- readRDS(paths$refinement_manifest)
  refinement <- read_complete_checkpoints(refinement_manifest, "refine")
  refinement$rows$in_refinement_target <- rep(FALSE, nrow(refinement$rows))
  selected_pairs <- classification_data$selected_pairs
  target_key <- paste(
    selected_pairs$solution_id,
    selected_pairs$cluster,
    selected_pairs$marker,
    sep = "\u001f"
  )
  refinement_key <- paste(
    refinement$rows$solution_id,
    refinement$rows$cluster,
    refinement$rows$marker,
    sep = "\u001f"
  )
  refinement$rows$in_refinement_target <- refinement_key %in% target_key

  expected_refinement_combinations <- nrow(refinement_manifest) *
    length(stage$markers)
  observed_refinement_combinations <- nrow(refinement$rows)
  if (observed_refinement_combinations != expected_refinement_combinations) {
    stop("Final refinement evidence is incomplete.", call. = FALSE)
  }
  if (nrow(refinement$rows) > 0L && anyDuplicated(
    refinement$rows[c(
      "solution_id",
      "cluster",
      "seed",
      "marker",
      "requested_sample_size"
    )]
  )) {
    stop("Final refinement evidence contains duplicate combinations.", call. = FALSE)
  }
  integrated <- integrate_refinement_status(
    classification_data$classifications,
    selected_pairs,
    refinement$rows
  )
  final_classifications <- integrated$classifications
  summaries <- status_summaries(final_classifications)
  node_ari <- pairwise_solution_ari(stage$partitions)
  partition_rows <- do.call(
    rbind,
    lapply(names(stage$partitions), function(solution_id) {
      data.frame(
        solution_id = solution_id,
        node = seq_along(stage$partitions[[solution_id]]),
        cluster = stage$partitions[[solution_id]],
        stringsAsFactors = FALSE
      )
    })
  )

  expected_base_combinations <- sum(stage$solution_metadata$k) *
    length(stage$markers) *
    length(seeds)
  observed_base_combinations <- nrow(classification_data$base_rows)
  if (observed_base_combinations != expected_base_combinations) {
    stop("Final base evidence is incomplete.", call. = FALSE)
  }
  if (anyDuplicated(
    classification_data$base_rows[
      c("solution_id", "cluster", "seed", "marker")
    ]
  )) {
    stop("Final base evidence contains duplicate combinations.", call. = FALSE)
  }

  plot_status_heatmaps(
    final_classifications,
    file.path(evidence_dir, "status-heatmaps.pdf")
  )
  detected <- final_classifications[
    final_classifications$status == "detected_multimodality",
    ,
    drop = FALSE
  ]
  if (nrow(detected) > 0L) {
    detected$strength <- -log10(
      pmax(detected$min_dip_q, .Machine$double.xmin)
    ) + -log10(pmax(detected$min_acr_q, .Machine$double.xmin))
    strongest <- head(detected[order(-detected$strength), ], 20L)
  } else {
    strongest <- detected
  }
  ambiguous <- final_classifications[
    final_classifications$status == "ambiguous",
    ,
    drop = FALSE
  ]
  audited_passes <- selected_pairs[
    selected_pairs$refinement_reason == "stratified_5_percent_pass_audit",
    ,
    drop = FALSE
  ]
  plot_pair_pages(
    strongest,
    stage$sample_manifest,
    classification_data$base_rows,
    file.path(evidence_dir, "strongest-detected-distributions.pdf"),
    "Detected multimodality"
  )
  plot_pair_pages(
    ambiguous,
    stage$sample_manifest,
    classification_data$base_rows,
    file.path(evidence_dir, "all-ambiguous-distributions.pdf"),
    "Ambiguous"
  )
  plot_pair_pages(
    audited_passes,
    stage$sample_manifest,
    classification_data$base_rows,
    file.path(evidence_dir, "audited-pass-distributions.pdf"),
    "Stratified pass audit"
  )

  solution_claims <- merge(
    stage$solution_metadata,
    summaries$solution_cluster,
    by = "solution_id",
    all.x = TRUE,
    sort = FALSE
  )
  status_counts <- reshape(
    summaries$status_fraction[c("solution_id", "status", "count")],
    idvar = "solution_id",
    timevar = "status",
    direction = "wide"
  )
  solution_claims <- merge(
    solution_claims,
    status_counts,
    by = "solution_id",
    all.x = TRUE,
    sort = FALSE
  )
  solution_claims$statement <- ifelse(
    solution_claims$strict_solution_no_detected_multimodality,
    paste0(
      "No detected multimodality under both tests across all three base seeds ",
      "and required sample-size refinements; no ambiguous or unresolved ",
      "cluster-marker entries."
    ),
    paste0(
      "Residual modality evidence remains; see exact detected, ambiguous, ",
      "and unresolved counts. No global unimodality claim is supported."
    )
  )

  peak_values <- c(
    stage$stage_peak_rss_kb,
    classification_data$base_peak_rss_kb,
    refinement$peak_rss_kb,
    if (file.exists(paths$pilot)) readRDS(paths$pilot)$peak_rss_kb else NA_real_
  )
  wall_values <- c(
    stage$stage_wall_seconds,
    classification_data$base_wall_seconds,
    refinement$wall_seconds,
    if (file.exists(paths$pilot)) readRDS(paths$pilot)$wall_seconds else 0
  )
  warning_values <- c(
    classification_data$base_rows$dip_warnings,
    classification_data$base_rows$acr_warnings,
    refinement$rows$dip_warnings,
    refinement$rows$acr_warnings
  )
  warning_values <- warning_values[
    !is.na(warning_values) & nzchar(warning_values)
  ]
  interpolation_warnings <- warning_values[grepl(
    "interpol|regulariz|collaps|duplicate.*grid|ties[[:space:]]*=",
    warning_values,
    ignore.case = TRUE
  )]
  bundle <- list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    scope = paste0(
      "Univariate marker marginals only; BH controls the 27-marker family ",
      "within each cluster, method, seed, and sample size."
    ),
    test_contract = list(
      dip = "diptest::dip.test",
      acr = "multimode::modetest(mod0=1, method='ACR', B=1999)",
      q_threshold = q_threshold,
      borderline_q = borderline_q,
      no_jitter = TRUE,
      discrete_rules = list(
        minimum_unique_values = 8L,
        minimum_unique_fraction = 0.01,
        maximum_tie_fraction = 0.50,
        maximum_boundary_fraction = 0.20
      )
    ),
    stage = stage,
    pilot = if (file.exists(paths$pilot)) readRDS(paths$pilot) else NULL,
    base_tests = classification_data$base_rows,
    base_classifications = classification_data$classifications,
    classifications = final_classifications,
    refinement_selected_pairs = selected_pairs,
    refinement_tests = refinement$rows,
    refinement_statuses = integrated$refinement_statuses,
    refinement_sensitivity = integrated$sensitivity,
    summaries = summaries,
    solution_claims = solution_claims,
    node_ari = node_ari,
    peak_rss_kb = safe_max(peak_values),
    component_wall_seconds = wall_values,
    component_wall_seconds_sum = sum(wall_values, na.rm = TRUE),
    completeness = list(
      expected_base_combinations = expected_base_combinations,
      observed_base_combinations = observed_base_combinations,
      expected_refinement_combinations = expected_refinement_combinations,
      observed_refinement_combinations = observed_refinement_combinations,
      max_materialized_k = max(stage$solution_metadata$k),
      no_scan_above_100 = max(stage$solution_metadata$k) <= 100L,
      interpolation_warning_count = length(interpolation_warnings),
      all_base_rows_accounted_for = all(
        final_classifications$status %in% c(
          "detected_multimodality",
          "ambiguous",
          "no_detected_multimodality",
          "unresolved_discrete"
        )
      ),
      every_partition_materialized = all(
        vapply(
          names(stage$partitions),
          function(id) {
            length(unique(stage$partitions[[id]])) ==
              stage$solution_metadata$k[
                stage$solution_metadata$solution_id == id
              ]
          },
          logical(1)
        )
      )
    ),
    source_provenance = source_provenance(),
    dependency_versions = dependency_versions(),
    session_info = utils::sessionInfo()
  )
  atomic_save_rds(bundle, paths$final_bundle)
  atomic_write_csv(partition_rows, file.path(evidence_dir, "partitions.csv"))
  atomic_write_csv(node_ari, file.path(evidence_dir, "node-level-ari.csv"))
  atomic_write_csv(
    refinement$rows,
    file.path(evidence_dir, "refinement-tests.csv")
  )
  atomic_write_csv(
    classification_data$classifications,
    file.path(evidence_dir, "base-classifications.csv")
  )
  atomic_write_csv(
    final_classifications,
    file.path(evidence_dir, "classifications.csv")
  )
  atomic_write_csv(
    integrated$refinement_statuses,
    file.path(evidence_dir, "refinement-status-by-size.csv")
  )
  atomic_write_csv(
    integrated$sensitivity,
    file.path(evidence_dir, "refinement-sensitivity.csv")
  )
  atomic_write_csv(
    summaries$status_fraction,
    file.path(evidence_dir, "status-fractions.csv")
  )
  atomic_write_csv(
    summaries$cluster_summary,
    file.path(evidence_dir, "cluster-status-summary.csv")
  )
  atomic_write_csv(
    solution_claims,
    file.path(evidence_dir, "solution-claims.csv")
  )

  report_lines <- c(
    "# Real-model modality audit",
    "",
    paste0("- Generated: ", bundle$generated_at),
    paste0("- Model events: ", format(stage$n_events, big.mark = ",")),
    paste0("- Peak RSS across stages: ", round(bundle$peak_rss_kb / 1024^2, 2), " GiB"),
    paste0(
      "- Summed component wall time: ",
      round(bundle$component_wall_seconds_sum / 3600, 2),
      " hours"
    ),
    paste0(
      "- BH family: 27 markers within each cluster, method, seed, and sample size; ",
      "this is not global control across clusters."
    ),
    paste0(
      "- Required refinement combinations: ",
      observed_refinement_combinations,
      "/",
      expected_refinement_combinations
    ),
    paste0(
      "- Interpolation warning count: ",
      length(interpolation_warnings)
    ),
    "",
    "## Solution statements",
    "",
    unlist(lapply(seq_len(nrow(solution_claims)), function(i) {
      paste0(
        "- ",
        solution_claims$solution_id[[i]],
        ": ",
        solution_claims$statement[[i]]
      )
    })),
    "",
    "The audit addresses univariate marker marginals only. It does not establish",
    "multivariate cluster homogeneity or biological validity."
  )
  writeLines(
    report_lines,
    file.path(evidence_dir, "README.md"),
    useBytes = TRUE
  )
  message("Evidence bundle complete: ", paths$final_bundle)
}

status_report <- function() {
  files <- c(
    stage = paths$stage,
    sample_manifest = paths$sample_manifest,
    pilot = paths$pilot,
    classifications = paths$classifications,
    refinement_manifest = paths$refinement_manifest,
    final_bundle = paths$final_bundle
  )
  print(data.frame(
    artifact = names(files),
    exists = file.exists(files),
    path = unname(files),
    row.names = NULL
  ))
  if (file.exists(paths$sample_manifest)) {
    sample_manifest <- readRDS(paths$sample_manifest)
    cat(
      "Sample checkpoints: ",
      sum(file.exists(sample_manifest$sample_path)),
      "/",
      nrow(sample_manifest),
      "\n",
      sep = ""
    )
    base_paths <- expected_checkpoint_paths(sample_manifest, "base")
    cat(
      "Base checkpoints: ",
      sum(file.exists(base_paths)),
      "/",
      length(base_paths),
      "\n",
      sep = ""
    )
  }
  if (file.exists(paths$refinement_manifest)) {
    refinement_manifest <- readRDS(paths$refinement_manifest)
    refinement_paths <- expected_checkpoint_paths(refinement_manifest, "refine")
    cat(
      "Refinement checkpoints: ",
      sum(file.exists(refinement_paths)),
      "/",
      length(refinement_paths),
      "\n",
      sep = ""
    )
  }
}

run_self_test <- function() {
  make_test_rows <- function(sample_size, dip_q, acr_q) {
    data.frame(
      solution_id = "test_solution",
      cluster = 1L,
      cluster_events = 10000L,
      seed = seeds,
      marker = "marker_1",
      requested_sample_size = sample_size,
      dip_q = rep(dip_q, length(seeds)),
      acr_q = rep(acr_q, length(seeds)),
      failure_reason = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  set.seed(91)
  expected_next <- stats::runif(1L)
  set.seed(91)
  invisible(with_preserved_seed(
    safe_seed(.Machine$integer.max, "overflow", "test"),
    stats::runif(3L)
  ))
  actual_next <- stats::runif(1L)
  stopifnot(identical(actual_next, expected_next))

  base_rows <- make_test_rows(2000L, 0.01, 0.01)
  base <- classify_pair(base_rows)
  base$size_stratum <- 1L
  base$refinement_reason <- "detected_multimodality"
  refinement <- rbind(
    make_test_rows(1000L, 0.01, 0.01),
    make_test_rows(5000L, 0.01, 0.01)
  )
  stable <- integrate_refinement_status(base, base, refinement)
  stopifnot(
    identical(stable$classifications$status, "detected_multimodality"),
    isTRUE(stable$classifications$refinement_stable)
  )

  sensitive <- refinement
  sensitive$dip_q[sensitive$requested_sample_size == 5000L] <- 0.5
  sensitive$acr_q[sensitive$requested_sample_size == 5000L] <- 0.5
  unstable <- integrate_refinement_status(base, base, sensitive)
  stopifnot(
    identical(unstable$classifications$status, "ambiguous"),
    identical(
      unstable$classifications$final_status_reason,
      "sample_size_instability"
    )
  )

  unresolved <- base
  unresolved$status <- "unresolved_discrete"
  unresolved$reason <- "excessive_ties"
  unresolved$refinement_reason <- NA_character_
  manifest <- build_refinement_manifest(
    unresolved,
    data.frame(
      solution_id = character(),
      cluster = integer(),
      seed = integer(),
      sample_path = character(),
      stringsAsFactors = FALSE
    )
  )
  stopifnot(
    nrow(manifest$selected_pairs) == 0L,
    nrow(manifest$tasks) == 0L,
    length(expected_checkpoint_paths(manifest$tasks, "refine")) == 0L
  )
  message("Workflow self-test passed.")
}

switch(
  mode,
  stage = run_stage(),
  pilot = run_pilot(),
  base = run_base(),
  classify = run_classify(),
  refine = run_refine(),
  summarize = run_summarize(),
  status = status_report(),
  `self-test` = run_self_test(),
  stop(
    "Unknown mode `",
    mode,
    "`. Use stage, pilot, base, classify, refine, summarize, status, or self-test.",
    call. = FALSE
  )
)
