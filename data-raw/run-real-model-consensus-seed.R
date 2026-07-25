#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
seed_arg <- grep("^--seed=", args, value = TRUE)
if (length(seed_arg) != 1L) {
  stop("Supply exactly one --seed=<integer> argument.", call. = FALSE)
}
seed <- suppressWarnings(as.integer(sub("^--seed=", "", seed_arg)))
if (is.na(seed) || seed < 1L) {
  stop("`--seed` must be a positive integer.", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not resolve the helper script path.", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
work_dir <- Sys.getenv(
  "INFLECT_SELECTION_WORKDIR",
  file.path(repo_root, "data-raw", ".real-model-selection-work")
)
full_path <- file.path(work_dir, "full-inflect-selection.rds")
output_dir <- file.path(work_dir, "consensus")
output_path <- file.path(
  output_dir,
  sprintf("consensus-seed-%04d.rds", seed)
)
schedule <- 25L:100L
consensus_reps <- 100L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(output_path)) {
  message("Consensus checkpoint already exists: ", output_path)
  quit(save = "no", status = 0L)
}
if (!file.exists(full_path)) {
  stop("Full INFLECT checkpoint does not exist: ", full_path, call. = FALSE)
}

read_status_kb <- function() {
  status <- readLines("/proc/self/status", warn = FALSE)
  values <- sub(".*:[[:space:]]*([0-9]+).*", "\\1", status)
  names(values) <- sub(":.*", "", status)
  requested <- c("VmHWM", "VmRSS")
  available <- requested[requested %in% names(values)]
  if (length(available) == 0L) {
    return(NA_real_)
  }
  as.numeric(values[[available[[1L]]]])
}

atomic_save_rds <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    stop("Failed to publish consensus checkpoint: ", path, call. = FALSE)
  }
  invisible(path)
}

consensus_cdf_area <- function(consensus_matrix) {
  values <- consensus_matrix[upper.tri(consensus_matrix)]
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(NA_real_)
  }
  1 - mean(values)
}

full <- readRDS(full_path)
if (!identical(as.integer(full$schedule), schedule)) {
  stop("The full checkpoint schedule is not the literal 25:100 audit.", call. = FALSE)
}
codes <- as.matrix(full$inflect_codes)
if (nrow(codes) != 900L || any(!is.finite(codes))) {
  stop("The full checkpoint code matrix is invalid.", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
results <- ConsensusClusterPlus::ConsensusClusterPlus(
  t(codes),
  maxK = max(schedule),
  reps = consensus_reps,
  pItem = 0.9,
  pFeature = 1,
  title = file.path(work_dir, sprintf("consensus-helper-%04d", seed)),
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
if (!all(vapply(
  seq_along(schedule),
  function(i) length(unique(classes[[i]])) == schedule[[i]],
  logical(1)
))) {
  stop("Consensus helper produced a partition with the wrong k.", call. = FALSE)
}

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
  peak_rss_kb = read_status_kb(),
  producer = list(
    implementation = "parallel exact-seed helper",
    script_file = script_file,
    script_md5 = unname(tools::md5sum(script_file)),
    full_checkpoint = normalizePath(full_path, mustWork = TRUE),
    full_checkpoint_md5 = unname(tools::md5sum(full_path)),
    parameters = list(
      max_k = max(schedule),
      reps = consensus_reps,
      p_item = 0.9,
      p_feature = 1,
      cluster_alg = "hc",
      distance = "euclidean"
    ),
    dependency_versions = c(
      R = as.character(getRversion()),
      ConsensusClusterPlus = as.character(
        utils::packageVersion("ConsensusClusterPlus")
      )
    ),
    session_info = utils::sessionInfo()
  )
)
atomic_save_rds(checkpoint, output_path)
message("Consensus checkpoint complete: ", output_path)
