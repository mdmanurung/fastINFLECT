## Stability benchmark for deterministic per-node event caps on a large SOM.
##
## HISTORICAL ENGINEERING REPRODUCTION ONLY. This script deliberately preserves
## the 2026-07-24 IQR-only/non-positive-excluding contract so its archived timing
## and cap-sensitivity result can be reproduced. It is invalid as modality
## evidence. Use validate-real-model-modality.R for the 2.0 scientific audit.
##
## Usage:
##   Rscript data-raw/benchmark-large-model-cap.R /tmp/inflect-cap-benchmark.rds
##
## Override the model path with INFLECT_LARGE_MODEL. The script loads the model
## once, reuses one hierarchy, and compares uncapped scoring with caps of 1000,
## 2000, and 5000 across seeds 1, 42, and 2026.

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args)) args[[1L]] else {
  file.path(tempdir(), "inflect-cap-benchmark.rds")
}
model_path <- Sys.getenv(
  "INFLECT_LARGE_MODEL",
  unset = paste0(
    "/exports/para-lipg-hpc/Xuran/results/bmv/aurora_exvivo/cluster/",
    "model/xyf_som_model_with_LD_filter.rds"
  )
)
if (!file.exists(model_path)) {
  stop("Large-model artifact not found: ", model_path, call. = FALSE)
}

pkgload::load_all(".", quiet = TRUE)
schedule <- seq.int(25L, 100L, by = 5L)
caps <- c(NA_integer_, 1000L, 2000L, 5000L)
seeds <- c(1L, 42L, 2026L)

read_status_kb <- function(field) {
  status <- readLines("/proc/self/status", warn = FALSE)
  hit <- grep(paste0("^", field, ":"), status, value = TRUE)
  if (!length(hit)) {
    return(NA_real_)
  }
  as.numeric(sub(".*?([0-9]+) kB.*", "\\1", hit[[1L]]))
}

message("Loading model: ", model_path)
load_time <- system.time(model <- readRDS(model_path))[["elapsed"]]
view <- as_inflect_som(model)
rm(model)
invisible(gc())
message(
  "Loaded ",
  nrow(view$data),
  " events, ",
  view$map$nNodes,
  " nodes, ",
  ncol(view$data),
  " data columns in ",
  round(load_time, 1),
  " s"
)

metaclustering <- iteration.metacluster(
  view,
  set.i = schedule,
  multicore = FALSE
)

runs <- list()
run_index <- 0L
for (cap in caps) {
  run_seeds <- if (is.na(cap)) 1L else seeds
  for (seed in run_seeds) {
    run_index <- run_index + 1L
    label <- if (is.na(cap)) "uncapped" else paste0("cap_", cap)
    message("Scoring ", label, " with seed ", seed)
    cap_value <- if (is.na(cap)) NULL else cap
    timing <- system.time({
      qc <- iteration.QC(
        FlowSOM.results = view,
        metaclustering.list = metaclustering,
        set.i = schedule,
        multicore = TRUE,
        cores = 2L,
        zeroes.in = FALSE,
        uniform.test = "spread",
        max.events.per.node = cap_value,
        seed = seed,
        verbose = FALSE
      )
      curve <- QC.to.curve(qc, basedata = "Curve")
      selection <- .inflect_selection(
        collection.U = curve$collection.U,
        lfunction = curve$lfunction,
        target = 0.95,
        fittedcurve = curve$fittedcurve
      )
    })
    runs[[run_index]] <- list(
      cap = cap,
      label = label,
      seed = seed,
      elapsed_seconds = unname(timing[["elapsed"]]),
      scores = qc$U.set,
      selection = selection,
      qc_provenance = qc$provenance,
      process_high_water_kb = read_status_kb("VmHWM")
    )
    rm(qc, curve, selection)
    invisible(gc())
  }
}

baseline <- runs[[1L]]
recommendation_methods <- c("inflection", "kneedle", "threshold")
summary_rows <- lapply(runs, function(run) {
  score_delta <- run$scores$Unimodality - baseline$scores$Unimodality
  run_k <- stats::setNames(run$selection$k, run$selection$method)
  baseline_k <- stats::setNames(
    baseline$selection$k,
    baseline$selection$method
  )
  k_delta <- abs(run_k[recommendation_methods] - baseline_k[recommendation_methods])
  data.frame(
    label = run$label,
    cap = run$cap,
    seed = run$seed,
    elapsed_seconds = run$elapsed_seconds,
    retained_events = run$qc_provenance$event_sampling$retained_events,
    max_abs_score_deviation_pp = max(abs(score_delta)),
    inflection_k = unname(run_k[["inflection"]]),
    kneedle_k = unname(run_k[["kneedle"]]),
    threshold_k = unname(run_k[["threshold"]]),
    max_abs_recommendation_delta_k = if (all(is.na(k_delta))) {
      NA_real_
    } else {
      max(k_delta, na.rm = TRUE)
    },
    run_gate = max(abs(score_delta)) <= 1 &&
      all(is.na(k_delta) | k_delta <= 5),
    process_high_water_kb = run$process_high_water_kb
  )
})
run_summary <- do.call(rbind, summary_rows)

cap_stability <- do.call(rbind, lapply(c(1000L, 2000L, 5000L), function(cap) {
  cap_runs <- runs[vapply(runs, function(x) identical(x$cap, cap), logical(1))]
  score_matrix <- do.call(
    cbind,
    lapply(cap_runs, function(x) x$scores$Unimodality)
  )
  selections <- do.call(
    cbind,
    lapply(cap_runs, function(x) {
      stats::setNames(x$selection$k, x$selection$method)[recommendation_methods]
    })
  )
  score_seed_range <- apply(score_matrix, 1L, function(x) diff(range(x)))
  recommendation_seed_range <- apply(selections, 1L, function(x) {
    if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
  })
  rows <- run_summary[run_summary$cap == cap, , drop = FALSE]
  data.frame(
    cap = cap,
    max_score_range_across_seeds_pp = max(score_seed_range),
    max_recommendation_range_across_seeds_k = max(
      recommendation_seed_range,
      na.rm = TRUE
    ),
    all_runs_pass_baseline_gate = all(rows$run_gate),
    stability_gate = max(score_seed_range) <= 1 &&
      all(is.na(recommendation_seed_range) | recommendation_seed_range <= 5)
  )
}))
cap_stability$acceptable <- with(
  cap_stability,
  all_runs_pass_baseline_gate & stability_gate
)
acceptable_caps <- cap_stability$cap[cap_stability$acceptable]
selected_cap <- if (length(acceptable_caps)) min(acceptable_caps) else NA_integer_

result <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
  model_path = model_path,
  model_size_bytes = file.info(model_path)$size,
  model_load_seconds = unname(load_time),
  n_events = nrow(view$data),
  n_nodes = view$map$nNodes,
  n_data_columns = ncol(view$data),
  schedule = schedule,
  uniform_test = "spread",
  zeroes_in = FALSE,
  cores = 2L,
  seeds = seeds,
  caps = caps,
  run_summary = run_summary,
  cap_stability = cap_stability,
  selected_cap = selected_cap,
  runs = runs,
  process_high_water_kb = read_status_kb("VmHWM"),
  session_info = utils::sessionInfo()
)
saveRDS(result, output_path)
print(run_summary)
print(cap_stability)
message("Selected cap: ", if (is.na(selected_cap)) "NULL" else selected_cap)
message("Saved benchmark: ", normalizePath(output_path, mustWork = TRUE))
