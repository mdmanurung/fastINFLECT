## data-raw/make-benchmark-cache.R
##
## HISTORICAL 1.0 CACHE GENERATOR. The package vignette uses only its
## engineering timing fields. Its zeroes.in=FALSE aggregate output is invalid as
## modality evidence under the 2.0 interpretation contract.
##
## Generates inst/extdata/benchmark_cache.rds for the vignette
## "Benchmarking fastINFLECT against consensus metaclustering".
##
## It benchmarks the common task "scan a range of k and score cluster quality"
## two ways on the bundled Levine32 SOM:
##   * fastINFLECT's memoised fast engine (hierarchical cut + unimodality QC), and
##   * FlowSOM consensus metaclustering (ConsensusClusterPlus) run per k.
## Because ConsensusClusterPlus is slow, this script is run ONCE offline and the
## vignette merely loads the cache. Timings are machine-specific and recorded
## here for illustration.
##
## Usage (from the package root, with the R4_51 / scale_r conda env):
##   Rscript data-raw/make-benchmark-cache.R

suppressPackageStartupMessages({
  library(FlowSOM)
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(fastINFLECT)
  }
})

set.seed(42)

existing_cache_path <- file.path("inst", "extdata", "benchmark_cache.rds")
existing_cache <- if (file.exists(existing_cache_path)) readRDS(existing_cache_path) else NULL
refresh_consensus <- identical(Sys.getenv("INFLECT_REFRESH_CONSENSUS"), "1")

# ── 1. Data ──────────────────────────────────────────────────────────────────
message("Loading Levine32sample dataset...")
data_path <- system.file("extdata", "Levine32sample.Rdata", package = "fastINFLECT")
stopifnot(nchar(data_path) > 0)
load(data_path)                     # object 'dataset' (class FlowSOM)
codes <- dataset$map$codes          # 375 x 32 SOM node prototypes

## Unimodality score, matching iteration.QC's aggregation exactly.
uni_score <- function(mc) {
  m <- FlowSOMQC(dataset, mc, zeroes.in = FALSE, verbose = FALSE)
  sum(m, na.rm = TRUE) * 100 / prod(dim(m))
}

## median wall-clock (seconds) of evaluating expr `reps` times.
time_median <- function(expr, reps = 1L) {
  expr <- substitute(expr)
  env <- parent.frame()
  stats::median(vapply(seq_len(reps), function(i) {
    as.numeric(system.time(eval(expr, env))[["elapsed"]])
  }, numeric(1)))
}

pkg_version <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(utils::packageVersion(pkg))
  } else {
    NA_character_
  }
}

object_hash <- function(object) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

file_hashes <- function(paths) {
  paths <- paths[file.exists(paths)]
  stats::setNames(unname(tools::md5sum(paths)), paths)
}

cache_provenance_matches <- function(cache, provenance) {
  is.list(cache) && identical(cache$cache_provenance, provenance)
}

make_qc_context <- function(som, only.clustering.markers = TRUE,
                            acquired_markers = NULL) {
  view <- getFromNamespace("as_inflect_som", "fastINFLECT")(som)
  prep <- getFromNamespace(".inflect_prepare_qc", "fastINFLECT")(
    view = view,
    only.clustering.markers = only.clustering.markers,
    acquired_markers = acquired_markers
  )
  list(
    view = view,
    prep = prep,
    fast_pvalue = getFromNamespace(".inflect_dip_pvalue", "fastINFLECT")
  )
}

marker_expression <- function(values, zeroes.in) {
  if (isFALSE(zeroes.in)) {
    me <- values[values > 0]
    if (length(me) < 5L) {
      me <- c(rep(0, 5L - length(me)), me)
    }
    return(me)
  }
  values
}

legacy_cell_score <- function(values,
                              zeroes.in = FALSE,
                              uniform.test = "both",
                              th.pvalue = 0.05,
                              th.IQR = 2) {
  me <- marker_expression(values, zeroes.in)
  p.value <- NA_real_
  marker.iqr <- NA_real_
  uniform <- TRUE
  if (uniform.test != "spread") {
    p.value <- diptest::dip.test(me)$p.value
    uniform <- uniform && p.value >= th.pvalue
  }
  if (uniform.test != "unimodality") {
    marker.quantile <- stats::quantile(me, names = FALSE)
    marker.iqr <- marker.quantile[[4L]] - marker.quantile[[2L]]
    uniform <- uniform && marker.iqr < th.IQR
  }
  list(
    values = me,
    qc_pass = uniform,
    dip_pvalue = p.value,
    iqr = marker.iqr,
    n_events = length(me)
  )
}

fast_cell_score <- function(values,
                            context,
                            zeroes.in = FALSE,
                            uniform.test = "both",
                            th.pvalue = 0.05,
                            th.IQR = 2) {
  me <- marker_expression(values, zeroes.in)
  p.value <- NA_real_
  marker.iqr <- NA_real_
  uniform <- TRUE
  if (uniform.test != "spread") {
    p.value <- context$fast_pvalue(diptest::dip(me), length(me))
    uniform <- uniform && p.value >= th.pvalue
  }
  if (uniform.test != "unimodality") {
    marker.quantile <- stats::quantile(me, names = FALSE)
    marker.iqr <- marker.quantile[[4L]] - marker.quantile[[2L]]
    uniform <- uniform && marker.iqr < th.IQR
  }
  list(
    values = me,
    qc_pass = uniform,
    dip_pvalue = p.value,
    iqr = marker.iqr,
    n_events = length(me)
  )
}

cluster_rows_for <- function(context, metaclustering) {
  event_cluster <- metaclustering[context$prep$mapping]
  split(
    seq_len(nrow(context$prep$data)),
    factor(event_cluster, levels = seq_len(max(metaclustering)))
  )
}

benchmark_legacy_inflect <- function(context,
                                     metaclustering.list,
                                     k_range,
                                     unimodality_by_k,
                                     sample_per_k = 8L,
                                     zeroes.in = FALSE,
                                     uniform.test = "both",
                                     th.pvalue = 0.05,
                                     th.IQR = 2) {
  set.seed(4201)
  markers <- context$prep$ordered.markers
  rows <- vector("list", length(k_range))

  for (idx in seq_along(k_range)) {
    k <- k_range[[idx]]
    mc <- metaclustering.list[[as.character(k)]]
    cluster.rows <- cluster_rows_for(context, mc)
    cells <- expand.grid(
      cluster = seq_len(max(mc)),
      marker = markers,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    sample_idx <- sort(sample(seq_len(nrow(cells)), min(sample_per_k, nrow(cells))))
    timings <- numeric(length(sample_idx))

    for (j in seq_along(sample_idx)) {
      cell <- cells[sample_idx[[j]], ]
      event_rows <- cluster.rows[[as.character(cell$cluster)]]
      values <- context$prep$data[event_rows, cell$marker]
      timings[[j]] <- as.numeric(system.time(invisible(legacy_cell_score(
        values,
        zeroes.in = zeroes.in,
        uniform.test = uniform.test,
        th.pvalue = th.pvalue,
        th.IQR = th.IQR
      )))[["elapsed"]])
    }

    rows[[idx]] <- data.frame(
      k = k,
      seconds = mean(timings) * nrow(cells),
      seconds_se = stats::sd(timings) / sqrt(length(timings)) * nrow(cells),
      sampled_seconds = sum(timings),
      sampled_tests = length(sample_idx),
      n_cluster_marker_tests = nrow(cells),
      mean_seconds_per_test = mean(timings),
      unimodality = unimodality_by_k[as.character(k)],
      method = "Original INFLECT (estimated)",
      timing_type = "estimated_from_sampled_diptest_cells",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

summarise_marker_pairs <- function(context,
                                   som,
                                   metaclustering,
                                   method,
                                   k,
                                   zeroes.in = FALSE,
                                   uniform.test = "both",
                                   th.pvalue = 0.05,
                                   th.IQR = 2) {
  accuracy <- FlowSOMQC(
    som,
    as.integer(metaclustering),
    zeroes.in = zeroes.in,
    uniform.test = uniform.test,
    th.pvalue = th.pvalue,
    th.IQR = th.IQR,
    verbose = FALSE
  )
  cluster.rows <- cluster_rows_for(context, metaclustering)
  markers <- colnames(accuracy)
  summaries <- vector("list", nrow(accuracy) * ncol(accuracy))
  expressions <- vector("list", length(summaries))
  event_ids <- vector("list", length(summaries))
  row_id <- 0L

  for (cluster in seq_len(nrow(accuracy))) {
    rows <- cluster.rows[[as.character(cluster)]]
    for (marker in markers) {
      row_id <- row_id + 1L
      if (length(rows) > 1L) {
        stats <- fast_cell_score(
          context$prep$data[rows, marker],
          context = context,
          zeroes.in = zeroes.in,
          uniform.test = uniform.test,
          th.pvalue = th.pvalue,
          th.IQR = th.IQR
        )
        qc_pass <- isTRUE(accuracy[as.character(cluster), marker])
      } else {
        stats <- list(
          values = numeric(0),
          dip_pvalue = NA_real_,
          iqr = NA_real_,
          n_events = 0L
        )
        qc_pass <- NA
      }
      expressions[[row_id]] <- stats$values
      event_ids[[row_id]] <- rows
      summaries[[row_id]] <- data.frame(
        row_id = row_id,
        method = method,
        k = k,
        cluster = cluster,
        marker = marker,
        qc_pass = qc_pass,
        dip_pvalue = stats$dip_pvalue,
        iqr = stats$iqr,
        n_events = stats$n_events,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    summary = do.call(rbind, summaries),
    expressions = expressions,
    event_ids = event_ids
  )
}

select_histogram_pairs <- function(consensus_pairs, inflect_pairs,
                                   max_inflect_panels = 2L) {
  consensus <- consensus_pairs$summary
  inflect <- inflect_pairs$summary
  failures <- consensus[
    !is.na(consensus$qc_pass) &
      !consensus$qc_pass &
      is.finite(consensus$dip_pvalue) &
      consensus$n_events >= 20L,
  ]
  failures <- failures[order(failures$dip_pvalue, -failures$n_events), ]

  for (i in seq_len(nrow(failures))) {
    standard <- failures[i, ]
    standard_events <- consensus_pairs$event_ids[[standard$row_id]]
    candidates <- inflect[
      inflect$marker == standard$marker &
        !is.na(inflect$qc_pass) &
        inflect$qc_pass,
    ]
    if (!nrow(candidates)) {
      next
    }
    overlaps <- vapply(candidates$row_id, function(id) {
      length(intersect(standard_events, inflect_pairs$event_ids[[id]]))
    }, integer(1))
    candidates$overlap_with_standard <- overlaps
    min_overlap <- max(10L, floor(0.05 * length(standard_events)))
    candidates <- candidates[candidates$overlap_with_standard >= min_overlap, ]
    if (!nrow(candidates)) {
      next
    }
    candidates <- candidates[order(-candidates$overlap_with_standard, -candidates$n_events), ]
    candidates <- utils::head(candidates, max_inflect_panels)
    standard$overlap_with_standard <- length(standard_events)
    standard$selected_role <- "standard residual bimodality"
    candidates$selected_role <- "overlapping fastINFLECT unimodal cluster"
    return(rbind(standard, candidates))
  }

  standard <- failures[1, ]
  standard$overlap_with_standard <- standard$n_events
  standard$selected_role <- "standard residual bimodality"
  candidates <- inflect[
    inflect$marker == standard$marker &
      !is.na(inflect$qc_pass) &
      inflect$qc_pass,
  ]
  candidates <- utils::head(candidates[order(-candidates$n_events), ], max_inflect_panels)
  candidates$overlap_with_standard <- NA_integer_
  candidates$selected_role <- "fastINFLECT unimodal cluster"
  rbind(standard, candidates)
}

build_marker_histogram <- function(context,
                                   som,
                                   inflect_mc,
                                   consensus_mc,
                                   k,
                                   zeroes.in = FALSE,
                                   uniform.test = "both",
                                   th.pvalue = 0.05,
                                   th.IQR = 2) {
  consensus_pairs <- summarise_marker_pairs(
    context, som, consensus_mc,
    method = "FlowSOM consensus",
    k = k,
    zeroes.in = zeroes.in,
    uniform.test = uniform.test,
    th.pvalue = th.pvalue,
    th.IQR = th.IQR
  )
  inflect_pairs <- summarise_marker_pairs(
    context, som, inflect_mc,
    method = "fastINFLECT threshold",
    k = k,
    zeroes.in = zeroes.in,
    uniform.test = uniform.test,
    th.pvalue = th.pvalue,
    th.IQR = th.IQR
  )

  selected <- select_histogram_pairs(consensus_pairs, inflect_pairs)
  selected$panel <- sprintf(
    "%s: cluster %s, %s (%s)",
    selected$method,
    selected$cluster,
    selected$marker,
    ifelse(selected$qc_pass, "passes QC", "fails QC")
  )

  histogram_rows <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, ]
    pair <- if (identical(row$method, "FlowSOM consensus")) consensus_pairs else inflect_pairs
    values <- pair$expressions[[row$row_id]]
    histogram_rows[[i]] <- data.frame(
      method = row$method,
      k = row$k,
      cluster = row$cluster,
      marker = row$marker,
      expression = values,
      qc_pass = row$qc_pass,
      dip_pvalue = row$dip_pvalue,
      iqr = row$iqr,
      n_events = row$n_events,
      overlap_with_standard = row$overlap_with_standard,
      selected_role = row$selected_role,
      panel = row$panel,
      stringsAsFactors = FALSE
    )
  }

  list(
    marker_histogram = do.call(rbind, histogram_rows),
    marker_histogram_summary = selected[, c(
      "method", "k", "cluster", "marker", "qc_pass", "dip_pvalue", "iqr",
      "n_events", "overlap_with_standard", "selected_role", "panel"
    )]
  )
}

k_range <- 5:25
context <- make_qc_context(dataset)
consensus_params <- list(
  reps = 100L,
  pItem = 0.9,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "euclidean",
  seed = 42L
)
cache_provenance <- list(
  dataset_hash = object_hash(list(
    data = dataset$data,
    mapping = dataset$map$mapping,
    colsUsed = dataset$map$colsUsed
  )),
  codes_hash = object_hash(codes),
  k_range = k_range,
  seed = 42L,
  consensus_params = consensus_params,
  package_versions = list(
    fastINFLECT = getFromNamespace("inflect_package_version", "fastINFLECT")(),
    FlowSOM = pkg_version("FlowSOM"),
    diptest = pkg_version("diptest"),
    ConsensusClusterPlus = pkg_version("ConsensusClusterPlus")
  ),
  source_hashes = file_hashes(c(
    "R/som-adapter.R",
    "R/inflect-qc-core.R",
    "R/FlowSOM-QC.R",
    "R/iteration-QC.R",
    "R/iteration-metacluster.R",
    "data-raw/make-benchmark-cache.R"
  ))
)

# ── 2. fastINFLECT: one memoised sweep over k_range ──────────────────────────
message("Timing fastINFLECT fast engine over k = 5:25 ...")
## Warm up (compile + load qDiptab) so timing reflects steady state.
invisible(iteration.QC(dataset,
                       iteration.metacluster(dataset, 5:6, multicore = FALSE),
                       5:6, multicore = FALSE, zeroes.in = FALSE, verbose = FALSE))

inflect_scan_seconds <- time_median({
  ml <- iteration.metacluster(dataset, set.i = k_range, multicore = FALSE)
  qc <- iteration.QC(dataset, ml, k_range, multicore = FALSE,
                     zeroes.in = FALSE, verbose = FALSE)
}, reps = 5L)

## Full INFLECT() result (curve + automatic k selection) for the same range.
inflect_res <- INFLECT(dataset, set.i = k_range, multicore = FALSE,
                       zeroes.in = FALSE, verbose = FALSE)

inflect_curve <- data.frame(
  k = inflect_res$collection.U$i,
  unimodality = inflect_res$collection.U$Unimodality,
  method = "fastINFLECT (hierarchical cut)",
  stringsAsFactors = FALSE
)
unimodality_by_k <- stats::setNames(inflect_curve$unimodality, inflect_curve$k)

# ── 3. Original INFLECT: sampled legacy cell timings ─────────────────────────
message("Estimating original INFLECT runtime from sampled dip.test cells ...")
legacy_inflect <- benchmark_legacy_inflect(
  context = context,
  metaclustering.list = inflect_res$metaclustering.list,
  k_range = k_range,
  unimodality_by_k = unimodality_by_k,
  sample_per_k = 8L,
  zeroes.in = FALSE
)
legacy_inflect_scan_seconds <- sum(legacy_inflect$seconds)
legacy_inflect_scan_seconds_se <- sqrt(sum(legacy_inflect$seconds_se ^ 2, na.rm = TRUE))
message(sprintf("  original INFLECT estimated scan: %.1f +/- %.1f s from %d sampled cells",
                legacy_inflect_scan_seconds, legacy_inflect_scan_seconds_se,
                sum(legacy_inflect$sampled_tests)))

# ── 4. Consensus metaclustering: one ConsensusClusterPlus per k ──────────────
if (!is.null(existing_cache) && !refresh_consensus &&
    cache_provenance_matches(existing_cache, cache_provenance) &&
    all(c("consensus_df", "amortized_df", "totals") %in% names(existing_cache))) {
  message("Reusing existing FlowSOM consensus timings with matching provenance. Set INFLECT_REFRESH_CONSENSUS=1 to retime.")
  consensus_df <- existing_cache$consensus_df
  consensus_scan_seconds <- existing_cache$totals$consensus_scan_seconds
  consensus_per_k_median <- existing_cache$totals$consensus_per_k_median
  amortized_seconds <- existing_cache$totals$amortized_seconds
  amortized_df <- existing_cache$amortized_df
} else {
  message("Timing FlowSOM consensus metaclustering plus QC scoring per k ...")
  consensus_rows <- vector("list", length(k_range))
  for (i in seq_along(k_range)) {
    k <- k_range[i]
    t_k <- as.numeric(system.time({
      mc <- suppressMessages(metaClustering_consensus(codes, k = k, seed = 42))
      unimodality <- uni_score(mc)
    })[["elapsed"]])
    consensus_rows[[i]] <- data.frame(
      k = k,
      unimodality = unimodality,
      method = "FlowSOM consensus",
      seconds = t_k,
      timing_type = "partition_plus_qc_scoring",
      stringsAsFactors = FALSE
    )
    message(sprintf("  k = %2d : %.2f s", k, t_k))
  }
  consensus_df <- do.call(rbind, consensus_rows)
  consensus_scan_seconds <- sum(consensus_df$seconds)
  consensus_per_k_median <- stats::median(consensus_df$seconds)

  ## Amortized: a single ConsensusClusterPlus run up to max(k) yields every k at
  ## once (metaClustering_consensus() just slices one such run). This is the
  ## fairest consensus baseline for scanning k; even so it does not pick k.
  message("Timing a single amortized ConsensusClusterPlus(maxK) run plus QC scoring ...")
  amortized_seconds <- as.numeric(system.time({
    ccp <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
      t(codes), maxK = max(k_range), reps = consensus_params$reps,
      pItem = consensus_params$pItem, pFeature = consensus_params$pFeature,
      title = tempfile("ccp"), plot = NULL, verbose = FALSE,
      clusterAlg = consensus_params$clusterAlg,
      distance = consensus_params$distance,
      seed = consensus_params$seed))
    amortized_df <- do.call(rbind, lapply(k_range, function(k) {
      data.frame(k = k, unimodality = uni_score(ccp[[k]]$consensusClass),
                 method = "FlowSOM consensus (single run)",
                 timing_type = "single_consensus_run_plus_qc_scoring",
                 stringsAsFactors = FALSE)
    }))
  })[["elapsed"]])
  message(sprintf("  single amortized run: %.1f s (vs %.1f s for %d separate calls)",
                  amortized_seconds, consensus_scan_seconds, length(k_range)))
}

# ── 5. Scaling with sweep density ────────────────────────────────────────────
## fastINFLECT's cost barely grows with the number of tested k (distinct SOM-node
## subtrees are bounded), whereas consensus runs one ConsensusClusterPlus per k.
message("Measuring fastINFLECT scaling across sweep sizes ...")
sweep_sets <- list(`5:25` = 5:25, `5:50` = 5:50, `5:100` = 5:100, `5:200` = 5:200)
scaling_rows <- lapply(names(sweep_sets), function(nm) {
  si <- sweep_sets[[nm]]
  s <- time_median({
    ml <- iteration.metacluster(dataset, set.i = si, multicore = FALSE)
    iteration.QC(dataset, ml, si, multicore = FALSE, zeroes.in = FALSE, verbose = FALSE)
  }, reps = 3L)
  data.frame(
    n_k = length(si),
    label = nm,
    inflect_seconds = s,
    ## consensus is linear in the number of tested k (one run each); projected
    ## from the measured median per-k cost.
    consensus_seconds_projected = length(si) * consensus_per_k_median,
    stringsAsFactors = FALSE
  )
})
scaling_df <- do.call(rbind, scaling_rows)

# ── 6. Recommended k (fastINFLECT's automatic selection) ─────────────────────
selection <- inflect_res$selection
threshold_k <- selection$k[selection$method == "threshold"]
if (!length(threshold_k) || is.na(threshold_k)) {
  threshold_k <- selection$k[selection$method == "kneedle"]
}
threshold_k <- as.integer(round(threshold_k[[1L]]))

# ── 7. Marker histograms: residual bimodality vs fastINFLECT unimodality ─────
message(sprintf("Building marker-expression histogram panels at k = %d ...", threshold_k))
consensus_threshold_mc <- suppressMessages(metaClustering_consensus(codes, k = threshold_k, seed = 42))
histogram_cache <- build_marker_histogram(
  context = context,
  som = dataset,
  inflect_mc = inflect_res$metaclustering.list[[as.character(threshold_k)]],
  consensus_mc = consensus_threshold_mc,
  k = threshold_k,
  zeroes.in = FALSE
)

# ── 8. Computational-efficiency counters ────────────────────────────────────
all_keys <- unique(unlist(lapply(inflect_res$metaclustering.list, function(mc) {
  vapply(seq_len(max(mc)), function(cl) paste0(which(mc == cl), collapse = ","), character(1))
}), use.names = FALSE))
n_markers <- length(context$prep$ordered.markers)
efficiency <- list(
  n_events = nrow(context$prep$data),
  n_nodes = as.integer(context$view$map$nNodes),
  n_markers = n_markers,
  total_cluster_marker_tests_legacy = sum(k_range) * n_markers,
  distinct_subtree_marker_tests_inflect = length(all_keys) * n_markers,
  distinct_subtrees_inflect = length(all_keys),
  memoization_work_reduction = (sum(k_range) * n_markers) / (length(all_keys) * n_markers),
  legacy_projection = list(
    sampled_tests = sum(legacy_inflect$sampled_tests),
    sample_per_k = 8L,
    seconds_se = legacy_inflect_scan_seconds_se,
    timing_type = "estimated_from_sampled_diptest_cells",
    note = paste(
      "Estimated from sampled original FlowSOMQC cluster-marker cells using",
      "diptest::dip.test(); the reported standard error reflects sampled cell",
      "timings and does not include small metaclustering or curve-fitting overhead."
    )
  )
)

# ── 9. Assemble and save ─────────────────────────────────────────────────────
cache <- list(
  cache_provenance = cache_provenance,
  k_range = k_range,
  comparison_df = rbind(
    inflect_curve[, c("k", "unimodality", "method")],
    consensus_df[, c("k", "unimodality", "method")]
  ),
  consensus_df = consensus_df,
  amortized_df = amortized_df,
  legacy_inflect = legacy_inflect,
  scaling_df = scaling_df,
  selection = selection,
  marker_histogram = histogram_cache$marker_histogram,
  marker_histogram_summary = histogram_cache$marker_histogram_summary,
  efficiency = efficiency,
  inflect_res = inflect_res,
  totals = list(
    inflect_scan_seconds = inflect_scan_seconds,
    consensus_scan_seconds = consensus_scan_seconds,
    consensus_per_k_median = consensus_per_k_median,
    amortized_seconds = amortized_seconds,
    legacy_inflect_scan_seconds = legacy_inflect_scan_seconds,
    legacy_inflect_scan_seconds_se = legacy_inflect_scan_seconds_se,
    speedup = consensus_scan_seconds / inflect_scan_seconds,
    speedup_amortized = amortized_seconds / inflect_scan_seconds,
    speedup_legacy_inflect = legacy_inflect_scan_seconds / inflect_scan_seconds
  ),
  machine = list(
    cores = parallel::detectCores(),
    sysname = Sys.info()[["sysname"]],
    r_version = R.version.string,
    inflect_version = cache_provenance$package_versions$fastINFLECT,
    flowsom_version = cache_provenance$package_versions$FlowSOM,
    diptest_version = cache_provenance$package_versions$diptest,
    consensusclusterplus_version = cache_provenance$package_versions$ConsensusClusterPlus
  )
)

out_path <- existing_cache_path
saveRDS(cache, file = out_path, compress = "xz")
message(sprintf("Saved %s (%.0f KB)", out_path, file.size(out_path) / 1024))

# ── 10. Sanity checks ────────────────────────────────────────────────────────
stopifnot(
  nrow(consensus_df) == length(k_range),
  all(c("k", "unimodality", "method") %in% names(cache$comparison_df)),
  is.finite(cache$totals$speedup),
  cache$totals$inflect_scan_seconds > 0,
  cache$totals$speedup_legacy_inflect > 1,
  all(c("method", "k", "cluster", "marker", "expression", "qc_pass") %in%
        names(cache$marker_histogram))
)
message(sprintf(paste0(
  "fastINFLECT scan: %.2fs | original INFLECT estimated: %.1fs | ",
  "consensus scan: %.1fs | speedups: original ~%.0fx, consensus ~%.0fx"
),
inflect_scan_seconds, legacy_inflect_scan_seconds, consensus_scan_seconds,
cache$totals$speedup_legacy_inflect, cache$totals$speedup))
message("All sanity checks passed.")
