## data-raw/make-vignette-cache.R
##
## Run this script ONCE with the R4_51 (or scale_r) conda environment to
## generate inst/extdata/vignette_cache.rds. The vignette loads this cache
## instead of re-running the heavy computations at build time.
##
## Usage (from the package root):
##   /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
##     data-raw/make-vignette-cache.R

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(fastINFLECT)
  }
  library(FlowSOM)
  library(ggplot2)
})

set.seed(42)

existing_cache_path <- file.path("inst", "extdata", "vignette_cache.rds")
existing_cache <- if (file.exists(existing_cache_path)) readRDS(existing_cache_path) else NULL
refresh_cache <- identical(Sys.getenv("FASTINFLECT_REFRESH_VIGNETTE_CACHE"), "1")

# ── 1. Load the bundled dataset ──────────────────────────────────────────────
message("Loading Levine32sample dataset...")
data_path <- system.file("extdata", "Levine32sample.Rdata", package = "fastINFLECT")
stopifnot(nchar(data_path) > 0)
load(data_path)   # creates object 'dataset' (class FlowSOM)

codes <- dataset$map$codes          # 375 x 32 SOM node prototypes
n_nodes <- dataset$map$nNodes       # 375

# ── 2. Unimodality scoring helper ─────────────────────────────────────────────
# Mirrors the exact aggregation in iteration.QC (R/iteration-QC.R line 55):
#   sum(accuracy.matrix, na.rm = TRUE) * 100 / prod(dim(accuracy.matrix))
uni_score <- function(mc) {
  m <- FlowSOMQC(
    FlowSOM.results = dataset,
    metaclustering   = mc,
    zeroes.in        = FALSE,
    verbose          = FALSE
  )
  sum(m, na.rm = TRUE) * 100 / prod(dim(m))
}

cluster_quality <- function(method, k, metaclustering, selection) {
  accuracy <- FlowSOMQC(
    FlowSOM.results = dataset,
    metaclustering = as.integer(metaclustering),
    zeroes.in = FALSE,
    verbose = FALSE
  )
  n_markers <- rowSums(!is.na(accuracy))
  passed_markers <- rowSums(accuracy, na.rm = TRUE)
  data.frame(
    method = method,
    selection = selection,
    k = as.integer(k),
    cluster = seq_len(nrow(accuracy)),
    passed_markers = passed_markers,
    n_markers = n_markers,
    cluster_unimodality = passed_markers * 100 / n_markers,
    stringsAsFactors = FALSE
  )
}

summarise_cluster_quality <- function(cluster_quality_df) {
  pieces <- split(
    cluster_quality_df,
    factor(cluster_quality_df$method, levels = unique(cluster_quality_df$method))
  )
  do.call(rbind, lapply(pieces, function(d) {
    data.frame(
      method = unique(d$method),
      selection = unique(d$selection),
      k = unique(d$k),
      clusters = nrow(d),
      median_cluster_unimodality = stats::median(d$cluster_unimodality),
      min_cluster_unimodality = min(d$cluster_unimodality),
      clusters_below_95 = sum(d$cluster_unimodality < 95),
      clusters_below_100 = sum(d$cluster_unimodality < 100),
      stringsAsFactors = FALSE
    )
  }))
}

bin_cluster_quality <- function(cluster_quality_df) {
  breaks <- seq(0, 100, by = 10)
  labels <- paste0(head(breaks, -1), "-", tail(breaks, -1), "%")
  bins <- cluster_quality_df
  bins$quality_bin <- cut(
    bins$cluster_unimodality,
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE,
    right = TRUE
  )
  pieces <- split(bins, factor(bins$method, levels = unique(bins$method)))
  do.call(rbind, lapply(pieces, function(d) {
    tab <- as.data.frame(table(d$quality_bin), stringsAsFactors = FALSE)
    names(tab) <- c("quality_bin", "cluster_count")
    tab$method <- unique(d$method)
    tab$k <- unique(d$k)
    tab$quality_midpoint <- breaks[-1] - 5
    tab$cluster_fraction <- tab$cluster_count / sum(tab$cluster_count)
    tab[, c("method", "k", "quality_bin", "quality_midpoint",
            "cluster_count", "cluster_fraction")]
  }))
}

# ── 3. Sweep k = 5:25 for consensus and hierarchical ─────────────────────────
k_range <- 5:25
if (!is.null(existing_cache) && !refresh_cache &&
    all(c("comparison_df", "autok_point", "inflect", "mp_plot") %in% names(existing_cache))) {
  message("Reusing existing vignette sweep cache. Set FASTINFLECT_REFRESH_VIGNETTE_CACHE=1 to retime.")
  comparison_df <- existing_cache$comparison_df
  autok_point <- existing_cache$autok_point
  inflect_res <- existing_cache$inflect
  mp <- list(plot = existing_cache$mp_plot)
} else {
  message(sprintf("Sweeping k = %d:%d for consensus and hierarchical clustering...",
                  min(k_range), max(k_range)))

  rows <- vector("list", length(k_range) * 2L)
  idx  <- 1L

  for (k in k_range) {
    message(sprintf("  k = %d", k))

    # FlowSOM consensus metaclustering
    mc_consensus <- metaClustering_consensus(codes, k = k, seed = 42)
    rows[[idx]] <- data.frame(
      k           = k,
      method      = "FlowSOM consensus",
      unimodality = uni_score(mc_consensus),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L

    # Ward.D2 hierarchical (fastINFLECT's internal engine)
    mc_hclust <- metaClusteringhclust(codes, nClus = k)
    rows[[idx]] <- data.frame(
      k           = k,
      method      = "Hierarchical (ward.D2)",
      unimodality = uni_score(mc_hclust),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  comparison_df <- do.call(rbind, rows)

  # ── 4. FlowSOM automatic-k metaclustering ───────────────────────────────────
  message("Running FlowSOM MetaClustering (auto-k)...")
  mc_autok     <- MetaClustering(codes,
                                 method = "metaClustering_consensus",
                                 max    = max(k_range),
                                 seed   = 42)
  autok_k      <- length(unique(mc_autok))
  autok_score  <- uni_score(mc_autok)
  autok_point  <- data.frame(
    k           = autok_k,
    method      = "FlowSOM auto-k",
    unimodality = autok_score,
    stringsAsFactors = FALSE
  )
  message(sprintf("  auto-k chose k = %d (unimodality = %.1f%%)", autok_k, autok_score))

  # ── 5. Run fastINFLECT ─────────────────────────────────────────────────────
  message("Running fastINFLECT (set.i = 5:25, multicore = FALSE)...")
  t0 <- proc.time()
  inflect_res <- INFLECT(
    FlowSOM.results = dataset,
    set.i           = k_range,
    multicore       = FALSE,
    zeroes.in       = FALSE,
    verbose         = FALSE
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]
  message(sprintf("  fastINFLECT complete in %.0f s; knee = %d",
                  elapsed, inflect_res$lfunction$knee))

  mp <- marker.performance(inflect_res, ggtitle = "Levine32 - marker performance")
}

# ── 6. Cluster-level quality histograms ──────────────────────────────────────
message("Computing all-cluster unimodality distributions for selected methods...")
fixed_k <- 10L
mc_consensus_fixed <- metaClustering_consensus(codes, k = fixed_k, seed = 42)
mc_hclust_fixed <- metaClusteringhclust(codes, nClus = fixed_k)
mc_autok <- MetaClustering(codes, method = "metaClustering_consensus",
                           max = max(k_range), seed = 42)
autok_k <- length(unique(mc_autok))
knee_k <- as.integer(round(inflect_res$lfunction$knee))
threshold_k <- inflect_res$selection$k[inflect_res$selection$method == "threshold"]
threshold_k <- if (length(threshold_k) && !is.na(threshold_k)) {
  as.integer(round(threshold_k[[1L]]))
} else {
  knee_k
}
cluster_quality_df <- do.call(rbind, list(
  cluster_quality("FlowSOM consensus", fixed_k, mc_consensus_fixed, "fixed k = 10"),
  cluster_quality("FlowSOM auto-k", autok_k, mc_autok, "automatic k"),
  cluster_quality("Hierarchical (ward.D2)", fixed_k, mc_hclust_fixed, "fixed k = 10"),
  cluster_quality("fastINFLECT inflection", knee_k,
                  inflect_res$metaclustering.list[[as.character(knee_k)]],
                  "automatic inflection"),
  cluster_quality("fastINFLECT threshold", threshold_k,
                  inflect_res$metaclustering.list[[as.character(threshold_k)]],
                  "automatic threshold")
))
cluster_quality_summary <- summarise_cluster_quality(cluster_quality_df)
cluster_quality_bins <- bin_cluster_quality(cluster_quality_df)

# ── 7. Save the cache ─────────────────────────────────────────────────────────
cache <- list(
  comparison_df = comparison_df,
  autok_point   = autok_point,
  inflect       = inflect_res,
  mp_plot       = mp$plot,
  cluster_quality_df = cluster_quality_df,
  cluster_quality_summary = cluster_quality_summary,
  cluster_quality_bins = cluster_quality_bins,
  k_range       = k_range,
  session_info  = sessionInfo()
)

out_path <- existing_cache_path
saveRDS(cache, file = out_path, compress = "xz")

size_kb <- file.size(out_path) / 1024
message(sprintf("Cache saved to %s (%.0f KB)", out_path, size_kb))

# ── 8. Quick sanity checks ────────────────────────────────────────────────────
stopifnot(
  nrow(comparison_df) == length(k_range) * 2L,
  all(c("k", "method", "unimodality") %in% names(comparison_df)),
  is.numeric(autok_point$k),
  inherits(inflect_res, "inflect.results"),
  is.numeric(inflect_res$lfunction$knee),
  inflect_res$lfunction$knee >= min(k_range),
  inflect_res$lfunction$knee <= max(k_range),
  all(c("method", "k", "cluster", "cluster_unimodality") %in% names(cluster_quality_df)),
  all(cluster_quality_df$cluster_unimodality >= 0),
  all(cluster_quality_df$cluster_unimodality <= 100),
  all(c("method", "quality_bin", "cluster_fraction") %in% names(cluster_quality_bins))
)
message("All sanity checks passed.")
