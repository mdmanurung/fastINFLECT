## POC: prove that memoizing QC per distinct dendrogram subtree (keyed by its
## sorted SOM-node set) reproduces the historical FlowSOMQC per-k score EXACTLY,
## and measure the speedup vs the current per-k recomputation.
##
## Run with the R4_51 conda env.
suppressPackageStartupMessages({library(fastINFLECT); library(FlowSOM); library(diptest)})

load(system.file("extdata", "Levine32sample.Rdata", package = "fastINFLECT"))
codes <- dataset$map$codes
nNodes <- nrow(codes)
mapping <- dataset$map$mapping[, 1]

## ---- prepare data exactly as FlowSOMQC does ------------------------------
data <- dataset$data
if (isTRUE(dataset$scale)) {
  for (j in seq_len(ncol(data)))
    data[, j] <- data[, j] * dataset$scaled.scale[j] + dataset$scaled.center[j]
}
colnames(data) <- dataset$prettyColnames
clustering.markers <- dataset$prettyColnames[dataset$map$colsUsed]
ordered.markers <- gtools::mixedsort(clustering.markers)   # only.clustering.markers=TRUE
nMarkers <- length(ordered.markers)

## events per SOM node (fixed, independent of k)
node_events <- split(seq_len(nrow(data)), mapping)         # names are node ids

th.pvalue <- 0.05; th.IQR <- 2

## accuracy row for one metacluster = union of member nodes' events.
## Historical 1.0 equivalence only: uniform.test = "both", zeroes.in = FALSE.
## This is not modality evidence under the 2.0 interpretation contract.
acc_row <- function(member_nodes) {
  rows <- unlist(node_events[as.character(member_nodes)], use.names = FALSE)
  out <- rep(NA, nMarkers); names(out) <- ordered.markers
  if (length(rows) > 1) {
    ex <- data[rows, ordered.markers, drop = FALSE]
    for (m in ordered.markers) {
      v <- ex[, m]; me <- v[v > 0]
      if (length(me) < 5) me <- c(rep(0, 5 - length(me)), me)
      p  <- diptest::dip.test(me)$p.value
      q  <- stats::quantile(me)
      out[m] <- (p >= th.pvalue) && ((q[4] - q[2]) < th.IQR)
    }
  }
  out
}

k_range <- 5:20
fit <- stats::hclust(stats::dist(codes, method = "minkowski"), method = "ward.D2")
ct  <- stats::cutree(fit, k = k_range)                     # nNodes x length(k_range)

## ---- ground truth: current FlowSOMQC path, per k -------------------------
t_ref <- system.time({
  ref_scores <- sapply(seq_along(k_range), function(i) {
    m <- FlowSOMQC(dataset, ct[, i], zeroes.in = FALSE, verbose = FALSE)
    sum(m, na.rm = TRUE) * 100 / prod(dim(m))
  })
})[["elapsed"]]

## ---- memoized: compute acc_row once per unique node-set ------------------
cache <- new.env(parent = emptyenv())
t_memo <- system.time({
  memo_scores <- sapply(seq_along(k_range), function(i) {
    lab <- ct[, i]
    clusters <- sort(unique(lab))
    pass <- 0
    for (cl in clusters) {
      nodes <- which(lab == cl)
      key <- paste0(nodes, collapse = ",")          # sorted (which() is increasing)
      row <- cache[[key]]
      if (is.null(row)) { row <- acc_row(nodes); assign(key, row, envir = cache) }
      pass <- pass + sum(row, na.rm = TRUE)
    }
    pass * 100 / (length(clusters) * nMarkers)
  })
})[["elapsed"]]

cat(sprintf("\nk range: %d..%d\n", min(k_range), max(k_range)))
cat(sprintf("distinct subtrees evaluated (cached): %d  (theoretical max 2*nNodes-1=%d)\n",
            length(ls(cache)), 2 * nNodes - 1))
cat(sprintf("reference (per-k FlowSOMQC) : %6.1f s\n", t_ref))
cat(sprintf("memoized (per-subtree once) : %6.1f s   -> %.1fx faster\n", t_memo, t_ref / t_memo))
cat("\nmax |memo - ref| score diff: ",
    format(max(abs(memo_scores - ref_scores)), digits = 3), "\n")
cat("EXACT MATCH: ", isTRUE(all.equal(memo_scores, ref_scores, tolerance = 1e-9)), "\n")
print(data.frame(k = k_range, ref = round(ref_scores, 4), memo = round(memo_scores, 4)))
