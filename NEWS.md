# fastINFLECT 1.0.0

* Package renamed from `INFLECT` to `fastINFLECT`. The exported `INFLECT()`
  function is retained for API continuity, and the documentation now explicitly
  frames this repository as a fast reimplementation of the original INFLECT
  method.

## Performance

* Rewrote the QC sweep to be dramatically faster while producing byte-identical
  results to previous versions. On the bundled Levine32 SOM (375 nodes, 32,288
  events, 32 markers) a full dense sweep of every k from 5 to 370 now runs in a
  few seconds; the previous per-k implementation would have taken hours.
  * **Sub-tree memoisation.** `iteration.QC()` now evaluates each distinct SOM-node
    set only once. Because `hclust`/`cutree` metaclusterings are nested, there are
    at most `2 * nNodes - 1` distinct clusters across the whole sweep regardless of
    how densely `set.i` samples k, so the cost no longer grows with the number of
    tested k values.
  * **Faster kernel.** The dip test now uses `diptest::dip()` for the statistic
    plus a direct, bit-exact reproduction of `diptest::dip.test()`'s p-value table
    interpolation, instead of the ~12x more expensive `dip.test()` per cell.
  * **Compiled accelerators (Rcpp).** The inter-quartile-range spread test and the
    p-value table interpolation are implemented in C++ (`src/inflect.cpp`), with
    the pure-R paths retained as validated, bit-identical fallbacks when the
    package is used without its compiled code.
  * Parallelism in `iteration.QC()` now fans out over distinct sub-trees via
    fork-based `parallel::mclapply()` (Unix), replacing the previous per-k
    `foreach`/`doParallel` loop. `doParallel` and `foreach` are no longer
    dependencies.

## Choosing k

* `INFLECT()` now reports the recommended number of metaclusters under three
  criteria in a new `selection` element (also shown by `print()` and `summary()`):
  * `inflection` — the original LL.4 log-logistic knee (unchanged default);
  * `kneedle` — a nonparametric knee (`inflect_kneedle()`, Satopää et al. 2011);
  * `threshold` — the smallest k reaching a target unimodality
    (`inflect_threshold_k()`, controlled by the new `target` argument, default
    0.95). This directly encodes "admit no cluster with residual bimodal marker
    expression, without over-clustering".

## Statistics

* New `max.n.diptest` argument (in `INFLECT()` and `iteration.QC()`) optionally
  caps the number of events used per (cluster, marker) dip test. This removes the
  sample-size sensitivity of the dip test — large clusters otherwise reject
  unimodality for negligible deviations — making the score size-robust. Subsampling
  is seeded (`seed`) and therefore deterministic. The default `NULL` reproduces the
  legacy score exactly.
* New exported `bimodality.coefficient()` computes Sarle's bimodality coefficient
  as a fast, size-robust diagnostic that complements the dip test (documented with
  its caveat: it can over-flag skewed-but-unimodal markers).

## Documentation

* New vignette "Benchmarking fastINFLECT against consensus metaclustering" times the
  cost of scanning k on the bundled Levine32 SOM: fastINFLECT's memoised sweep versus
  the original INFLECT loop and FlowSOM consensus metaclustering, both as a single
  ConsensusClusterPlus run and as the per-k convenience pattern. It also includes
  marker-expression histograms showing residual bimodality in the standard
  consensus partition versus overlapping fastINFLECT unimodal clusters. (Timings and
  histogram data are precomputed by `data-raw/make-benchmark-cache.R`.)
* The metaclustering comparison vignette now includes all-cluster histograms on
  the Levine32 dataset, showing the distribution of per-metacluster
  marker-unimodality pass rates for FlowSOM consensus, FlowSOM auto-k,
  hierarchical cuts, and fastINFLECT recommendations.

# INFLECT 0.2.1

* Added `inflect.results` S3 class with `print`, `summary`, `plot`, and `as.data.frame` methods.
* Added kohonen SOM support via `as_inflect_som` adapter.
* Added provenance tracking to `INFLECT()` output.
* Added vignette "Comparing metaclustering strategies with INFLECT" demonstrating the package
  on the bundled Levine32 CyTOF dataset.
