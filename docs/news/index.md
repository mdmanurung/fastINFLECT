# Changelog

## fastINFLECT 2.0.0

### Breaking API changes

  - `set.i` is now required and literal. It must contain at least five
    unique, strictly increasing integer cluster counts within the
    SOM-node range. Two-value inputs such as `c(25, 100)` now fail
    instead of silently expanding past 100. Use `seq.int(25L, 100L, by
    = 5L)` for that exact scan, or `inflect_adaptive_set_i(n_nodes,
    max_k)` for the former adaptive spacing with an explicit upper
    bound.
  - Parallel QC is now opt-in: `multicore = FALSE` is the default.
    Enabling it requires `cores >= 2`; failed worker/subtree results are
    validated and reported before result matrices are assembled.
    Platforms without fork support use a provenance-recorded serial
    fallback after applying the same validation.
  - `zeroes.in = TRUE` is now the scientific default. It retains
    negative, zero, and positive finite transformed values. The
    compatibility setting `zeroes.in = FALSE` excludes every
    non-positive value, warns when negatives are present, and records
    per-marker exclusion counts.
  - Aggregate output is now named `qc_pass_rate` and carries
    criterion-aware labels. Canonical output includes separate
    `dip_pass`, `iqr_pass`, `combined_pass`, and `criterion_pass`
    matrices. The historical `U.set`, `collection.U`, `Unimodality`,
    `Accuracy.sets`, and `Accuracy.matrixes` names remain as documented
    deprecated aliases for migration.

### Memory and statistical controls

  - Subtrees are scored by row index, one marker at a time, instead of
    copying a full events-by-markers matrix for every subtree.
  - Added `max.events.per.node`. When set, every SOM node is sampled
    once using `seed`, and that shared event sample is used for both dip
    and IQR scoring in every worker. Original and retained event counts
    are recorded in provenance.
  - `qc.details` now retains dip p-values, IQR values, original and
    testable event counts, non-finite and non-positive exclusions, and
    explicit failure reasons. Both component criteria are retained
    regardless of the aggregate selected by `uniform.test`.
  - Simulated dip p-values and optional subsampling use overflow-safe
    seeds derived per subtree and marker. Results preserve the caller
    RNG state and are identical across serial and two-worker execution.
  - Small-sample dip interpolation now passes `ties = mean` explicitly,
    avoiding duplicate-grid warnings without changing p-values.
  - Recommendation tables now identify whether each k was directly
    tested and report the fitted score for fitted-curve recommendations.
    An unscheduled fitted k is labelled `fitted_estimate_no_partition`
    and is never silently evaluated.
  - QC-boundary provenance now consolidates the selected criterion,
    thresholds, marker panel, zero handling, sampling order and seeds,
    effective workers, model layer weights, retained events, and
    per-criterion summaries.

### Validation and interpretation

  - Added a checkpointed, Slurm-array-compatible real-model modality
    workflow in `data-raw/validate-real-model-modality.R`. It
    materialises the requested Ward and FlowSOM controls, uses
    deterministic nested samples, combines the dip test with the
    optional `multimode` ACR test, and records base and sample-size
    refinement evidence under `inst/benchmarks/`.
  - Documentation now uses `detected_multimodality`, `ambiguous`,
    `no_detected_multimodality`, and `unresolved_discrete`. A criterion
    pass or a failure to reject unimodality is not described as proof
    that a cluster is truly unimodal.

## fastINFLECT 1.0.0

  - Package renamed from `INFLECT` to `fastINFLECT`. The exported
    `INFLECT()` function is retained for API continuity, and the
    documentation now explicitly frames this repository as a fast
    reimplementation of the original INFLECT method.

### Performance

  - Rewrote the QC sweep to be dramatically faster while preserving the
    original QC decisions. On the bundled Levine32 SOM (375 nodes,
    32,288 events, 32 markers) a full dense sweep of every k from 5 to
    370 now runs in a few seconds; the previous per-k implementation
    would have taken hours.
      - **Sub-tree memoisation.** `iteration.QC()` now evaluates each
        distinct SOM-node set only once. Because `hclust`/`cutree`
        metaclusterings are nested, there are at most `2 * nNodes - 1`
        distinct clusters across the whole sweep regardless of how
        densely `set.i` samples k, so the cost no longer grows with the
        number of tested k values.
      - **Faster kernel.** The dip test now uses `diptest::dip()` for
        the statistic plus legacy-matched p-value interpolation, with
        exact small-sample fallback, instead of the \~12x more expensive
        `dip.test()` per cell.
      - **Compiled accelerators (Rcpp).** The inter-quartile-range
        spread test and the hot-path p-value interpolation are
        implemented in C++ (`src/inflect.cpp`), with the pure-R paths
        retained as validated fallbacks when the package is used without
        its compiled code.
      - Parallelism in `iteration.QC()` now fans out over distinct
        sub-trees via fork-based `parallel::mclapply()` (Unix),
        replacing the previous per-k `foreach`/`doParallel` loop.
        `doParallel` and `foreach` are no longer dependencies.

### Choosing k

  - `INFLECT()` now reports the recommended number of metaclusters under
    three criteria in a new `selection` element (also shown by `print()`
    and `summary()`):
      - `inflection` — the original LL.4 log-logistic knee (unchanged
        default);
      - `kneedle` — a nonparametric knee (`inflect_kneedle()`, Satopää
        et al. 2011);
      - `threshold` — the smallest k reaching a target unimodality
        (`inflect_threshold_k()`, controlled by the new `target`
        argument, default 0.95). This directly encodes “admit no cluster
        with residual bimodal marker expression, without
        over-clustering”.

### Statistics

  - New `max.n.diptest` argument (in `INFLECT()` and `iteration.QC()`)
    optionally caps the number of events used per (cluster, marker) dip
    test. This removes the sample-size sensitivity of the dip test —
    large clusters otherwise reject unimodality for negligible
    deviations — making the score size-robust. Subsampling is seeded
    (`seed`) and therefore deterministic. The default `NULL` preserves
    the legacy scoring path.
  - New exported `bimodality.coefficient()` computes Sarle’s bimodality
    coefficient as a fast, size-robust diagnostic that complements the
    dip test (documented with its caveat: it can over-flag
    skewed-but-unimodal markers).

### Documentation

  - New vignette “Benchmarking fastINFLECT against consensus
    metaclustering” times the cost of scanning k on the bundled Levine32
    SOM: fastINFLECT’s memoised sweep versus the original INFLECT loop
    and FlowSOM consensus metaclustering, both as a single
    ConsensusClusterPlus run and as the per-k convenience pattern. It
    also includes marker-expression histograms showing residual
    bimodality in the standard consensus partition versus overlapping
    fastINFLECT unimodal clusters. (Timings and histogram data are
    precomputed by `data-raw/make-benchmark-cache.R`.)
  - The metaclustering comparison vignette now includes all-cluster
    histograms on the Levine32 dataset, showing the distribution of
    per-metacluster marker-unimodality pass rates for FlowSOM consensus,
    FlowSOM auto-k, hierarchical cuts, and fastINFLECT recommendations.
