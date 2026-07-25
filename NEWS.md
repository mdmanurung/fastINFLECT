# fastINFLECT 1.0.0

fastINFLECT is the renamed, faster implementation of the original INFLECT
method. The exported `INFLECT()` function remains available for API
continuity.

## Package interface

* `set.i` is required and literal. It must contain at least five unique,
  strictly increasing integer cluster counts within the SOM-node range.
  Two-value inputs such as `c(25, 100)` fail because their spacing is
  ambiguous. Use `seq.int(25L, 100L, by = 5L)` for that exact scan, or
  `inflect_adaptive_set_i(n_nodes, max_k)` for an explicitly bounded adaptive
  schedule.
* Serial QC is the default. Parallel QC requires `multicore = TRUE` and at
  least two effective workers. Worker and subtree results are validated before
  result matrices are assembled. Platforms without fork support use a recorded
  serial fallback after the same validation.
* `zeroes.in = TRUE` retains negative, zero, and positive finite transformed
  values. Setting `zeroes.in = FALSE` excludes every non-positive value, warns
  when negatives are present, and records per-marker exclusion counts.
* Canonical output uses `scores` and `qc_pass_rate`. Separate `dip_pass`,
  `iqr_pass`, `combined_pass`, and `criterion_pass` matrices retain the
  component evidence. The historical `U.set`, `collection.U`, `Unimodality`,
  `Accuracy.sets`, and `Accuracy.matrixes` names remain as deprecated
  compatibility aliases.

## Performance and memory

* `iteration.QC()` evaluates each distinct SOM-node set once. A nested
  `hclust` and `cutree` sweep contains at most `2 * nNodes - 1` distinct
  clusters, regardless of schedule density.
* On the bundled Levine32 SOM with 375 nodes, 32,288 events, and 32 markers, a
  historical engineering benchmark completed the dense `k = 5:370` sweep in
  seconds. The projected repeated per-k implementation required hours on the
  same machine.
* The dip-test kernel uses `diptest::dip()` with legacy-matched p-value
  interpolation and an exact small-sample fallback. The sampled legacy
  benchmark measured `diptest::dip.test()` at about 12 times the per-call cost.
* Rcpp accelerates the IQR-spread test and p-value interpolation. Validated
  pure-R paths remain available when compiled code is unavailable.
* Subtrees are scored by row index, one marker at a time, instead of copying a
  full events-by-markers matrix for every subtree.
* `max.events.per.node` can cap event retention for sensitivity analyses. Each
  SOM node is sampled once with `seed`, and the same retained rows feed every
  marker, criterion, metaclustering, and worker. Provenance records original
  and retained event counts.
* Fork-based parallel scoring operates over distinct subtrees on Unix-like
  systems. This replaces the repeated per-k `foreach` and `doParallel`
  implementation, and those packages are no longer dependencies.

## Statistical evidence and reproducibility

* `qc.details` retains dip p-values, IQR values, original and testable event
  counts, non-finite and non-positive exclusions, and explicit failure reasons.
  Both component criteria remain available regardless of `uniform.test`.
* Simulated dip p-values and optional event samples use overflow-safe seeds
  derived for each subtree and marker. The implementation preserves the caller
  RNG state and returns identical tested results with serial and two-worker
  execution.
* Small-sample dip interpolation uses `ties = mean`, which avoids
  duplicate-grid interpolation warnings without changing the tested p-values.
* `max.n.diptest` provides a deterministic marker-level sample cap for
  sensitivity analysis of dip-test sample-size effects. The default `NULL`
  retains all testable events.
* `bimodality.coefficient()` provides Sarle's bimodality coefficient as an
  additional diagnostic. Its documentation notes that skewed unimodal
  distributions can be flagged.
* Provenance records the selected criterion, thresholds, marker panel, zero
  handling, sampling order and seeds, effective workers, model layer weights,
  retained events, and per-criterion summaries.

## Selecting k

* `INFLECT()` reports three candidate endpoints in `selection`:
  * `inflection`, the original LL.4 log-logistic knee;
  * `kneedle`, a nonparametric knee from `inflect_kneedle()` following
    Satopää et al. (2011);
  * `threshold`, the smallest tested k that reaches the requested
    criterion-specific pass rate through `inflect_threshold_k()`, with a
    default target of 0.95.
* Each recommendation reports whether k was directly tested, whether a
  partition is available, and whether its score came from tested data or the
  fitted curve. An unscheduled fitted k is labelled
  `fitted_estimate_no_partition` and is not evaluated implicitly.

## Real-model validation

* The checkpointed, Slurm-array-compatible workflow in
  `data-raw/validate-real-model-modality.R` materialises Ward and FlowSOM
  controls, uses deterministic nested samples, and combines the dip test with
  the optional `multimode` ACR test. Base and sample-size refinement evidence
  is stored under `inst/benchmarks/`.
* The frozen BMV model contains 39,050,953 events, 900 SOM nodes, and 27
  markers. The full-event comparison gives a combined INFLECT inflection at
  k = 27, a combined and dip Kneedle at k = 44, a consensus plateau beginning
  at k = 29, and an upper spread sensitivity at k = 55.
* The operational BMV partition is nominal k = 44, with 43 event-populated
  clusters and one empty nominal cluster. Its independent audit contains 56
  detected, 118 ambiguous, 28 unresolved, and 986 no-detected cluster-marker
  entries. These results do not support a global unimodality claim.
* Documentation uses `detected_multimodality`, `ambiguous`,
  `no_detected_multimodality`, and `unresolved_discrete`. A criterion pass or
  failure to reject the unimodal null is not treated as proof of unimodality.

## Documentation

* The vignette "Using fastINFLECT 1.0.0" defines the literal schedule, value
  handling, criterion-specific outputs, parallel execution, sampling controls,
  selection metadata, and provenance contract.
* "Benchmarking fastINFLECT against consensus metaclustering" separates the
  historical Levine32 performance cache from the full BMV selection and
  modality evidence. The historical cache is retained for engineering
  context, not modality claims.
* "Comparing metaclustering strategies with fastINFLECT" covers Ward.D2 and
  FlowSOM consensus partitions, common marker panels, cluster-size-aware
  summaries, and independent modality validation.

# INFLECT 0.2.1

* Added the `inflect.results` S3 class with `print`, `summary`, `plot`, and
  `as.data.frame` methods.
* Added kohonen SOM support through the `as_inflect_som` adapter.
* Added provenance tracking to `INFLECT()` output.
* Added the original metaclustering comparison vignette for the bundled
  Levine32 CyTOF dataset.
