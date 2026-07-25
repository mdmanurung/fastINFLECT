# fastINFLECT

fastINFLECT performs marker-level quality control across literal
metaclustering schedules for FlowSOM and kohonen self-organising maps.
It retains separate Hartigan dip-test and IQR-spread evidence, fits a
diagnostic pass-rate curve, and reports candidate endpoints without
pretending that an unscheduled fitted value has a partition.

A dip or IQR pass is a criterion-specific screening result. It is not
proof that a marker distribution, cluster, or biological population is
truly unimodal.

## Installation

fastINFLECT is not on CRAN. The simplest installation path uses `pak`,
which also resolves the Bioconductor dependency on FlowSOM:

``` r
# install.packages("pak")
pak::pak("mdmanurung/fastINFLECT")
```

For kohonen SOM objects, install the optional `kohonen` package:

``` r
install.packages("kohonen")
```

The independent real-model modality audit additionally uses the optional
`multimode` package. It is not required for ordinary fastINFLECT use.

## Quick start

The package ships a downsampled Levine32 FlowSOM object:

``` r
library(fastINFLECT)

load(system.file(
  "extdata",
  "Levine32sample.Rdata",
  package = "fastINFLECT"
))

schedule <- 5:20
result <- INFLECT(
  FlowSOM.results = dataset,
  set.i = schedule,
  uniform.test = "both",
  multicore = FALSE
)
```

Version 2.0 retains all finite transformed values, including negative
values and zero, by default. The aggregate score is named for what it
is:

``` r
result$scores
# k, qc_pass_rate, criterion

result$dip_pass[["10"]]
result$iqr_pass[["10"]]
result$combined_pass[["10"]]

details <- result$qc.details[["10"]]
details$dip_p_value
details$iqr
details$failure_reason
```

The diagnostic curve is available as `result$ggplot`. Its selection
table distinguishes tested partitions from fitted estimates:

``` r
result$selection[c(
  "method",
  "k",
  "qc_pass_rate_at_k",
  "partition_available",
  "k_status"
)]
```

Only rows with `partition_available = TRUE` have a corresponding entry
in `result$metaclustering.list`. A fitted inflection with `k_status =
"fitted_estimate_no_partition"` remains an estimate until that literal k
is explicitly evaluated.

The pass-rate curve alone should not be used to declare a selected
partition unimodal. Inspect the separate criterion matrices and, for
consequential claims, validate complete transformed marker distributions
with independent tests and sensitivity analyses.

## Migrating to 2.0.0

`set.i` is mandatory, literal, unique, and strictly increasing. This
call evaluates exactly the printed values and never scans above 100:

``` r
schedule <- seq.int(25L, 100L, by = 5L)
print(schedule)
result <- INFLECT(dataset, set.i = schedule)
```

To request the former adaptive spacing, construct it explicitly with a
hard upper bound:

``` r
schedule <- inflect_adaptive_set_i(
  n_nodes = dataset$map$nNodes,
  max_k = 100L
)
```

Parallel scoring is opt-in and requires an explicit valid worker count:

``` r
result <- INFLECT(
  dataset,
  set.i = schedule,
  multicore = TRUE,
  cores = 2L,
  seed = 42L
)
```

On platforms without fork support, the same arguments are validated and
a recorded serial fallback is used. Simulated dip p-values and optional
subsampling are deterministic per subtree and marker, preserve the
caller RNG state, and are invariant to worker count.

`zeroes.in = FALSE` is now an explicit compatibility mode. It excludes
every non-positive value, warns when negative values are present, and
records per-marker exclusion counts:

``` r
result$provenance$zero_handling
```

Canonical 2.0 names are `scores`, `qc_pass_rate`, `dip_pass`,
`iqr_pass`, `combined_pass`, and `criterion_pass`. `collection.U`,
`U.set`, `Unimodality`, `Accuracy.sets`, and `Accuracy.matrixes` are
deprecated aliases retained for migration.

See `vignette("migrating-to-fastINFLECT-2")` for the complete contract.

## Real-model modality audit

The repository includes a resumable audit for the 39,050,953-event BMV
model:

``` bash
data-raw/submit-real-model-modality.sh
```

It materialises every requested Ward and FlowSOM control partition,
stages nested deterministic event samples, runs dip and ACR tests with
Benjamini-Hochberg correction within each 27-marker cluster family,
repeats flagged and audited pairs at 1,000 and 5,000 events, and writes
checkpointed evidence under `inst/benchmarks/real-model-modality/`.

For the frozen BMV model (`md5 784b3a6048c379ca64bbe8d4c0fb063d`), the
literal full-event k = 25,…,100 comparison gives a combined INFLECT
inflection at k = 27 and combined/dip Kneedle at k = 44. The first
five-k consensus-stability plateau starts at k = 29, while IQR and
event-weighted dispersion place the upper spread sensitivity at k = 55.
Accordingly, the fitted-model operational recommendation is nominal k =
44, with k = 27–29 as a lower-resolution structural sensitivity and k =
55 as a spread sensitivity.

Nominal k = 44 contains 43 event-populated clusters because one SOM-only
cluster has zero events. Its independent final audit contains 56
detected, 118 ambiguous, 28 unresolved, and 986 no-detected
cluster-marker entries. This supports neither a global unimodality claim
nor biological optimality. Exact comparisons and provenance are in
`inst/benchmarks/real-model-selection/`.

The allowed outcome labels are `detected_multimodality`, `ambiguous`,
`no_detected_multimodality`, and `unresolved_discrete`. The workflow
never jitters tied measurements and never converts “no detected
multimodality” into a claim that a cluster is truly unimodal.

## Relationship to the original INFLECT package

fastINFLECT is derived from the INFLECT implementation developed by Jan
Verhoeff in the lab of JJ. Garcia-Vallejo and available from the
[GarciaVallejoLab
repository](https://github.com/jnverhoeff/GarciaVallejoLab). The
exported `INFLECT()` name is retained for API continuity. Memoised
hierarchy cuts, marker-at-a-time indexed scoring, deterministic
sampling, and compiled accelerators provide the faster implementation.
