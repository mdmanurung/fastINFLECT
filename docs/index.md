# fastINFLECT

fastINFLECT performs marker-level quality control across literal
metaclustering schedules for FlowSOM and kohonen self-organising maps.
It retains separate Hartigan dip-test and IQR-spread evidence, fits a
diagnostic pass-rate curve, and labels unscheduled fitted endpoints as
estimates without partitions.

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

The real-model modality audit uses the optional `multimode` package.
Ordinary fastINFLECT use does not require it.

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

fastINFLECT retains all finite transformed values, including negative
values and zero, by default. The aggregate score is named for what it
measures:

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

## Package interface

`set.i` is mandatory, literal, unique, and strictly increasing. This
call evaluates exactly the printed values and never scans above 100:

``` r
schedule <- seq.int(25L, 100L, by = 5L)
print(schedule)
result <- INFLECT(dataset, set.i = schedule)
```

For an adaptive schedule with a hard upper bound, construct the values
explicitly:

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

Setting `zeroes.in = FALSE` excludes every non-positive value, warns
when negative values are present, and records per-marker exclusion
counts:

``` r
result$provenance$zero_handling
```

Canonical result names are `scores`, `qc_pass_rate`, `dip_pass`,
`iqr_pass`, `combined_pass`, and `criterion_pass`. `collection.U`,
`U.set`, `Unimodality`, `Accuracy.sets`, and `Accuracy.matrixes` are
deprecated aliases retained for compatibility.

See `vignette("using-fastINFLECT-1")` for the complete contract.

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
The fitted-model operational recommendation is nominal k = 44. Values
from k = 27 to 29 provide a lower-resolution structural sensitivity, and
k = 55 provides a spread sensitivity.

Nominal k = 44 contains 43 event-populated clusters because one SOM-only
cluster has zero events. Its independent final audit contains 56
detected, 118 ambiguous, 28 unresolved, and 986 no-detected
cluster-marker entries. These residual statuses preclude claims of
global unimodality or biological optimality. Exact comparisons and
provenance are in `inst/benchmarks/real-model-selection/`.

The allowed outcome labels are `detected_multimodality`, `ambiguous`,
`no_detected_multimodality`, and `unresolved_discrete`. The workflow
leaves tied measurements unchanged and reports “no detected
multimodality” without upgrading it to confirmed unimodality.

## Relationship to the original INFLECT package

fastINFLECT is derived from the INFLECT implementation developed by Jan
Verhoeff in the lab of JJ. Garcia-Vallejo and available from the
[GarciaVallejoLab
repository](https://github.com/jnverhoeff/GarciaVallejoLab). The
exported `INFLECT()` name is retained for API continuity. Memoised
hierarchy cuts, marker-at-a-time indexed scoring, deterministic
sampling, and compiled accelerators provide the faster implementation.
