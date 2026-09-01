# Run a fastINFLECT metaclustering scan

Evaluates every cluster count supplied in `set.i`, scores each
cluster-marker pair with separate dip and IQR criteria, and reports
candidate values of `k`. The result retains every tested partition, the
component QC evidence, and the settings needed to interpret the run. A
candidate is a screening result, not proof that a partition is
biologically valid.

## Usage

``` r
INFLECT(
  FlowSOM.results,
  set.i,
  workers = 1L,
  markers = NULL,
  uniform.test = c("both", "spread", "unimodality"),
  th.pvalue = 0.05,
  th.IQR = 2,
  max.n.diptest = NULL,
  max.events.per.node = NULL,
  seed = 1L,
  target = 0.95,
  progress = interactive()
)
```

## Arguments

- FlowSOM.results:

  A supported SOM object with completed SOM clustering. Supports FlowSOM
  objects and kohonen objects returned by
  [`som`](https://rdrr.io/pkg/kohonen/man/supersom.html) or
  [`xyf`](https://rdrr.io/pkg/kohonen/man/supersom.html).

- set.i:

  Required vector of at least five unique, strictly increasing integer
  cluster counts within the SOM-node range. Values are used literally.
  See
  [`inflect_adaptive_set_i`](https://mdmanurung.github.io/fastINFLECT/reference/inflect_adaptive_set_i.md)
  to construct an explicitly bounded adaptive schedule.

- workers:

  Number of QC workers. `1L` is serial; values above one use fork-based
  parallelism where supported.

- markers:

  Marker names to score. `NULL` uses the SOM clustering markers.

- uniform.test:

  Criterion used for the aggregate pass-rate curve: `"both"` requires
  dip and IQR to pass, `"unimodality"` uses dip, and `"spread"` uses
  IQR. Component evidence is retained in every run.

- th.pvalue:

  Dip-test pass threshold. A pair passes the dip criterion when its
  p-value is at least this value. Default is `0.05`.

- th.IQR:

  Threshold for rejecting marker distribution based on inter-quartile
  range. Default is arc-sinh transformed value of `2`.

- max.n.diptest:

  Optional integer cap on the number of events used per cluster-marker
  dip test. The dip test's power depends on sample size, so a cap can be
  used for a sensitivity analysis. Subsampling is seeded and
  deterministic; the default `NULL` retains all testable events.

- max.events.per.node:

  Optional positive integer. If set, each SOM node is sampled once to at
  most this many events. That shared event sample is used for both dip
  and IQR scoring across all metaclusterings and workers. Default `NULL`
  uses all events.

- seed:

  Non-negative integer base seed for optional event sampling. Default
  `1L`.

- target:

  Target QC pass rate (fraction in `(0,1]` or percentage in `(1,100]`)
  used to report the smallest tested k reaching the threshold. Default
  `0.95`. See
  [`inflect_threshold_k`](https://mdmanurung.github.io/fastINFLECT/reference/inflect_threshold_k.md).

- progress:

  Show stage messages and a serial progress bar. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

## Value

An S3 `inflect.results` object. Use
[`print()`](https://rdrr.io/r/base/print.html) for candidate values,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for the
pass-rate curve,
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) for
tested scores, `result$selection` for candidate metadata, and the
criterion matrices plus `qc.details` for marker-level interpretation.

## See also

[`iteration.metacluster`](https://mdmanurung.github.io/fastINFLECT/reference/iteration.metacluster.md),
[`iteration.QC`](https://mdmanurung.github.io/fastINFLECT/reference/iteration.QC.md),
[`FlowSOMQC`](https://mdmanurung.github.io/fastINFLECT/reference/FlowSOMQC.md),
[`QC.to.curve`](https://mdmanurung.github.io/fastINFLECT/reference/QC.to.curve.md),
[`leastError`](https://mdmanurung.github.io/fastINFLECT/reference/leastError.md),
[`Lfunction`](https://mdmanurung.github.io/fastINFLECT/reference/Lfunction.md),
[`marker.performance`](https://mdmanurung.github.io/fastINFLECT/reference/marker.performance.md)

## Examples

``` r
# Load the bundled, downsampled Levine32 FlowSOM object.
flowsom <- system.file(
  "extdata",
  "Levine32sample.Rdata",
  package = "fastINFLECT"
)
load(flowsom)

result <- INFLECT(
  FlowSOM.results = dataset,
  set.i = 5:12,
  uniform.test = "both"
)

result
#> fastINFLECT result
#>   knee: 9
#>   range: 12
#>   angle: 12.88413
#>   tested k values: 8
#>   tested k range: 5-12
#>   markers: 32
#>   provenance: uniform.test=both
#>   k estimates and tested thresholds:
#>     inflection k=9 (71.2% QC pass) [tested_partition]
#>     kneedle    k=7 (69.2% QC pass) [tested_partition]
#>     threshold  k=NA [not_available]
plot(result)

as.data.frame(result)
#>    k qc_pass_rate criterion
#> 1  5     56.25000  combined
#> 2  6     62.50000  combined
#> 3  7     69.19643  combined
#> 4  8     68.75000  combined
#> 5  9     71.18056  combined
#> 6 10     74.06250  combined
#> 7 11     77.55682  combined
#> 8 12     79.16667  combined
result$selection
#>       method  k qc_pass_rate_at_k directly_tested partition_available
#> 1 inflection  9          71.18056            TRUE                TRUE
#> 2    kneedle  7          69.19643            TRUE                TRUE
#> 3  threshold NA                NA           FALSE               FALSE
#>           k_status score_source target unimodality_at_k
#> 1 tested_partition       tested     NA         71.18056
#> 2 tested_partition       tested     NA         69.19643
#> 3    not_available         <NA>     95               NA
```
