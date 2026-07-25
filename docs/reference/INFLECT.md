# Run the fastINFLECT computation

Runs iterative metaclustering, separated marker-level dip/IQR quality
control, and diagnostic inflection-point estimation. The exported
function is named `INFLECT()` for continuity with the original INFLECT
API; the package uses a memoised engine that avoids the original
repeated per-k QC loop. `set.i` is a required, literal schedule:
fastINFLECT evaluates and materializes exactly the supplied cluster
counts. Complete finite transformed distributions, including negative
and zero values, are scored by default. A fitted inflection can estimate
an unscheduled k, but no partition is silently created for that
estimate.

## Usage

``` r
INFLECT(
  FlowSOM.results,
  set.i,
  multicore = FALSE,
  cores = NULL,
  zeroes.in = TRUE,
  only.clustering.markers = TRUE,
  acquired_markers = NULL,
  basedata = "Curve",
  ggtitle = NULL,
  uniform.test = c("both", "spread", "unimodality"),
  th.pvalue = 0.05,
  th.IQR = 2,
  verbose = FALSE,
  max.n.diptest = NULL,
  max.events.per.node = NULL,
  seed = 1L,
  target = 0.95,
  ...
)
```

## Arguments

  - FlowSOM.results:

    A supported SOM object with completed SOM clustering. Supports
    FlowSOM objects and kohonen objects returned by `som` or `xyf`.

  - set.i:
    
    Required vector of at least five unique, strictly increasing integer
    cluster counts within the SOM-node range. Values are used literally.
    See `inflect_adaptive_set_i` to construct an explicitly bounded
    adaptive schedule.

  - multicore:
    
    logical, should the QC sweep be run in parallel (fork-based
    `mclapply` over distinct SOM-node subtrees on Unix). Ignored on
    Windows. Default is `FALSE`.

  - cores:
    
    If `multicore == TRUE`, number of cores to be used. Must be at least
    2. If `NULL`, the number of detected cores minus one is used.

  - zeroes.in:
    
    If `TRUE` (default), retain negative, zero, and positive finite
    transformed values. If `FALSE`, exclude all non-positive values,
    record per-marker counts, and warn when negatives are present.

  - only.clustering.markers:
    
    If `TRUE` only evaluates markers specified as clustering markers.
    For kohonen objects this is the first data layer.

  - acquired\_markers:
    
    Vector of column names with marker data to be evaluated by
    fastINFLECT. Ignored if `only.clustering.markers == TRUE`

  - basedata:
    
    Data to be used to calculate inflection point, given as a string.
    Options are `Curve` and `Points`

  - ggtitle:
    
    Optional title for resulting diagnostic graph. Default `NULL`

  - uniform.test:
    
    What tests are performed per marker per cluster. Options are
    `"both"`, `"spread"`, or `"unimodality"`.

  - th.pvalue:
    
    Dip-test pass threshold. A pair passes the dip criterion when its
    p-value is at least this value. Default is `0.05`.

  - th.IQR:
    
    Threshold for rejecting marker distribution based on inter-quartile
    range. Default is arc-sinh transformed value of `2`.

  - verbose:
    
    `logical`, default is `FALSE`.

  - max.n.diptest:
    
    Optional integer cap on the number of events used per cluster-marker
    dip test. The dip test's power depends on sample size, so a cap can
    be used for a sensitivity analysis. Subsampling is seeded and
    deterministic; the default `NULL` retains all testable events.

  - max.events.per.node:

    Optional positive integer. If set, each SOM node is sampled once to
    at most this many events. That shared event sample is used for both
    dip and IQR scoring across all metaclusterings and workers. Default
    `NULL` uses all events.

  - seed:
    
    Non-negative integer base seed for optional event sampling. Default
    `1L`.

  - target:
    
    Target QC pass rate (fraction in `(0,1]` or percentage in `(1,100]`)
    used to report the smallest tested k reaching the threshold. Default
    `0.95`. See `inflect_threshold_k`.

  - ...:
    
    Arguments to pass to `dip.test` through `FlowSOMQC`.

## Value

An S3 `inflect.results` object containing criterion-specific `scores`,
separated dip/IQR/combined matrices and `qc.details`, the fitted curve,
literal metaclustering partitions, selection estimates, and consolidated
provenance. Deprecated score and accuracy aliases are retained for
migration. Running the individual function `iteration.metacluster` and
`QC.to.curve` might provide more options and flexibility.

## See also

`iteration.metacluster`, `iteration.QC`, `FlowSOMQC`, `QC.to.curve`,
`leastError`, `Lfunction`, `marker.performance`

## Examples

``` r
# Read in FlowSOM object from file. Downsampled clustering result of Levine32 dataset clustering.
# SOM-clustered to 375 clusters.
flowsom <- system.file("extdata", "Levine32sample.Rdata", package="fastINFLECT")
load(flowsom)
set.i<- 5:12
inflect.results<- INFLECT(FlowSOM.results= dataset, set.i= set.i, multicore=FALSE)

# Display diagnostic graph
inflect.results$ggplot

```
