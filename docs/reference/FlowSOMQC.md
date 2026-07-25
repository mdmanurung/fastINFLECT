# Score one SOM metaclustering with dip and IQR criteria

Computes separate Hartigan dip-test and IQR-spread evidence for every
cluster-marker pair. The returned matrix contains the pass decision
selected by `uniform.test`; its `qc.details` attribute contains
`dip_pass`, `iqr_pass`, `combined_pass`, p-values, IQRs, counts,
exclusions, and failure reasons. A pass is criterion-specific and does
not prove true unimodality.

## Usage

``` r
FlowSOMQC(
  FlowSOM.results,
  metaclustering,
  zeroes.in = TRUE,
  only.clustering.markers = TRUE,
  acquired_markers = NULL,
  uniform.test = c("both", "spread", "unimodality"),
  th.pvalue = 0.05,
  th.IQR = 2,
  max.n.diptest = NULL,
  seed = 1L,
  verbose = TRUE,
  ...
)
```

## Arguments

  - FlowSOM.results:

    A supported FlowSOM or kohonen SOM object.

  - metaclustering:

    Integer vector with one metacluster label per SOM node.

  - zeroes.in:

    If `TRUE` (default), retain all finite transformed values. If
    `FALSE`, exclude every non-positive value and warn when negative
    values are present.

  - only.clustering.markers:

    Evaluate only clustering markers.

  - acquired\_markers:

    Marker names used when `only.clustering.markers = FALSE`.

  - uniform.test:

    Aggregate criterion: `"both"` (dip and IQR), `"spread"` (IQR), or
    `"unimodality"` (dip).

  - th.pvalue:

    Dip-test pass threshold.

  - th.IQR:

    IQR pass threshold.

  - max.n.diptest:

    Optional dip-test sample cap of at least four.

  - seed:

    Non-negative seed. Simulated dip p-values and optional subsampling
    use deterministic subtree-marker streams and preserve the caller's
    RNG state.

  - verbose:

    Logical.

  - ...:

    Additional arguments passed to `dip.test`.

## Value

Invisibly, the selected-criterion logical matrix. Attributes
`qc.details` and `provenance` retain the separated evidence.

## See also

`INFLECT`, `iteration.QC`
