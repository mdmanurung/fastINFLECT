# Score one SOM metaclustering with dip and IQR criteria

Computes separate Hartigan dip-test and IQR-spread evidence for every
cluster-marker pair. The returned matrix contains the pass decision
selected by `uniform.test`; its `qc.details` attribute contains
`dip_pass`, `iqr_pass`, `combined_pass`, p-values, IQRs, counts, and
failure reasons. All selected marker values must be finite; negative and
zero values are retained unchanged. A pass is criterion-specific and
does not prove true unimodality.

## Usage

``` r
FlowSOMQC(
  FlowSOM.results,
  metaclustering,
  markers = NULL,
  uniform.test = c("both", "spread", "unimodality"),
  th.pvalue = 0.05,
  th.IQR = 2,
  max.n.diptest = NULL,
  seed = 1L,
  progress = interactive()
)
```

## Arguments

- FlowSOM.results:

  A supported FlowSOM or kohonen SOM object.

- metaclustering:

  Integer vector with one metacluster label per SOM node.

- markers:

  Marker names to score. `NULL` uses the SOM clustering markers.

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

  Non-negative seed. Simulated dip p-values and optional subsampling use
  deterministic subtree-marker streams and preserve the caller's RNG
  state.

- progress:

  Show scoring progress. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

## Value

Invisibly, the selected-criterion logical matrix. Attributes
`qc.details` and `provenance` retain the separated evidence.

## See also

[`INFLECT`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md),
[`iteration.QC`](https://mdmanurung.github.io/fastINFLECT/reference/iteration.QC.md)
