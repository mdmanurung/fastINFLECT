# Score marker-level QC criteria across metaclustering schedules

Scores every requested metaclustering using separate Hartigan dip-test
and IQR-spread criteria. Each distinct SOM-node subtree is evaluated
once and the per-k matrices are assembled from that cache. The aggregate
`qc_pass_rate` is explicitly tied to the criterion selected by
`uniform.test`; it is not evidence that a distribution is truly
unimodal.

## Usage

``` r
iteration.QC(
  FlowSOM.results,
  metaclustering.list,
  set.i,
  multicore = FALSE,
  cores = NULL,
  zeroes.in = TRUE,
  only.clustering.markers = TRUE,
  acquired_markers = NULL,
  uniform.test = c("both", "spread", "unimodality"),
  th.pvalue = 0.05,
  th.IQR = 2,
  verbose = FALSE,
  max.n.diptest = NULL,
  max.events.per.node = NULL,
  seed = 1L,
  ...
)
```

## Arguments

- FlowSOM.results:

  A supported SOM object with completed SOM clustering. Supports FlowSOM
  objects and kohonen objects returned by
  [`som`](https://rdrr.io/pkg/kohonen/man/supersom.html) or
  [`xyf`](https://rdrr.io/pkg/kohonen/man/supersom.html).

- metaclustering.list:

  Named list containing exactly one node-label vector for every value in
  `set.i`.

- set.i:

  Literal, unique, strictly increasing integer cluster counts.

- multicore:

  Logical. On Unix, use fork-based
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html) over distinct
  SOM-node subtrees. Other platforms validate the same arguments and use
  a recorded serial fallback. Default `FALSE`.

- cores:

  Worker count when `multicore = TRUE`; must be at least two.

- zeroes.in:

  Logical. If `TRUE` (the default), retain negative, zero, and positive
  finite transformed values. If `FALSE`, every non-positive value is
  excluded, per-marker counts are recorded, and negative inputs trigger
  a warning.

- only.clustering.markers:

  If `TRUE`, evaluate only clustering markers.

- acquired_markers:

  Marker names to evaluate when `only.clustering.markers = FALSE`.

- uniform.test:

  Aggregate criterion: `"both"` selects the combined dip and IQR pass,
  `"spread"` selects IQR only, and `"unimodality"` selects the dip test
  only. Both component tests are always retained in `qc.details`.

- th.pvalue:

  Dip-test pass threshold. A cell passes when `p_value >= th.pvalue`.

- th.IQR:

  IQR-spread pass threshold. A cell passes when `IQR < th.IQR`.

- verbose:

  Logical.

- max.n.diptest:

  Optional positive dip-test sample cap of at least four. Sampling is
  deterministic per subtree and marker.

- max.events.per.node:

  Optional positive integer. Each SOM node is sampled once before
  subtree assembly; the retained rows feed every marker and both tests.

- seed:

  Non-negative base seed.

- ...:

  Additional arguments passed to
  [`dip.test`](https://rdrr.io/pkg/diptest/man/dip.test.html).

## Value

A list containing canonical `scores`, separate named lists of
`dip_pass`, `iqr_pass`, `combined_pass`, and selected `criterion_pass`
matrices, full `qc.details`, and `provenance`. Deprecated `U.set` and
`Accuracy.matrixes` aliases are retained for compatibility.

## See also

[`INFLECT`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md),
[`iteration.metacluster`](https://mdmanurung.github.io/fastINFLECT/reference/iteration.metacluster.md),
[`FlowSOMQC`](https://mdmanurung.github.io/fastINFLECT/reference/FlowSOMQC.md)
