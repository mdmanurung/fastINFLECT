# Cluster quality control using diptest and IQR check

Unimodality score is calculated for cluster results. Per marker per
cluster [dip.test](https://rdrr.io/pkg/diptest/man/dip.test.html) is
applied and inter-quartile range is assessed. This function preserves
the original INFLECT marker-level QC criterion so users can inspect or
reproduce the base statistic directly; the faster package-level sweep
reuses the same criterion through memoised helper code.

## Usage

``` r
FlowSOMQC(
  FlowSOM.results,
  metaclustering,
  zeroes.in = FALSE,
  only.clustering.markers = TRUE,
  acquired_markers = NULL,
  uniform.test = c("both", "spread", "unimodality"),
  th.pvalue = 0.05,
  th.IQR = 2,
  verbose = TRUE,
  ...
)
```

## Arguments

  - FlowSOM.results:
    
    A supported SOM object with completed SOM clustering. Supports
    FlowSOM objects and kohonen objects returned by `som` or `xyf`.

  - metaclustering:
    
    Vector with metacluster codes for all SOM-clusters.

  - zeroes.in:
    
    Should be values at and below `0` be included. Recommended default
    for mass cytometry data is `FALSE`

  - only.clustering.markers:
    
    If `TRUE` only evaluates markers specified as clustering markers.
    For kohonen objects this is the first data layer.

  - acquired\_markers:
    
    Vector of column names with marker data to be evaluated by
    fastINFLECT. Ignored if `only.clustering.markers == TRUE`

  - uniform.test:
    
    What tests are performed per marker per cluster. Options are "both",
    "spread" , or "unimodality" as a string.

  - th.pvalue:
    
    Threshold for rejecting Unimodality dip.test result. Default is
    `0.05`. For more information see
    [dip.test](https://rdrr.io/pkg/diptest/man/dip.test.html)

  - th.IQR:
    
    Threshold for rejecting marker distribution based on inter-quartile
    range. Default is arc-sinh transformed value of `2`.

  - verbose:
    
    `logical` , default is `TRUE`

  - ...:
    
    Additional arguments to pass to `dip.test`.

## Value

A `matrix` with evaluated markers in columns and clusters in rows. Each
position in the matrix is `logical` indicating a pass or a fail.

## See also

`INFLECT` , `iteration.QC`
