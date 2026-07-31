# Plot marker-level QC pass rates

Plots criterion-specific marker QC across metaclusterings. The score is
the percentage of clusters where a marker passed the aggregate criterion
selected in
[`INFLECT()`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md).
It is a screening pass rate, not proof of unimodality. The colour
gradient denotes the tested cluster count.

## Usage

``` r
marker.performance(inflect.results, ggtitle = NULL, markers = NULL)
```

## Arguments

- inflect.results:

  An `inflect.results` object returned by
  [`INFLECT()`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md).

- ggtitle:

  Optional plot title.

- markers:

  Evaluated marker names to include. A `NULL` value displays all
  markers.

## Value

A list with `marker.dataframe` and `plot`. The data frame contains
canonical `k`, `marker`, and `qc_pass_rate` columns. Deprecated `i`,
`Marker`, and `Performance` aliases are retained for compatibility.

## Examples

``` r
flowsom <- system.file(
  "extdata",
  "Levine32sample.Rdata",
  package = "fastINFLECT"
)
load(flowsom)
result <- INFLECT(
  FlowSOM.results = dataset,
  set.i = 5:12
)

marker_qc <- marker.performance(result, ggtitle = "Levine32 sample")
marker_qc$plot
```
