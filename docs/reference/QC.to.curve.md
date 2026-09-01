# Plot a criterion-specific fastINFLECT QC curve

Fits the original four-parameter log-logistic diagnostic curve to
criterion-specific `qc_pass_rate` values and locates its geometric
inflection estimate. The estimate can fall outside the literal schedule;
it does not create a corresponding metaclustering partition.

## Usage

``` r
QC.to.curve(collection.U)
```

## Arguments

- collection.U:

  Result from
  [`iteration.QC`](https://mdmanurung.github.io/fastINFLECT/reference/iteration.QC.md),
  a canonical score data frame with `k` and `qc_pass_rate`, or the
  deprecated score aliases.

## Value

A list containing canonical `scores`, `fittedcurve`, `lfunction`, and
`ggplot`. Deprecated `collection.U` is retained as a compatibility
alias.

## See also

[`INFLECT`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md),
[`iteration.QC`](https://mdmanurung.github.io/fastINFLECT/reference/iteration.QC.md),
[`Lfunction`](https://mdmanurung.github.io/fastINFLECT/reference/Lfunction.md)
