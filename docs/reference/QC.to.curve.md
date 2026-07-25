# Plot a criterion-specific fastINFLECT QC curve

Fits the original four-parameter log-logistic diagnostic curve to
criterion-specific `qc_pass_rate` values and locates its geometric
inflection estimate. The estimate can fall outside the literal schedule;
it does not create a corresponding metaclustering partition.

## Usage

``` r
QC.to.curve(collection.U, basedata, ggtitle = NULL)
```

## Arguments

  - collection.U:
    
    Result from `iteration.QC`, a canonical score data frame with `k`
    and `qc_pass_rate`, or the deprecated score aliases.

  - basedata:
    
    `"Curve"` or `"Points"`.

  - ggtitle:
    
    Optional plot title.

## Value

A list containing canonical `scores`, `fittedcurve`, `lfunction`, and
`ggplot`. Deprecated `collection.U` is retained as a compatibility
alias.

## See also

`INFLECT`, `iteration.QC`, `Lfunction`
