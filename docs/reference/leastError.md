# Determine kneepoint by minimizing fitting errors

For each candidate split of a point collection, fits two lines and
computes the weighted combined RMSE. Returns the split index with the
smallest combined error. This is the split criterion used by
`Lfunction`.

## Usage

``` r
leastError(dataframe)
```

## Arguments

  - dataframe:
    
    Data frame with cluster counts and corresponding QC pass rates,
    typically from `iteration.QC`.

## Value

Integer denoting the row of dataframe whose coordinate provides the
kneepoint

## See also

`INFLECT` , `QC.to.curve`,`Lfunction`
