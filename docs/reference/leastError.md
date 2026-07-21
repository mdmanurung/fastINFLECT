# Determine kneepoint based on minimizing fitting errors

Evaluates the errors for two fitted lines on both parts of a collection
of points. The helper calculates the error for each possible split and
selects the split with the smallest combined error. This is the split
criterion used by `Lfunction`.

## Usage

``` r
leastError(dataframe)
```

## Arguments

  - dataframe:
    
    Dataframe with set.i and corresponding unimodality scores, resulting
    from `iteration.QC`

## Value

Integer denoting the row of dataframe whose coordinate provides the
kneepoint

## See also

`INFLECT` , `QC.to.curve`,`Lfunction`
