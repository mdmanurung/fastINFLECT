# Determine kneepoint based on minimizing fitting errors

Evaluates the errors for two fitted lines on both parts of a collection
of points. Function calculates the errors for a split on each element
and selects the one with least amount of error. Internal function of
`Lfunction`

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
