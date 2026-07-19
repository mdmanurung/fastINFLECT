# Lfunction for determing Inflection Point

Wrapper function for `leastError` with an option for refining the
Inflection Point, by limiting the range of of points to evaluate for
fitting line errors. Default is to not use the refinement.

## Usage

``` r
Lfunction(totaldataframe, cutoff = 1000, plot = FALSE)
```

## Arguments

  - totaldataframe:
    
    Dataframe with set.i and corresponding unimodality scores, resulting
    from `iteration.QC`. This can be the calculcated points collection.U
    or the fitted curve.

  - cutoff:
    
    Integer. Initial Inflection Point calculated on the entire curve is
    multiplied by `cutoff` to determine the new range of the curve that
    is used as input for `leastError`

  - plot:
    
    Logical. If `TRUE`, draw the legacy diagnostic base plot.

## Value

`list` with 3 items: Inflection Point, the endpoint of the second
touchline, and the angle between the two touchlines

## See also

`INFLECT` , `QC.to.curve`,`leastError`
