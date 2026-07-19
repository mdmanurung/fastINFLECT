# Smallest k that reaches a target unimodality

Returns the smallest number of metaclusters whose unimodality score
meets or exceeds `target`. This fastINFLECT addition operationalises the
goal of admitting no cluster with residual bimodal marker expression
while avoiding over-clustering: it is the first k at which (nearly)
every (cluster, marker) pair is unimodal.

## Usage

``` r
inflect_threshold_k(collection.U, target = 0.95)
```

## Arguments

  - collection.U:
    
    Data frame with columns `i` and `Unimodality` (the percentage of
    unimodal (cluster, marker) pairs), as returned in
    `inflect.results$collection.U`.

  - target:
    
    Desired unimodality. Values in `(0, 1]` are read as a fraction and
    values in `(1, 100]` as a percentage. Default `0.95`.

## Value

The smallest `i` meeting the target, or `NA_integer_` if no tested k
reaches it.

## See also

`inflect_kneedle`, `INFLECT`
