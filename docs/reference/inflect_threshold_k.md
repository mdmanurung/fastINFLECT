# Smallest k that reaches a target QC pass rate

Returns the smallest tested number of metaclusters whose
criterion-specific QC pass rate meets or exceeds `target`. This
aggregate is a screening metric; it does not establish that clusters are
truly unimodal.

## Usage

``` r
inflect_threshold_k(collection.U, target = 0.95)
```

## Arguments

  - collection.U:
    
    Canonical score data frame with `k` and `qc_pass_rate`, or the
    deprecated `i` and `Unimodality` aliases.

  - target:
    
    Desired pass rate. Values in `(0, 1]` are read as a fraction and
    values in `(1, 100]` as a percentage. Default `0.95`.

## Value

The smallest `i` meeting the target, or `NA_integer_` if no tested k
reaches it.

## See also

`inflect_kneedle`, `INFLECT`
