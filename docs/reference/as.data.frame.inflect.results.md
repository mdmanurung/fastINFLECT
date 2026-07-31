# Extract tested fastINFLECT QC scores

Returns one row per tested `k`, with the aggregate `qc_pass_rate` and
the criterion selected by `uniform.test`.

## Usage

``` r
# S3 method for class 'inflect.results'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `inflect.results` object.

- row.names:

  `row.names` passed to `as.data.frame`.

- optional:

  `optional` passed to `as.data.frame`.

- ...:

  Additional arguments passed to `as.data.frame`.

## Value

Canonical `x$scores` as a `data.frame`.
