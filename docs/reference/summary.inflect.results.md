# Summarize a fastINFLECT result

Returns one row with the tested range, QC pass-rate range, run settings,
and the inflection, Kneedle, and threshold candidates.

## Usage

``` r
# S3 method for class 'inflect.results'
summary(object, ...)
```

## Arguments

- object:

  An `inflect.results` object.

- ...:

  Unused.

## Value

A one-row `data.frame` of scan settings and candidate values.
