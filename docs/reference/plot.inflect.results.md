# Plot the fastINFLECT QC pass-rate curve

Returns the diagnostic plot stored during
[`INFLECT()`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md).
Use the curve to see how the selected criterion changes across the
literal `set.i` schedule.

## Usage

``` r
# S3 method for class 'inflect.results'
plot(x, ...)
```

## Arguments

- x:

  An `inflect.results` object.

- ...:

  Unused.

## Value

The stored diagnostic `ggplot` object.
