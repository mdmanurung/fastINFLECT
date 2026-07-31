# Lfunction for determining Inflection Point

Wrapper around
[`leastError`](https://mdmanurung.github.io/fastINFLECT/reference/leastError.md)
with an option to refine the Inflection Point by limiting the range of
points evaluated for fitting line errors. Refinement is disabled by
default.

## Usage

``` r
Lfunction(totaldataframe, cutoff = 1000, plot = FALSE)
```

## Arguments

- totaldataframe:

  Data frame with cluster counts and corresponding QC pass rates, either
  observed points or a fitted curve.

- cutoff:

  Integer. Initial Inflection Point calculated on the entire curve is
  multiplied by `cutoff` to determine the new range of the curve that is
  used as input for
  [`leastError`](https://mdmanurung.github.io/fastINFLECT/reference/leastError.md)

- plot:

  Logical. If `TRUE`, draw the legacy diagnostic base plot.

## Value

`list` with 3 items: Inflection Point, the endpoint of the second
touchline, and the angle between the two touchlines

## See also

[`INFLECT`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md)
,
[`QC.to.curve`](https://mdmanurung.github.io/fastINFLECT/reference/QC.to.curve.md),[`leastError`](https://mdmanurung.github.io/fastINFLECT/reference/leastError.md)
