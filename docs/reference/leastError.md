# Determine kneepoint by minimizing fitting errors

For each candidate split of a point collection, fits two lines and
computes the weighted combined RMSE. Returns the split index with the
smallest combined error. This is the split criterion used by
[`Lfunction`](https://mdmanurung.github.io/fastINFLECT/reference/Lfunction.md).

## Usage

``` r
leastError(dataframe)
```

## Arguments

- dataframe:

  Data frame with cluster counts and corresponding QC pass rates,
  typically from
  [`iteration.QC`](https://mdmanurung.github.io/fastINFLECT/reference/iteration.QC.md).

## Value

Integer denoting the row of dataframe whose coordinate provides the
kneepoint

## See also

[`INFLECT`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md)
,
[`QC.to.curve`](https://mdmanurung.github.io/fastINFLECT/reference/QC.to.curve.md),[`Lfunction`](https://mdmanurung.github.io/fastINFLECT/reference/Lfunction.md)
