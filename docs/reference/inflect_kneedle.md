# Locate the knee of a diminishing-returns curve (Kneedle)

Nonparametric knee/elbow detector following Satopää et al. (2011). For a
concave, increasing curve (as the fastINFLECT unimodality sweep
typically is) the knee is the x at which the normalised curve is
farthest above the diagonal joining its endpoints. Unlike the
four-parameter log-logistic fit used by `QC.to.curve`, it assumes no
functional form, which makes it a robust cross-check on dense sweeps.

## Usage

``` r
inflect_kneedle(x, y, concave = TRUE, increasing = TRUE)
```

## Arguments

  - x:
    
    Numeric vector of x coordinates (e.g. number of metaclusters).

  - y:
    
    Numeric vector of y coordinates (e.g. unimodality score).

  - concave:
    
    `logical`; `TRUE` (default) for a concave curve (diminishing
    returns). Set `FALSE` for a convex curve.

  - increasing:
    
    `logical`; `TRUE` (default) if `y` rises with `x`.

## Value

The `x` value at the detected knee, or `NA` if it cannot be determined
(fewer than three finite points or a degenerate range).

## References

Satopää, V., Albrecht, J., Irwin, D., & Raghavan, B. (2011). Finding a
"kneedle" in a haystack. *31st ICDCS Workshops*, 166-171.

## See also

`inflect_threshold_k`, `QC.to.curve`
