# Sarle's bimodality coefficient

Computes Sarle's bimodality coefficient (BC) for a numeric sample, a
fast moment-based screen for bimodality: $$BC = (g^2 + 1) / (k + 3
(n-1)^2 / ((n-2)(n-3)))$$ where \\(g\\) is the sample skewness and
\\(k\\) the sample excess kurtosis. BC lies in `(0, 1]`; the benchmark
`5/9` (about 0.555) is the value for a uniform distribution, and larger
values indicate a more bimodal shape.

BC is offered as a lightweight, size-robust diagnostic that complements
the dip test used for the fastINFLECT unimodality score inherited from
the original INFLECT method. It is deliberately *not* the default QC
statistic: because it is driven by skewness and kurtosis it can flag
heavy-tailed or strongly skewed but genuinely unimodal markers (common
in cytometry) as bimodal. Use it to rank markers/clusters for follow-up,
and confirm with the dip test.

## Usage

``` r
bimodality.coefficient(x, na.rm = TRUE)
```

## Arguments

  - x:
    
    A numeric vector (missing values are dropped). At least four finite
    values are required.

  - na.rm:
    
    `logical`; drop missing values before computing. Default `TRUE`.

## Value

A single numeric bimodality coefficient, or `NA_real_` when there are
too few points or the sample has zero variance.

## See also

`FlowSOMQC`, `computemode`

## Examples

``` r
bimodality.coefficient(rnorm(500))                       # ~0.35, unimodal
#> [1] 0.3675891
bimodality.coefficient(c(rnorm(250), rnorm(250, 6)))     # high, bimodal
#> [1] 0.7405664
```
