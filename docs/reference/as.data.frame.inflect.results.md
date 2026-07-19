# Coerce a fastINFLECT result to a data frame

Coerce a fastINFLECT result to a data frame

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

`x$collection.U` as a `data.frame`.
