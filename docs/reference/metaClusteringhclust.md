# Hierarchical (Ward) metaclustering of SOM codes

Performs hierarchical metaclustering of SOM cluster codes using Ward's
D2 linkage, adapted from FlowSOM.

## Usage

``` r
metaClusteringhclust(data, nClus)
```

## Arguments

  - data:
    
    Matrix with the median values for each clustering parameter for all
    SOM clusters.

  - nClus:
    
    Amount of metaclusters to be obtained.

## Value

A `scalar` of `length(SOMclusters)` with metacluster identity.

## See also

`INFLECT`
