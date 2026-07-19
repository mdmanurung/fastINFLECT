# Metaclustering function based on hierarchical clustering

Orphaned function from FlowSOM, imported here.

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
