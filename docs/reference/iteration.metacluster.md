# Run iterative metaclustering on SOM-clustered objects

Computes metaclusterings for all cluster numbers in `set.i`. The Ward.D2
hierarchy is built once from the SOM codebook and then cut at each
requested k, rather than recomputing a separate clustering per k.

## Usage

``` r
iteration.metacluster(FlowSOM.results, set.i, multicore = FALSE, cores = NULL)
```

## Arguments

  - FlowSOM.results:
    
    A supported SOM object with completed SOM clustering. Supports
    FlowSOM objects and kohonen objects returned by `som` or `xyf`.

  - set.i:
    
    Literal, unique, strictly increasing integer cluster counts within
    the SOM-node range.

  - multicore:
    
    Retained for backward compatibility. Hierarchical clustering is now
    computed once and is not parallelized.

  - cores:
    
    Retained for backward compatibility.

## Value

metaclustering.list A `list` of `arrays` with metacluster-codes for the
SOM nodes within the input object.

## See also

`INFLECT`
