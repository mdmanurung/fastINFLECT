# Run iterative metaclustering on SOM-clustered objects

Metaclustering runs are determined for all cluster numbers in `set.i`.
In the original INFLECT workflow, each requested k was handled as a
separate scan point. fastINFLECT computes the Ward.D2 hierarchy once and
cuts it for all requested cluster numbers, preserving the hierarchical
strategy while avoiding repeated work.

## Usage

``` r
iteration.metacluster(FlowSOM.results, set.i, multicore = TRUE, cores = NULL)
```

## Arguments

  - FlowSOM.results:
    
    A supported SOM object with completed SOM clustering. Supports
    FlowSOM objects and kohonen objects returned by `som` or `xyf`.

  - set.i:
    
    Vector containing either the desired iterations to be tested

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
