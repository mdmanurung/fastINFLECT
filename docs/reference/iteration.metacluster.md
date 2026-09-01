# Run iterative metaclustering on SOM-clustered objects

Computes metaclusterings for all cluster numbers in `set.i`. The Ward.D2
hierarchy is built once from the SOM codebook and then cut at each
requested k, rather than recomputing a separate clustering per k.

## Usage

``` r
iteration.metacluster(FlowSOM.results, set.i)
```

## Arguments

- FlowSOM.results:

  A supported SOM object with completed SOM clustering. Supports FlowSOM
  objects and kohonen objects returned by
  [`som`](https://rdrr.io/pkg/kohonen/man/supersom.html) or
  [`xyf`](https://rdrr.io/pkg/kohonen/man/supersom.html).

- set.i:

  Literal, unique, strictly increasing integer cluster counts within the
  SOM-node range.

## Value

metaclustering.list A `list` of `arrays` with metacluster-codes for the
SOM nodes within the input object.

## See also

[`INFLECT`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md)
