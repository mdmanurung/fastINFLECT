# Package index

## Main workflow

<!-- end list -->

  - `INFLECT()` : Run the fastINFLECT computation
  - `FlowSOMQC()` : Cluster quality control using diptest and IQR check
  - `marker.performance()` : Plot marker performance across
    metaclustering results

## Metaclustering and QC internals

<!-- end list -->

  - `iteration.metacluster()` : Run iterative metaclustering on
    SOM-clustered objects
  - `iteration.QC()` : Run unimodality quality control on
    metaclusterings in parallel
  - `metaClusteringhclust()` : Metaclustering function based on
    hierarchical clustering
  - `QC.to.curve()` : Plot diagnostic fastINFLECT curve and find
    inflection point
  - `Lfunction()` : Lfunction for determing Inflection Point
  - `leastError()` : Determine kneepoint based on minimizing fitting
    errors

## Selection helpers and diagnostics

<!-- end list -->

  - `computemode()` : Estimate the mode of a numeric vector
  - `inflect_kneedle()` : Locate the knee of a diminishing-returns curve
    (Kneedle)
  - `inflect_threshold_k()` : Smallest k that reaches a target
    unimodality
  - `bimodality.coefficient()` : Sarle's bimodality coefficient

## Result methods

<!-- end list -->

  - `print(<inflect.results>)` : Print a fastINFLECT result
  - `summary(<inflect.results>)` : Summarize a fastINFLECT result
  - `plot(<inflect.results>)` : Plot a fastINFLECT result
  - `as.data.frame(<inflect.results>)` : Coerce a fastINFLECT result to
    a data frame
