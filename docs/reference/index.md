# Package index

## Main workflow

<!-- end list -->

  - `INFLECT()` : Run the fastINFLECT computation
  - `FlowSOMQC()` : Score one SOM metaclustering with dip and IQR
    criteria
  - `marker.performance()` : Plot marker performance across
    metaclustering results

## Metaclustering and QC internals

<!-- end list -->

  - `iteration.metacluster()` : Run iterative metaclustering on
    SOM-clustered objects
  - `iteration.QC()` : Score marker-level QC criteria across
    metaclustering schedules
  - `metaClusteringhclust()` : Hierarchical (Ward) metaclustering of SOM
    codes
  - `QC.to.curve()` : Plot a criterion-specific fastINFLECT QC curve
  - `Lfunction()` : Lfunction for determining Inflection Point
  - `leastError()` : Determine kneepoint by minimizing fitting errors

## Selection helpers and diagnostics

<!-- end list -->

  - `computemode()` : Estimate the mode of a numeric vector
  - `inflect_adaptive_set_i()` : Build an explicitly bounded adaptive k
    schedule
  - `inflect_kneedle()` : Locate the knee of a diminishing-returns curve
    (Kneedle)
  - `inflect_threshold_k()` : Smallest k that reaches a target QC pass
    rate
  - `bimodality.coefficient()` : Sarle's bimodality coefficient

## Result methods

<!-- end list -->

  - `print(<inflect.results>)` : Print a fastINFLECT result
  - `summary(<inflect.results>)` : Summarize a fastINFLECT result
  - `plot(<inflect.results>)` : Plot a fastINFLECT result
  - `as.data.frame(<inflect.results>)` : Coerce a fastINFLECT result to
    a data frame
