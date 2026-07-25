# Build an explicitly bounded adaptive k schedule

Constructs the legacy adaptive schedule explicitly: dense single-k steps
at low k, five-k steps at medium k, and ten-k steps at high k. `max_k`
is an explicit upper bound, and the returned schedule never exceeds it.

## Usage

``` r
inflect_adaptive_set_i(n_nodes, max_k, dense_until = 25L, medium_until = 100L)
```

## Arguments

  - n\_nodes:

    Number of SOM nodes.

  - max\_k:

    Largest metacluster count to include. Must not exceed `n_nodes`.

  - dense\_until:

    Last k in the dense, one-k part of the schedule.

  - medium\_until:

    Last k in the medium, five-k part of the schedule.

## Value

A strictly increasing integer vector suitable for `INFLECT(set.i=)`.
