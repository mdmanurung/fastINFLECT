# fastINFLECT: Fast discovery of the metaclustering endpoint for SOM objects

fastINFLECT is a fast reimplementation of the original INFLECT package
for selecting a FlowSOM or kohonen metaclustering endpoint from marker
unimodality. The public `INFLECT()` workflow and the original diagnostic
interpretation are retained; the implementation replaces repeated per-k
QC work with memoised SOM-node subtree scoring, faster dip-test p-value
lookup, and small Rcpp accelerators.

The original INFLECT method established the package's core idea: choose
k by scanning metaclusterings, measuring whether marker expression is
unimodal within each cluster, and locating the inflection point where
additional clusters stop improving that score. fastINFLECT acknowledges
that work and focuses this repository on making the same method
practical for dense sweeps and benchmarked comparisons.

## See also

Useful links:

  - <https://mdmanurung.github.io/INFLECT>

  - <https://github.com/mdmanurung/INFLECT>

  - Report bugs at <https://github.com/mdmanurung/INFLECT/issues>
