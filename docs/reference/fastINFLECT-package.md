# fastINFLECT: Fast discovery of the metaclustering endpoint for SOM objects

fastINFLECT reimplements the INFLECT method for selecting a FlowSOM or
kohonen metaclustering endpoint from criterion-specific marker QC. The
public `INFLECT()` workflow is retained; the implementation replaces
repeated per-k QC work with memoised SOM-node subtree scoring, faster
dip-test p-value lookup, and Rcpp accelerators.

The core INFLECT idea is to scan metaclusterings, calculate a
marker-level QC pass rate, and locate the inflection where additional
clusters stop improving that aggregate. fastINFLECT retains separate
dip-test and IQR-spread evidence and distinguishes fitted estimates from
materialised partitions. A pass rate is a screening result, not proof of
unimodality.

The original INFLECT implementation was developed by Jan Verhoeff in the
lab of JJ. Garcia-Vallejo and is available at
<https://github.com/jnverhoeff/GarciaVallejoLab>.

## See also

Useful links:

  - <https://mdmanurung.github.io/fastINFLECT>

  - <https://github.com/mdmanurung/fastINFLECT>

  - Report bugs at <https://github.com/mdmanurung/fastINFLECT/issues>
