# fastINFLECT: Fast discovery of the metaclustering endpoint for SOM objects

fastINFLECT evaluates a literal schedule of metacluster counts for a
FlowSOM or kohonen self-organising map. Start with
[`INFLECT()`](https://mdmanurung.github.io/fastINFLECT/reference/INFLECT.md),
then use [`print()`](https://rdrr.io/r/base/print.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to read
the candidate values and tested QC pass rates.

The result retains the tested partitions, separate dip and IQR
decisions, detailed cluster-marker measurements, and run provenance.
Candidate values narrow the search for `k`; they do not establish
biological validity or prove unimodality.

## See also

Useful links:

- <https://mdmanurung.github.io/fastINFLECT>

- <https://github.com/mdmanurung/fastINFLECT>

- Report bugs at <https://github.com/mdmanurung/fastINFLECT/issues>
