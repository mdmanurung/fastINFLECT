# fastINFLECT

Fast quality control and determination of optimal k for FlowSOM and kohonen
high-dimensional SOM clustering data.

fastINFLECT reimplements the original INFLECT package for speed. It scans
metaclusterings for marker unimodality, locates the inflection point on the
fitted unimodality curve, and replaces the expensive repeated per-k QC loop
with memoised hierarchy cuts, faster dip-test p-value lookup, and small
compiled accelerators. On the fitted curve of this unimodality score,
fastINFLECT determines the point where unimodality stops increasing for greater
numbers of clusters.

### Installation

fastINFLECT is not on CRAN; install it from GitHub. It depends on FlowSOM,
which is distributed through Bioconductor.

The simplest path uses `pak`, which resolves the Bioconductor dependency for
you:

``` r
# install.packages("pak")
pak::pak("mdmanurung/fastINFLECT")
```

To use `remotes` instead, install FlowSOM first:

``` r
# install.packages(c("BiocManager", "remotes"))
BiocManager::install("FlowSOM")
remotes::install_github("mdmanurung/fastINFLECT")
```

The `kohonen` package is optional and only needed to run fastINFLECT on
kohonen SOM objects:

``` r
install.packages("kohonen")
```

### Relationship to original INFLECT

fastINFLECT is derived from the original INFLECT implementation developed by
Jan Verhoeff in the lab of JJ. Garcia-Vallejo at the
[GarciaVallejoLab](https://github.com/jnverhoeff/GarciaVallejoLab). The package
keeps the original marker-unimodality idea and its public `INFLECT()` workflow
as the main user entry point. The new package name reflects the focus of this
reimplementation: making the original algorithm fast enough for dense k sweeps,
explicit benchmark comparisons, and direct inspection of residual bimodal marker
expression.

### What's new in 1.0.0

- **Much faster, legacy-matched results.** The QC sweep memoises each distinct SOM-node
  set (hierarchical metaclusterings are nested, so there are only ~`2·nNodes`
  distinct clusters across *any* sweep) and uses `diptest::dip()` with
  legacy-matched p-value interpolation plus small Rcpp accelerators. A full dense
  sweep of every k on the bundled Levine32 SOM (375 nodes) runs in a few seconds
  instead of hours while preserving the original QC decisions. Dense sweeps being
  cheap means the diagnostic curve and its inflection point are far better resolved.
- **Clearer answer to "how many clusters".** `INFLECT()` reports three recommended
  k values: the LL.4 `inflection` (default), a nonparametric `kneedle` knee, and a
  `threshold`, the smallest k reaching a target unimodality (`target=`), which
  directly encodes "no cluster with residual bimodal marker expression, without
  over-clustering".
- **Size-robust scoring.** `max.n.diptest=` caps the events per dip test (seeded,
  deterministic) so large clusters no longer make the dip test over-powered.
- **`bimodality.coefficient()`** is exported as a fast complementary diagnostic.
- **Benchmark evidence.** The benchmark vignette compares fastINFLECT 1.0 against the
  original INFLECT loop and FlowSOM consensus, and includes marker-expression
  histograms for residual bimodality versus overlapping fastINFLECT unimodal clusters.
- **Levine metacluster quality audit.** The comparison vignette includes
  all-cluster histograms of per-metacluster marker-unimodality pass rates for
  FlowSOM consensus, FlowSOM auto-k, hierarchical cuts, and fastINFLECT.
