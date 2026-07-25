# BMV cluster-number synthesis

- Generated: 2026-07-25T16:58:42+0200
- Full model: 39,050,953 events, 900 SOM nodes, 27 markers.
- Model MD5: `784b3a6048c379ca64bbe8d4c0fb063d`.

## Recommendation

- Use nominal **k=44** as the operational fitted-model partition. It is the combined and dip Kneedle and was directly materialised and independently audited.
- Treat k=27 as the lower-resolution parametric elbow and k=55 as the spread/variance sensitivity endpoint; the defensible sensitivity range is therefore 44-55.
- k=44 contains 43 event-populated clusters and 1 empty nominal cluster; its largest cluster contains 23.22% of events.
- Final modality audit at k=44: 56 detected, 118 ambiguous, 28 unresolved, and 986 no-detected cluster-marker entries.
- No directly audited k supports a global unimodality claim.

## Method comparison

- Combined INFLECT: inflection k=27; Kneedle k=44.
- Dip-only: inflection k=28; Kneedle k=44.
- IQR spread: inflection k=38; Kneedle k=55.
- Consensus metaclustering: The first five-k consensus plateau starts at k=29.
- At k=44, consensus minimum seed ARI=0.999 and median relative delta area=0.0012.
- Fitted-model internal metrics: silhouette k=29, Calinski-Harabasz k=25, Davies-Bouldin k=100, unweighted dispersion k=54, event-weighted dispersion k=55.
- Historical intended 80/20 code space: silhouette k=26, Calinski-Harabasz k=25, Davies-Bouldin k=92, unweighted dispersion k=54, event-weighted dispersion k=61.
- The captured `NaNs produced` warning was reproducibly classified as temporary LL.4 optimizer exploration: all final fitted values were finite and the rerun inflections exactly matched 28/38/27.
- Silhouette, Calinski-Harabasz, Davies-Bouldin, and explained-dispersion recommendations are retained in `metric-recommendations.csv`; they diagnose different trade-offs and are not interchangeable votes.

## Interpretation limits

- IQR spread is a dispersion screen, not a modality test.
- Dip plus ACR evidence addresses univariate marker marginals with BH correction within each 27-marker cluster/method/seed/sample-size family.
- Consensus stability measures reproducibility conditional on k; it does not establish biological correctness or select k by itself.
- The fitted model's effective layer weighting differs from the historical intended 80/20 code space, so the historical k=50/62 partitions remain sensitivity comparators rather than substitutes for the fitted-model recommendation.
- The slow 5,000-event primary refinement tail was resumed from immutable per-task checkpoints with a wider logical sharding; hashes and job mappings are retained in `refinement-reshard-provenance.tsv`.

See `audited-inflect-k-comparison.csv` and `comparator-solution-quality.csv` for exact counts and metrics.
