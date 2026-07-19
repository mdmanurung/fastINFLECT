// Compiled accelerators for the INFLECT QC engine.
//
// Scope is deliberately narrow. Profiling the (already memoised, pure-R) engine
// showed the dominant remaining cost is the inter-quartile-range spread test:
// stats::quantile() fully sorts every (cluster, marker) sample. The dip statistic
// itself is left to diptest's own tested C routine (diptest::dip); it is NOT
// re-implemented here.
//
//  * inflect_iqr_cpp()          - type-7 IQR in O(n) via std::nth_element,
//                                 bit-compatible with stats::quantile(type = 7).
//  * inflect_dip_pvalue_cpp()   - vectorised, allocation-light reproduction of the
//                                 table-interpolation p-value of diptest::dip.test().
//
// Both have pure-R equivalents in R/inflect-qc-core.R; the R paths are the
// reference implementation and the fallback when the package is used without its
// compiled code. Equivalence is asserted in the test suite.

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

// One type-7 sample quantile from a mutable buffer, using selection rather than a
// full sort. Matches stats::quantile(x, probs = p, type = 7):
//   h = (n - 1) * p; q = x[floor(h)] + (h - floor(h)) * (x[floor(h)+1] - x[floor(h)]).
static double quantile7(std::vector<double>& v, double p) {
  const int n = static_cast<int>(v.size());
  if (n == 0) return NA_REAL;
  if (n == 1) return v[0];
  const double h = (n - 1) * p;
  const int lo = static_cast<int>(std::floor(h));
  const double frac = h - lo;
  std::nth_element(v.begin(), v.begin() + lo, v.end());
  const double xlo = v[lo];
  if (frac == 0.0 || lo + 1 >= n) {
    return xlo;
  }
  // Everything to the right of position lo is >= xlo, so its minimum is the
  // (lo+1)-th order statistic. Use the same weighting expression as
  // stats::quantile(type = 7) so the result is bit-identical, not merely equal.
  const double xhi = *std::min_element(v.begin() + lo + 1, v.end());
  return (1.0 - frac) * xlo + frac * xhi;
}

// Inter-quartile range via selection (internal)
//
// @param x numeric vector
// @return the type-7 inter-quartile range (Q3 - Q1)
// @keywords internal
// [[Rcpp::export(.inflect_iqr_cpp)]]
double inflect_iqr_cpp(NumericVector x) {
  std::vector<double> v(x.begin(), x.end());
  const double q25 = quantile7(v, 0.25);
  const double q75 = quantile7(v, 0.75);
  return q75 - q25;
}

// Dip-test p-value by table interpolation (internal)
//
// Vectorised reproduction of the table branch of \code{diptest::dip.test()}.
//
// @param D numeric vector of dip statistics
// @param n integer vector of sample sizes
// @param qd the qDiptab matrix (rows = n, cols = Pr)
// @param nn integer vector of tabulated sample sizes (row labels of qd)
// @param Ps numeric vector of tabulated probabilities (column labels of qd)
// @return numeric vector of p-values
// @keywords internal
// [[Rcpp::export(.inflect_dip_pvalue_cpp)]]
NumericVector inflect_dip_pvalue_cpp(NumericVector D, IntegerVector n,
                                     NumericMatrix qd, IntegerVector nn,
                                     NumericVector Ps) {
  const int L = nn.size();
  const int M = Ps.size();
  const int N = D.size();
  const int max_n = nn[L - 1];
  NumericVector out(N);

  std::vector<double> grid(M);
  for (int idx = 0; idx < N; idx++) {
    const int ni = n[idx];
    const double Di = D[idx];
    // ISNAN matches R's is.na(): true for both NA and NaN, so this agrees with
    // the pure-R .inflect_dip_pvalue() reference on every input.
    if (ISNAN(Di) || ni <= 3) {
      out[idx] = 1.0;
      continue;
    }
    int i_n, i2, n0, n1;
    double f_n;
    if (ni >= max_n) {
      n0 = n1 = max_n;
      i_n = i2 = L - 1;
      f_n = 0.0;
    } else {
      // findInterval: largest i0 (0-based) with nn[i0] <= ni.
      int i0 = 0;
      while (i0 + 1 < L && nn[i0 + 1] <= ni) i0++;
      i_n = i0;
      i2 = i0 + 1;
      n0 = nn[i_n];
      n1 = nn[i2];
      f_n = static_cast<double>(ni - n0) / (n1 - n0);
    }
    const double s0 = std::sqrt(static_cast<double>(n0));
    const double s1 = std::sqrt(static_cast<double>(n1));
    for (int j = 0; j < M; j++) {
      const double y0 = s0 * qd(i_n, j);
      grid[j] = y0 + f_n * (s1 * qd(i2, j) - y0);
    }
    const double sD = std::sqrt(static_cast<double>(ni)) * Di;

    double pv;
    if (sD <= grid[0]) {
      pv = Ps[0];
    } else if (sD >= grid[M - 1]) {
      pv = Ps[M - 1];
    } else {
      int j = 0;
      while (j + 1 < M && grid[j + 1] <= sD) j++;
      const double gj = grid[j];
      const double t = (sD - gj) / (grid[j + 1] - gj);
      pv = Ps[j] + t * (Ps[j + 1] - Ps[j]);
    }
    out[idx] = 1.0 - pv;
  }
  return out;
}
