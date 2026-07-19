## Verify that a hand-rolled p-value from dip() reproduces diptest::dip.test()'s
## table-based p-value BIT-EXACTLY, so we can swap the 12x-cheaper kernel in
## without changing any pass/fail decision.
suppressPackageStartupMessages(library(diptest))
data(qDiptab, package = "diptest")

## Replicated exactly from diptest::dip.test (table branch).
fast_dip_p <- function(x) {
  x <- sort(x[stats::complete.cases(x)])
  n <- length(x)
  D <- diptest::dip(x)
  if (n <= 3) return(1)
  dn <- dimnames(qDiptab)
  max.n <- max(nn <- as.integer(dn[["n"]]))
  P.s <- as.numeric(dn[["Pr"]])
  if (n >= max.n) {
    n1 <- n0 <- max.n; i2 <- i.n <- length(nn); f.n <- 0
  } else {
    n0 <- nn[i.n <- findInterval(n, nn)]
    n1 <- nn[(i2 <- i.n + 1)]
    f.n <- (n - n0) / (n1 - n0)
  }
  y.0 <- sqrt(n0) * qDiptab[i.n, ]
  y.1 <- sqrt(n1) * qDiptab[i2, ]
  sD  <- sqrt(n) * D
  1 - approx(y.0 + f.n * (y.1 - y.0), P.s, rule = 2, xout = sD)[["y"]]
}

set.seed(11)
maxdiff <- 0; maxdecdiff <- 0L; ncmp <- 0L
gens <- list(
  function(n) rnorm(n),
  function(n) c(rnorm(n %/% 2), rnorm(n - n %/% 2, 6)),   # bimodal
  function(n) rexp(n),                                    # skewed
  function(n) runif(n),
  function(n) round(rnorm(n), 1),                         # ties
  function(n) c(rep(0, n %/% 2), rlnorm(n - n %/% 2))     # zero-inflated
)
for (n in c(4,5,6,10,25,50,133,500,2000, 32288)) {
  for (g in gens) {
    for (rep in 1:12) {
      x <- g(n)
      p_ref  <- diptest::dip.test(x)$p.value
      p_fast <- fast_dip_p(x)
      maxdiff <- max(maxdiff, abs(p_ref - p_fast))
      for (th in c(0.01, 0.05, 0.1)) {
        maxdecdiff <- max(maxdecdiff, as.integer((p_ref >= th) != (p_fast >= th)))
      }
      ncmp <- ncmp + 1L
    }
  }
}
cat(sprintf("comparisons: %d\n", ncmp))
cat(sprintf("max |p_fast - p_ref|     : %.3e\n", maxdiff))
cat(sprintf("max pass/fail disagreement: %d  (0 == perfect)\n", maxdecdiff))

## speed
x <- rnorm(8000)
t1 <- system.time(for (i in 1:300) diptest::dip.test(x))[["elapsed"]]
t2 <- system.time(for (i in 1:300) fast_dip_p(x))[["elapsed"]]
cat(sprintf("dip.test: %.3f s   fast_dip_p: %.3f s   -> %.1fx\n", t1, t2, t1/t2))
