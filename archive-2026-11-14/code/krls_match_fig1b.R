# Which settings reproduce the features of the paper's Figure 1?
# Study 1 stand-in built from Table 1: N 474; safety mean 3.335, SD .751, on a 0.05 grid (as digitised); r = .26;
# tenure mean 5.5, SD 7.3 (skewed), r with safety -.17. Authors' package; penalty by leave-one-out.
suppressMessages(library(KRLS))
a <- commandArgs(trailingOnly = TRUE); sigs <- as.numeric(strsplit(a[1], ",")[[1]]); reps <- as.integer(a[2])
N <- 474; z01 <- function(v) (v - mean(v)) / sd(v)
truths <- list(plateau = function(x) pmin((x - 2.48)^2, (4.3 - 2.48)^2), dip = function(x) ifelse(x < 4.3, (x - 2.48)^2, (4.3 - 2.48)^2 - 2 * (x - 4.3)^2))
cross <- function(x, y, up) { s <- which(diff(sign(y)) != 0); if (!length(s)) return(NA)
  xs <- x[s] - y[s] * (x[s + 1] - x[s]) / (y[s + 1] - y[s]); if (up) xs[1] else xs[length(xs)] }
for (sg in sigs) for (nm in names(truths)) {
  set.seed(2026); out <- NULL
  for (i in seq_len(reps)) {
    z <- rnorm(N); x <- pmin(pmax(round((3.335 + .751 * z) / .05) * .05, 1), 5)
    u <- -.17 * z + sqrt(1 - .17^2) * rnorm(N); ten <- qgamma(pnorm(u), shape = (5.5 / 7.3)^2, scale = 7.3^2 / 5.5)
    y <- pmin(pmax(round((3.695 + .752 * z01(.27 * z01(truths[[nm]](x)) + sqrt(1 - .27^2) * rnorm(N))) * 4) / 4, 1), 5)
    f <- krls(X = cbind(x, ten), y = y, sigma = sg, approx = "none", print.level = 0)
    pw <- f$derivatives[, 1]; l <- lowess(x, pw, f = .8, iter = 0); k <- which.max(l$y)
    lo <- if (k > 1) cross(l$x[1:k], l$y[1:k], FALSE) else NA; hi <- if (k < length(l$x)) cross(l$x[k:length(l$x)], l$y[k:length(l$y)], TRUE) else NA
    out <- rbind(out, c(start = l$y[1], locross = lo, pk = l$y[k], pkx = l$x[k], hicross = hi, end = tail(l$y, 1), dmax = max(pw), dmin = min(pw), r = cor(x, y)))
  }
  q <- function(v) sprintf("%.2f [%.2f, %.2f]", median(v, na.rm = TRUE), quantile(v, .1, na.rm = TRUE), quantile(v, .9, na.rm = TRUE))
  cat(sprintf("sigma2=%g | %s | r %.2f | start %s | low crossing %s (%.0f%% have one) | peak %s at %s | high crossing %s (%.0f%%) | end %s | dots max %s min %s\n",
      sg, nm, mean(out[, "r"]), q(out[, "start"]), q(out[, "locross"]), 100 * mean(!is.na(out[, "locross"])), q(out[, "pk"]), q(out[, "pkx"]),
      q(out[, "hicross"]), 100 * mean(!is.na(out[, "hicross"])), q(out[, "end"]), q(out[, "dmax"]), q(out[, "dmin"])))
}
