# Settings grid: does the result depend on the kernel width or the penalty?
# Study 1 stand-in (474 people, 1-5 scale, r about .26). Same datasets in every cell (same seed).
suppressMessages(library(KRLS))
a <- commandArgs(trailingOnly = TRUE); mode <- a[1]; reps <- as.integer(a[2])
MU <- 3.335; SD <- .751; N <- 474; RHO <- -.17; R <- .26
z01 <- function(v) (v - mean(v)) / sd(v)
truths <- list("straight line" = function(x) x, "paper's reading (peak at 4.75)" = function(x) -(x - 4.75)^2)
cells <- if (mode == "sigma") lapply(c(1, 2, 4, 8), function(s) list(sigma = s, lambda = NULL)) else
         lapply(c(0.1, 1, 5, 15, 50), function(l) list(sigma = 2, lambda = l))
for (cl in cells) for (nm in names(truths)) {
  set.seed(777); out <- NULL
  for (i in seq_len(reps)) {
    z <- rnorm(N); x <- pmin(pmax(round((MU + SD * z) * 7) / 7, 1), 5)
    t <- RHO * z + sqrt(1 - RHO^2) * rnorm(N)
    y <- pmin(pmax(round((3.695 + .752 * z01(R * z01(truths[[nm]](x)) + sqrt(1 - R^2) * rnorm(N))) * 4) / 4, 1), 5)
    f <- krls(X = cbind(x, t), y = y, sigma = cl$sigma, lambda = cl$lambda, approx = "none", print.level = 0)
    pw <- f$derivatives[, 1]; hi <- x >= quantile(x, .975); lo <- x <= quantile(x, .025)
    l8 <- lowess(x, pw, f = .8, iter = 0)
    out <- rbind(out, c(top = mean(pw[hi]), bot = mean(pw[lo]), e8 = tail(l8$y, 1), pk = mean(x < l8$x[which.max(l8$y)]), lam = f$lambda,
                        mid = mean(pw[x > quantile(x, .25) & x < quantile(x, .75)])))
  }
  cat(sprintf("sigma2=%g | lambda=%s | %s | top<0: %.0f%% | line(.8)<0 at top: %.0f%% | bottom<0: %.0f%% | peak pct %.0f | mid %+.3f top %+.3f | median lambda %.1f | %d runs\n",
      cl$sigma, if (is.null(cl$lambda)) "LOO" else format(cl$lambda), nm, 100 * mean(out[, "top"] < 0), 100 * mean(out[, "e8"] < 0), 100 * mean(out[, "bot"] < 0),
      100 * median(out[, "pk"]), mean(out[, "mid"]), mean(out[, "top"]), median(out[, "lam"]), reps))
}
