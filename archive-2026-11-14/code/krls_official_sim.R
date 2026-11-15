# Independent replication with the method authors' own package (KRLS, CRAN).
# Settings pinned to the pre-1.7 defaults in force when Eldor et al. (2023) was written:
#   sigma = number of covariates, exact solution, lambda by leave-one-out.
suppressMessages(library(KRLS))
args <- commandArgs(trailingOnly = TRUE); reps <- as.integer(args[1]); mode <- args[2]
MU <- 3.335; SD <- 0.751; N <- 474; RHO <- -0.17
truths <- list("Straight line" = list(f = function(x) x, r = .26),
               "Study 1's own quadratic (min 2.48)" = list(f = function(x) (x - 2.48)^2, r = .28),
               "Rises, then flat above the mean" = list(f = function(x) pmin(x, MU), r = .26),
               "A real downturn (peak at 4.0)" = list(f = function(x) -(x - 4.0)^2, r = .26))
set.seed(4711)
for (nm in names(truths)) {
  tr <- truths[[nm]]; out <- NULL
  for (i in seq_len(reps)) {
    z <- rnorm(N); x <- pmin(pmax(round((MU + SD * z) * 7) / 7, 1), 5)
    t <- RHO * z + sqrt(1 - RHO^2) * rnorm(N)
    m <- tr$f(x); m <- (m - mean(m)) / sd(m); y <- tr$r * m + sqrt(1 - tr$r^2) * rnorm(N)
    X <- cbind(x, t)
    f <- if (mode == "old") krls(X = X, y = y, sigma = ncol(X), approx = "none", print.level = 0)
         else suppressMessages(krls(X = X, y = y, print.level = 0))
    pw <- f$derivatives[, 1]
    hi <- x >= quantile(x, .975); lo <- x <= quantile(x, .025)
    l8 <- lowess(x, pw, f = .8, iter = 0); l5 <- lowess(x, pw, f = .5, iter = 0); l3 <- lowess(x, pw, f = .3, iter = 0)
    out <- rbind(out, c(top = mean(pw[hi]), bot = mean(pw[lo]), e8 = tail(l8$y, 1), e5 = tail(l5$y, 1), e3 = tail(l3$y, 1),
                        s8 = l8$y[1], pk = mean(x < l8$x[which.max(l8$y)]), lam = f$lambda, sig = f$sigma, ntop = sum(hi)))
  }
  cat(sprintf("%s | top<0: %.0f%% | line<0 at top (.8/.5/.3): %.0f%% / %.0f%% / %.0f%% | bottom<0: %.0f%% | peak pct: %.0f | median lambda %.1f, sigma %.2f | n top %d-%d\n",
      nm, 100 * mean(out[, "top"] < 0), 100 * mean(out[, "e8"] < 0), 100 * mean(out[, "e5"] < 0), 100 * mean(out[, "e3"] < 0),
      100 * mean(out[, "bot"] < 0), 100 * median(out[, "pk"]), median(out[, "lam"]), median(out[, "sig"]), min(out[, "ntop"]), max(out[, "ntop"])))
  flush.console()
}
