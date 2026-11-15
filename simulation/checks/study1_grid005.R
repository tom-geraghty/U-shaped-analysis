# Study 1 re-run in the official KRLS package, with a third truth (rises then flat) and every run saved to CSV.
# Derived from krls_official_ext.R (14 Nov); only the truths list and the CSV output were added.
# Extension: each study on its own scale, lumpy and skewed ratings, fixed effects in the kernel.
# Authors' package, defaults as when the paper was written (sigma = number of covariates, exact, LOO).
suppressMessages(library(KRLS))
a <- commandArgs(trailingOnly = TRUE); study <- a[1]; reps <- as.integer(a[2]); set.seed(as.integer(a[3]))
z01 <- function(v) (v - mean(v)) / sd(v)
S <- list(
 S1 = list(N = 474, mu = 3.335, sd = .751, lo = 1, hi = 5, step = .05, r = .26, rho = -.17, orgs = 0, peak = 4.75,
           yfun = function(l) pmin(pmax(round((3.695 + .752 * l) * 4) / 4, 1), 5)),               # 4-item mean, 1 to 5
 S2 = list(N = 301, mu = 4.020, sd = 1.306, lo = 1, hi = 7, step = 1/7, r = .19, rho = .01, orgs = 0, peak = 6.3,
           yfun = function(l) as.numeric(cut(rank(l) / length(l), c(0, .2, .4, .6, .8, 1)))),      # five flat categories, as printed
 S5 = list(N = 800, mu = 3.390, sd = .425, lo = 1, hi = 5, step = 0, r = .18, rho = -.08, orgs = 0, peak = 4.09,
           yfun = function(l) { s <- exp(-.6 * l); pmin(10, 10 - .88 * s / mean(s)) }),            # skewed towards a ceiling of 10
 S3 = list(N = 365, mu = 3.35, sd = .30, lo = 1, hi = 5, step = 1/7, r = .15, rho = 0, orgs = 4, peak = 3.7,
           yfun = function(l) pmin(pmax(round((3.8 + .7 * l) * 3) / 3, 1), 5)),                    # descriptives not published: assumed
 S4 = list(N = 82, mu = 5.286, sd = .824, lo = 1, hi = 7, step = 0, r = .48, rho = .06, orgs = 4, peak = 6.69,
           yfun = function(l) pmin(pmax(5.145 + .589 * l, 1), 7)))
p <- S[[study]]
truths <- list("straight line" = function(x) x, "rises then flat above the mean" = function(x) pmin(x, p$mu), "paper's reading (performance peaks then falls)" = function(x) -(x - p$peak)^2)
allrows <- NULL
for (nm in names(truths)) {
  out <- NULL
  for (i in seq_len(reps)) {
    z <- rnorm(p$N); x <- p$mu + p$sd * z
    if (p$step > 0) x <- round(x / p$step) * p$step
    x <- pmin(pmax(x, p$lo), p$hi)
    t <- p$rho * z + sqrt(1 - p$rho^2) * rnorm(p$N)
    lat <- p$r * z01(truths[[nm]](x)) + sqrt(1 - p$r^2) * rnorm(p$N)
    X <- cbind(x, t)
    if (p$orgs > 0) { g <- sample(p$orgs, p$N, TRUE); lat <- lat + c(-.3, -.1, .1, .3)[g]
                      X <- cbind(X, sapply(2:p$orgs, function(k) as.numeric(g == k))) }
    y <- p$yfun(z01(lat))
    f <- krls(X = X, y = y, sigma = ncol(X), approx = "none", print.level = 0)
    pw <- f$derivatives[, 1]; hi <- x >= quantile(x, .975); lo <- x <= quantile(x, .025)
    l8 <- lowess(x, pw, f = .8, iter = 0); l5 <- lowess(x, pw, f = .5, iter = 0); l3 <- lowess(x, pw, f = .3, iter = 0)
    out <- rbind(out, c(top = mean(pw[hi]), bot = mean(pw[lo]), e8 = tail(l8$y, 1), e5 = tail(l5$y, 1), e3 = tail(l3$y, 1),
                        pk = mean(x < l8$x[which.max(l8$y)]), cor = cor(x, y), ntop = sum(hi), mid = mean(pw[x > quantile(x, .25) & x < quantile(x, .75)])))
  }
  allrows <- rbind(allrows, data.frame(truth = nm, run = seq_len(reps), out))
  cat(sprintf("%s | %s | r(x,y) %.2f | D=%d | top<0: %.0f%% | line<0 at top (.8/.5/.3): %.0f%% / %.0f%% / %.0f%% | bottom<0: %.0f%% | peak pct %.0f | mid effect %+.3f, top %+.3f | n top %d-%d | %d runs\n",
      study, nm, mean(out[, "cor"]), ncol(X), 100 * mean(out[, "top"] < 0), 100 * mean(out[, "e8"] < 0), 100 * mean(out[, "e5"] < 0), 100 * mean(out[, "e3"] < 0),
      100 * mean(out[, "bot"] < 0), 100 * median(out[, "pk"]), mean(out[, "mid"]), mean(out[, "top"]), min(out[, "ntop"]), max(out[, "ntop"]), reps))
}
write.csv(allrows, sprintf("grid05_perrun_%s_seed%s.csv", study, a[3]), row.names = FALSE)
