# Size of the implied fall at the top and bottom of the footnote-4 smoothed line (lowess f = .8 of the pointwise
# derivatives), for the stand-ins of Studies 1 to 5. Data generation and KRLS settings as study1_official_perrun.R
# (authors' package defaults: sigma = number of columns, exact, LOO lambda). Fall = area under zero of the smoothed
# line beyond its last zero crossing (top) or before its first (bottom), in outcome-scale points.
suppressMessages(library(KRLS))
a <- commandArgs(trailingOnly = TRUE); study <- a[1]; reps <- as.integer(a[2]); set.seed(as.integer(a[3]))
z01 <- function(v) (v - mean(v)) / sd(v)
S <- list(
 S1 = list(N = 474, mu = 3.335, sd = .751, lo = 1, hi = 5, step = 1/7, r = .26, rho = -.17, orgs = 0, peak = 4.75,
           yfun = function(l) pmin(pmax(round((3.695 + .752 * l) * 4) / 4, 1), 5)),
 S2 = list(N = 301, mu = 4.020, sd = 1.306, lo = 1, hi = 7, step = 1/7, r = .19, rho = .01, orgs = 0, peak = 6.3,
           yfun = function(l) as.numeric(cut(rank(l) / length(l), c(0, .2, .4, .6, .8, 1)))),
 S5 = list(N = 800, mu = 3.390, sd = .425, lo = 1, hi = 5, step = 0, r = .18, rho = -.08, orgs = 0, peak = 4.09,
           yfun = function(l) { s <- exp(-.6 * l); pmin(10, 10 - .88 * s / mean(s)) }),
 S3 = list(N = 365, mu = 3.35, sd = .30, lo = 1, hi = 5, step = 1/7, r = .15, rho = 0, orgs = 4, peak = 3.7,
           yfun = function(l) pmin(pmax(round((3.8 + .7 * l) * 3) / 3, 1), 5)),
 S4 = list(N = 82, mu = 5.286, sd = .824, lo = 1, hi = 7, step = 0, r = .48, rho = .06, orgs = 4, peak = 6.69,
           yfun = function(l) pmin(pmax(5.145 + .589 * l, 1), 7)))
p <- S[[study]]
truths <- list("straight line" = function(x) x, "rises then flat above the mean" = function(x) pmin(x, p$mu),
               "paper's reading (peaks then falls)" = function(x) -(x - p$peak)^2)
integ <- function(xx, yy) sum(diff(xx) * (head(yy, -1) + tail(yy, -1)) / 2)
falls <- function(lx, ly) {
  k <- !duplicated(lx); lx <- lx[k]; ly <- ly[k]; n <- length(ly)
  top <- 0; if (ly[n] < 0) { j <- max(which(ly >= 0), 0); i0 <- max(j, 1); top <- -integ(lx[i0:n], pmin(ly[i0:n], 0)) }
  bot <- 0; if (ly[1] < 0) { j <- min(c(which(ly >= 0), n)); bot <- -integ(lx[1:j], pmin(ly[1:j], 0)) }
  c(top = top, bot = bot, rise = integ(lx, pmax(ly, 0)))
}
rows <- NULL
for (nm in names(truths)) for (i in seq_len(reps)) {
  z <- rnorm(p$N); x <- p$mu + p$sd * z
  if (p$step > 0) x <- round(x / p$step) * p$step
  x <- pmin(pmax(x, p$lo), p$hi)
  t <- p$rho * z + sqrt(1 - p$rho^2) * rnorm(p$N); if (Sys.getenv("SKEWT") == "1") t <- qgamma(pnorm(t), shape = (5.504/7.306)^2, scale = 7.306^2/5.504)
  lat <- p$r * z01(truths[[nm]](x)) + sqrt(1 - p$r^2) * rnorm(p$N)
  X <- cbind(x, t)
  if (p$orgs > 0) { g <- sample(p$orgs, p$N, TRUE); lat <- lat + c(-.3, -.1, .1, .3)[g]
                    X <- cbind(X, sapply(2:p$orgs, function(k) as.numeric(g == k))) }
  y <- p$yfun(z01(lat))
  f <- krls(X = X, y = y, sigma = ncol(X), approx = "none", print.level = 0)
  l8 <- lowess(x, f$derivatives[, 1], f = .8, iter = 0)
  rows <- rbind(rows, data.frame(study = study, truth = nm, run = i, t(falls(l8$x, l8$y)), end = tail(l8$y, 1)))
}
write.csv(rows, sprintf("fallsize_%s_skewt%s_seed%s.csv", study, Sys.getenv("SKEWT"), a[3]), row.names = FALSE)
