# Size of the implied fall at the top (and bottom) of the footnote-4 smoothed line, Study 1 stand-in.
# Same data generation and KRLS settings as study1_official_perrun.R; adds the integral of the smoothed line.
suppressMessages(library(KRLS))
a <- commandArgs(trailingOnly = TRUE); reps <- as.integer(a[1]); set.seed(as.integer(a[2]))
z01 <- function(v) (v - mean(v)) / sd(v)
N<-474; mu<-3.335; s<-.751; step<-1/7; r<-.26; rho<- -.17; peak<-4.75
yfun <- function(l) pmin(pmax(round((3.695 + .752 * l) * 4) / 4, 1), 5)
truths <- list("straight line"=function(x) x, "rises then flat above the mean"=function(x) pmin(x, mu),
               "paper's reading (peaks then falls)"=function(x) -(x-peak)^2)
integ <- function(xx, yy) sum(diff(xx) * (head(yy,-1) + tail(yy,-1)) / 2)
falls <- function(lx, ly) {
  k <- !duplicated(lx); lx <- lx[k]; ly <- ly[k]; n <- length(ly)
  top <- 0; if (ly[n] < 0) { j <- max(which(ly >= 0), 0); i0 <- max(j, 1); top <- -integ(lx[i0:n], pmin(ly[i0:n], 0)) }
  bot <- 0; if (ly[1] < 0) { j <- min(c(which(ly >= 0), n)); bot <- -integ(lx[1:j], pmin(ly[1:j], 0)) }
  rise <- integ(lx, pmax(ly, 0)); c(top = top, bot = bot, rise = rise)
}
rows <- NULL
for (nm in names(truths)) for (i in seq_len(reps)) {
  z <- rnorm(N); x <- pmin(pmax(round((mu + s*z)/step)*step, 1), 5)
  t <- rho*z + sqrt(1-rho^2)*rnorm(N)
  lat <- r*z01(truths[[nm]](x)) + sqrt(1-r^2)*rnorm(N); y <- yfun(z01(lat))
  f <- krls(X = cbind(x, t), y = y, sigma = 2, approx = "none", print.level = 0)
  l8 <- lowess(x, f$derivatives[,1], f = .8, iter = 0)
  rows <- rbind(rows, data.frame(truth = nm, run = i, t(falls(l8$x, l8$y)), end = tail(l8$y,1), xmax = max(x)))
}
write.csv(rows, sprintf("fallsize_seed%s.csv", a[2]), row.names = FALSE)
