suppressMessages(library(KRLS))
for (nm in c("line", "quad")) {
  d <- read.csv(sprintf("rcheck/%s.csv", nm)); X <- as.matrix(d[, c("safety", "tenure")])
  f <- krls(X = X, y = d$y, sigma = 2, approx = "none", print.level = 0)
  write.csv(data.frame(d_safety = f$derivatives[, 1], fitted = as.vector(f$fitted)), sprintf("rcheck/%s_old.csv", nm), row.names = FALSE)
  writeLines(sprintf("%.10f", f$lambda), sprintf("rcheck/%s_old.lambda", nm))
  cat(nm, "official, sigma = 2: lambda", f$lambda, "\n")
}
