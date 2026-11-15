suppressMessages(library(KRLS)); set.seed(7)
z01 <- function(v) (v - mean(v)) / sd(v)
for (r in c(.15, .26, .5)) { res <- NULL
 for (i in 1:20) { N<-474; z<-rnorm(N); x<-pmin(pmax(round((3.335+.751*z)*7)/7,1),5); t<- -.17*z+sqrt(1-.17^2)*rnorm(N)
  lat<-r*z01(x)+sqrt(1-r^2)*rnorm(N); y<-pmin(pmax(round((3.695+.752*z01(lat))*4)/4,1),5)
  f<-krls(X=cbind(x,t),y=y,sigma=2,print.level=0); pw<-f$derivatives[,1]
  hi<-x>=quantile(x,.975); mid<-x>quantile(x,.25)&x<quantile(x,.75)
  res<-rbind(res,c(lambda=f$lambda, ratio=mean(pw[hi])/mean(pw[mid]), top_neg=mean(pw[hi])<0)) }
 cat(sprintf("r=%.2f  median lambda %.3f  median top/mid slope ratio %.2f  top negative %.0f%%\n", r, median(res[,1]), median(res[,2]), 100*mean(res[,3]))) }
# Output, 15 Nov 2026 (seed 7, 20 datasets per level):
# r=0.15  median lambda 45.946  median top/mid slope ratio -0.38  top negative 95%
# r=0.26  median lambda 10.982  median top/mid slope ratio -0.28  top negative 75%
# r=0.50  median lambda 2.151  median top/mid slope ratio -0.21  top negative 75%
