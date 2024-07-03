# FUNCS
#func.mle
func.mle <- function(x) 1/(1+mean(x))
#func.remove_geo_outliers
func.remove_geo_outliers <- function(x, a) {
  1 - (sum(x>=a & x<(2*a)) / sum(x<a))^(1/a)
}
#func.geo_trimmer
func.geo_trimmer <- function(x, a) {
  require(MASS)
  geom.tm <- function(p, xbar, cutoff) {(sum(dgeom(1:cutoff, p)*c(1:cutoff))-xbar)^2}
  cutoff <- quantile(x, a)
  xbar <- mean(x[x <= cutoff])
  optim(0.5, geom.tm, lower=1.0e-05, upper=1-1.0e-05, method="L-BFGS-B",
        xbar=xbar, c=cutoff)$par
}

# BUILD GEOMETRIC SAMPLING
# ?Geometric
# p(x|p)= p * (1−p)^x
# therefore, p(x+1|p)/p(x|p) = 1−p
a100 <- matrix(rgeom(100000, 0.3), 1000, 100)
b1000 <- matrix(rgeom(1000000, 0.3), 1000, 1000)
c10000 <- matrix(rgeom(1000000, 0.3), 1000, 10000)

#mle100 <- apply(a100, 1, func.remove_geo_outliers, a=3)
#mle100 <- apply(a100, 1, func.geo_trimmer, a=0.75)
mle100 <- apply(a100, 1, func.mle)
mle100 <- sort(mle100)
plot(mle100, pch=20, col="coral3", ylim=c(0.22, 0.41), log = "y")
# standard error of estimate = error sum of squares (SSE)
# sigma[i=1..n](() yi - mean(y) )^2)
mn = mean(mle100)
mle100diff = (mle100 - mn)^2
sse100 = sum(mle100diff)
cmpr100 = 1/(sqrt(100))

#mle1000 <- apply(b1000, 1, func.remove_geo_outliers, a=3)
#mle1000 <- apply(b1000, 1, func.geo_trimmer, a=0.75)
mle1000 <- apply(b1000, 1, func.mle)
mle1000 <- sort(mle1000)
par(new=T)
plot(mle1000, pch=20, col="blue", ylim=c(0.22, 0.41), log = "y")
mn = mean(mle1000)
mle1000diff = (mle1000 - mn)^2
sse1000 = sum(mle1000diff)
cmpr1000 = 1/(sqrt(1000))

#mle10000 <- apply(c10000, 1, func.remove_geo_outliers, a=3)
#mle10000 <- apply(c10000, 1, func.geo_trimmer, a=0.75)
mle10000 <- apply(c10000, 1, func.mle)
mle10000 <- sort(mle10000)
par(new=T)
plot(mle10000, pch=20, col="gold3", ylim=c(0.22, 0.41), log = "y")
mn = mean(mle10000)
mle10000diff = (mle10000 - mn)^2
sse10000 = sum(mle10000diff)
cmpr10000 = 1/(sqrt(10000))

t <- c(cmpr100, cmpr1000, cmpr10000)
s <- c(sse100, sse1000, sse10000)
plot(t, pch=10, col="seagreen", ylim=c(0, 0.7))
abline(lsfit(1:3,t))
par(new=T)
plot(s, pch=20, col="mediumpurple3", ylim=c(0, 0.7))
abline(lsfit(1:3,s))
