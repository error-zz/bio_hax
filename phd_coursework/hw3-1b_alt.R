#Mle
func.mle <- function(x) 1/(1+mean(x))
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

X <- rgeom(100000, 0.3)

a100 <- matrix(sample(x = X, size = 100, replace = TRUE), 1000, 100)
mle100 <- apply(a100, 1, func.remove_geo_outliers, a=3)
mle100 <- apply(a100, 1, func.geo_trimmer, a=0.75)
mle100 <- apply(a100, 1, func.mle)
plot(mle100, pch=20, col="coral3")
mn = mean(mle100)
mle100diff = (mle100 - mn)^2
sse100 = sum(mle100diff)
mse100 = var(mle100)-(sse100/(100-2))
cmpr100 = 1/(sqrt(100))

b1000 <- matrix(sample(x = X, size = 1000, replace = TRUE), 1000, 1000)
mle1000 <- apply(b1000, 1, func.remove_geo_outliers, a=3)
mle1000 <- apply(b1000, 1, func.geo_trimmer, a=0.75)
mle1000 <- apply(b1000, 1, func.mle)
mle1000 <- sort(mle1000)
par(new=T)
plot(mle1000, pch=20, col="blue", ylim=c(0.2, 0.5), log = "y")
mn = mean(mle1000)
mle1000diff = (mle1000 - mn)^2
sse1000 = sum(mle1000diff)
mse1000 = var(mle1000)-(0/(1000-2))
cmpr1000 = 1/(sqrt(1000))

c10000 <- matrix(sample(x = X, size = 10000, replace = TRUE), 1000, 100000)
mle10000 <- apply(c10000, 1, func.remove_geo_outliers, a=3)
mle10000 <- apply(c10000, 1, func.geo_trimmer, a=0.75)
mle10000 <- apply(c10000, 1, func.mle)
mle10000 <- sort(mle10000)
par(new=T)
plot(mle10000, pch=20, col="gold3", ylim=c(0.2, 0.5), log = "y")
mn = mean(mle10000)
mle10000diff = (mle10000 - mn)^2
sse10000 = sum(mle10000diff)
mse10000 = var(mle10000)-(0/(10000-2))
cmpr10000 = 1/(sqrt(10000))

t <- c(cmpr100, cmpr1000, cmpr10000)
s <- c(mse100, mse1000, mse10000)
plot(t, pch=10, col="seagreen", ylim=c(-1, 0.12))
abline(lsfit(1:3,t))
par(new=T)
plot(s, pch=20, col="mediumpurple3", ylim=c(-1, 0.12))
abline(lsfit(1:3,s))