## Simple 1D Gaussian Process Regression example
## use equation: y=x*sin(x)

## generate n data points as known data set (from equation assumed unknown)
n <- 10
x <- matrix(seq(-2*pi,2*pi,length=n), ncol=1)
y <- x*sin(x)

## calculate Euclidean distance
library(plgp)
d <- distance(x)

## use inverse exponential squared Euclidean distance as Covariance matrix
eps <- sqrt(.Machine$double.eps)
Sigma <- exp(-d) + diag(eps, ncol(d))


## define interested new data points interval
x_new <- matrix(seq(-0.5-2*pi, 2*pi+0.5, length=100), ncol=1)
d_new <- distance(x_new)
Sigma_new <- exp(-d_new) + diag(eps, ncol(d_new))

## calculate new covariance matrix
# for later matrix calculation, use distance(x_new, x) rather than distance(x, x_new)
D <- distance(x_new, x)
S <- exp(-D)

# calculate the distribution characteristics using equations
S_inverse <- solve(Sigma)
mup <- S %*% S_inverse %*% y
Sigmap <- Sigma_new - S %*% S_inverse %*% t(S)


## plot the figure to compare results
y_predict <- rmvnorm(100, mup, Sigmap)

q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))


matplot(x_new, t(y_predict), type="l", col="gray", lty=1, xlab="x", ylab="y")
# initial given data points
points(x, y, pch=20, cex=2)
# real function line
lines(x_new, x_new*sin(x_new), col="blue")
# mean
lines(x_new, mup, lwd=2)
# boundaries
lines(x_new, q1, lwd=2, lty=2, col=2)
lines(x_new, q2, lwd=2, lty=2, col=2)

