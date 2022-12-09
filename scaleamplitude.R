# Parameter: Scale

## Simple 1D Gaussian Process Regression example with amplitude change


## Case1: 0.1*x*sin(X)

## training data points
n <- 10
X <- matrix(seq(-2*pi, 2*pi, length=n), ncol=1)
y <- 0.1*X*sin(X)

## GPR
D <- distance(X)
Sigma <- exp(-D)
XX <- matrix(seq(-0.5-2*pi, 2*pi + 0.5, length=100), ncol=1)
DXX <- distance(XX)
eps <- sqrt(.Machine$double.eps)
SXX <- exp(-DXX) + diag(eps, ncol(DXX))
DX <- distance(XX, X)
SX <- exp(-DX)
Si <- solve(Sigma);
mup <- SX %*% Si %*% y
Sigmap <- SXX - SX %*% Si %*% t(SX)

YY <- rmvnorm(100, mup, Sigmap)
q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))

par(mfrow=c(3,2))
matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y=0.1*x*sin(x)")
points(X, y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, 0.1*XX*sin(XX), col="blue")
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)


## estimated scale
CX <- SX
Ci <- Si
CXX <- SXX
# estimated choice of scale
tau2hat <- drop(t(y) %*% Ci %*% y / length(y))

2*sqrt(tau2hat)

mup2 <- CX %*% Ci %*% y
# difference: multiply tau2hat before old Sigmap (scaled)
Sigmap2 <- tau2hat*(CXX - CX %*% Ci %*% t(CX))

# generate plot under new covariance matrix: Sigmap2 (scaled)
YY <- rmvnorm(100, mup2, Sigmap2)
q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap2)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap2)))

matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y=0.1*x*sin(x)")
points(X, y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, 0.1*XX*sin(XX), col="blue")
lines(XX, q1, lwd=2, lty=2, col=2); lines(XX, q2, lwd=2, lty=2, col=2)

## Case 2: 4*x*sin(X)

## training data points
n <- 10
X <- matrix(seq(-2*pi, 2*pi, length=n), ncol=1)
y <- 4*X*sin(X)

## GPR
D <- distance(X)
Sigma <- exp(-D)
XX <- matrix(seq(-0.5-2*pi, 2*pi + 0.5, length=100), ncol=1)
DXX <- distance(XX)
SXX <- exp(-DXX) + diag(eps, ncol(DXX))
DX <- distance(XX, X)
SX <- exp(-DX)
Si <- solve(Sigma);
mup <- SX %*% Si %*% y
Sigmap <- SXX - SX %*% Si %*% t(SX)

YY <- rmvnorm(100, mup, Sigmap)
q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))
matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y=4*x*sin(x)")
points(X, y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, 4*XX*sin(XX), col="blue")
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)


## estimated scale
CX <- SX
Ci <- Si
CXX <- SXX
# estimated choice of scale
tau2hat <- drop(t(y) %*% Ci %*% y / length(y))

2*sqrt(tau2hat)

mup2 <- CX %*% Ci %*% y
# difference: multiply tau2hat before old Sigmap (scaled)
Sigmap2 <- tau2hat*(CXX - CX %*% Ci %*% t(CX))

# generate plot under new covariance matrix: Sigmap2 (scaled)
YY <- rmvnorm(100, mup2, Sigmap2)
q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap2)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap2)))


matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y=4*x*sin(x)")
points(X, y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, 4*XX*sin(XX), col="blue")
lines(XX, q1, lwd=2, lty=2, col=2); lines(XX, q2, lwd=2, lty=2, col=2)



## Case 3: x*sin(X)

## training data points
n <- 10
X <- matrix(seq(-2*pi, 2*pi, length=n), ncol=1)
y <- X*sin(X)

## GPR
D <- distance(X)
Sigma <- exp(-D)
XX <- matrix(seq(-0.5-2*pi, 2*pi + 0.5, length=100), ncol=1)
DXX <- distance(XX)
SXX <- exp(-DXX) + diag(eps, ncol(DXX))
DX <- distance(XX, X)
SX <- exp(-DX)
Si <- solve(Sigma);
mup <- SX %*% Si %*% y
Sigmap <- SXX - SX %*% Si %*% t(SX)

YY <- rmvnorm(100, mup, Sigmap)
q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))

#par(mfcol=c(2,1))
matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y=x*sin(x)")
points(X, y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, XX*sin(XX), col="blue")
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)


## estimated scale
CX <- SX
Ci <- Si
CXX <- SXX
# estimated choice of scale
tau2hat <- drop(t(y) %*% Ci %*% y / length(y))

2*sqrt(tau2hat)

mup2 <- CX %*% Ci %*% y
# difference: multiply tau2hat before old Sigmap (scaled)
Sigmap2 <- tau2hat*(CXX - CX %*% Ci %*% t(CX))

# generate plot under new covariance matrix: Sigmap2 (scaled)
YY <- rmvnorm(100, mup2, Sigmap2)
q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap2)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap2)))

matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y=x*sin(x)")
points(X, y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, XX*sin(XX), col="blue")
lines(XX, q1, lwd=2, lty=2, col=2); lines(XX, q2, lwd=2, lty=2, col=2)






