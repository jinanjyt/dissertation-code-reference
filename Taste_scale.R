## Taste: parameter scale
## initialisation
n <- 100
X <- matrix(seq(0, 10, length=n), ncol=1)
D <- distance(X)

C <- exp(-D) + diag(eps, n)

# generate plot
Y <- rmvnorm(10, sigma=C)
Y0 <- rmvnorm(10, sigma=0.01*C)
Y1 <- rmvnorm(10, sigma=25*C)
Y2 <- rmvnorm(10, sigma=100*C)

par(mfcol=c(4,1))
matplot(X, t(Y), type="l")
matplot(X, t(Y0), type="l")
matplot(X, t(Y1), type="l")
matplot(X, t(Y2), type="l")