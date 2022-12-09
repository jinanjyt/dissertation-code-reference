## Simple 1D Gaussian Process Regression example
## use equation: y=x*sin(x)
## scale parameter added
## nugget parameter added

## generate n data points as known data set (from equation assumed unknown)
n <- 10
x <- matrix(seq(-2*pi,2*pi,length=n), ncol=1)
# add noise
y <- x*sin(x) + rnorm(n, sd=1)


## calculate Euclidean distance
library(plgp)
d <- distance(x)

##-------------------------------------------------------------------
## nugget part

# function creation in Surrogates book
nlg <- function(g, D, Y)
{
  n <- length(Y)
  K <- exp(-D) + diag(g, n)
  Ki <- solve(K)
  # gain logarithm
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  # calculate negative log-likelihood for later minimisation
  ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
  counter <<- counter + 1
  return(-ll)
}


# initialise counter
counter <- 0
eps <- sqrt(.Machine$double.eps)
# method: maximum log-likelihood
#       = minimising negative log-likelihood using optimize build-in function
g <- optimize(nlg, interval=c(eps, var(y)), D=d, Y=y)$minimum
g


## end of nugget part
##-------------------------------------------------------------------
## use inverse exponential squared Euclidean distance as Covariance matrix
Sigma <- exp(-d) + diag(g, ncol(d)) # here use g instead of eps


## define interested new data points interval
x_new <- matrix(seq(-0.5-2*pi, 2*pi+0.5, length=100), ncol=1)
d_new <- distance(x_new)
Sigma_new <- exp(-d_new) + diag(g, ncol(d_new)) # again use g instead of eps

## calculate new covariance matrix
# for later matrix calculation, use distance(x_new, x) 
# rather than distance(x, x_new)
D <- distance(x_new, x)
S <- exp(-D)

# scale parameter-------------------------------------------------------------
S_inverse <- solve(Sigma)
tau2hat <- drop(t(y) %*% S_inverse %*% y / length(y))


## kernel function with scale and nugget parameters----------------------------
mup <- S %*% S_inverse %*% y
Sigmap <- tau2hat*(Sigma_new - S %*% S_inverse %*% t(S))


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
