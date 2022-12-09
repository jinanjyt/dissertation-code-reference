### use SIR model to generate Infected number of a particular time under a few initial guesses of (beta, gamma)

## for example, Infected number after one week (day 7);
nday = 7;

## (beta, gamma), say at the beginning, gamma is very low, assume zero;
## beta choice: 0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5;
beta_trail = c(0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
n = length(beta_trail)

## Set parameters (not changed the population proportion in the original demo code)
# Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init <- c(S = 1 - 1e-6 , I = 1e-6, R = 0.0)

# Time frame
times      <- seq(0, nday, by = 1)

## call SIR model (similar)
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

## trails
library(deSolve)

# initialise an empty vector to store results of infected people number at day 7
day7 <- rep(0,n)

# calculate the infected proportion at day 7
for (i in 1:n) {
  parameters = c(beta = beta_trail[i], gamma = 0)
  
  ## Solve using ode (General Solver for Ordinary Differential Equations)
  out <- ode(y = init, times = times, func = sir, parms = parameters)
  out <- as.data.frame(out)
  day7[i] <- out[nday+1,3]
}

# show results
print(day7)

### apply GPR to SIR model

library(plgp)

## 1d since gamma is zero here
# take the input: x and y
x <- beta_trail
y <- day7

# calculate distance and hence covariance matrix
d <- distance(x)
eps <- sqrt(.Machine$double.eps)
Sigma <- exp(-d) + diag(eps, ncol(d))


## set new coordinate interval of beta where you wanted to explore
x_new <- matrix(seq(beta_trail[1]-0.5, beta_trail[n]+0.5, length=100), ncol=1)
explore_length=length(x_new)

library(lhs)

## apply GPR
d_new <- distance(x_new)
Sigma_new <- exp(-d_new) + diag(eps, ncol(d_new))

DX <- distance(x_new, x)
SX <- exp(-DX)

Si <- solve(Sigma)
mup <- SX %*% Si %*% y
Sigmap <- Sigma_new - SX %*% Si %*% t(SX)

## generate predictions
y_predict <- rmvnorm(explore_length, mup, Sigmap)

# confidence interval
q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))


## generate the plot of GPR-SIR predictions
matplot(x_new, t(y_predict), type="l", col="gray", lty=1, xlab="beta", ylab="Infected (I_7)")
points(x, y, pch=20, cex=2)
lines(x_new, mup, lwd=2)
lines(x_new, q1, lwd=2, lty=2, col=2)
lines(x_new, q2, lwd=2, lty=2, col=2)
