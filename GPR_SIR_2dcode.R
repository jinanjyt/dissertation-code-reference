### use SIR model to generate Infected number of a particular time under a few initial guesses of (beta, gamma)

## for example, Infected number after one week (day 7);
nday = 7;

## Set parameters (I have not changed the population proportion in the original demo code)
# Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init <- c(S = 1 - 1e-6 , I = 1e-6, R = 0.0)

# beta: infection parameter; gamma: recovery parameter
# Here, beta and gamma are not single data parameter, we want a range of beta and gamma
beta_start = 0;
beta_end = 3;
gamma_start = 0;
gamma_end = 2;
# set the number of data points wanted
n_beta <- 10
n_gamma <- 10

# generate pair of (beta, gamma) and gain all possible combinations
beta_seq <- seq(beta_start, beta_end, length.out = n_beta)
gamma_seq <- seq(gamma_start, gamma_end, length.out = n_gamma)
generate_points <- expand.grid(x = beta_seq, y = gamma_seq)
# Time frame
times      <- seq(0, nday, by = 1)

## call SIR model
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

# to store results of infected people number at day 7
nn = n_beta*n_gamma
day7 <- rep(0, nn)

# compute results at day 7
for (i in 1:nn) {
   parameters = c(beta = generate_points[i,1], gamma = generate_points[i,2])
   
   ## Solve using ode (General Solver for Ordinary Differential Equations)
   out <- ode(y = init, times = times, func = sir, parms = parameters)
   out <- as.data.frame(out)
   day7[i] <- out[nday+1, 3]
}

print(day7)

### use GPR to generate a plot of first few trails of (beta, gamma) of SIR model

library(plgp)

X <- generate_points
y <- day7 # later used in mup
D <- distance(X)
eps <- sqrt(.Machine$double.eps)
Sigma <- exp(-D) + diag(eps, ncol(D))



# 2d: create matrix of where you wanted to explore beta, gamma
library(lhs)
explore_length <- 20

X1 <- seq(beta_start, beta_end, length = explore_length)
X2 <- seq(gamma_start, gamma_end, length = explore_length)
XX <- expand.grid(X1, X2)

DXX <- distance(XX)
SXX <- exp(-DXX) + diag(0.001, ncol(DXX))

DX <- distance(XX, X)
SX <- exp(-DX)

Si <- solve(Sigma)
mup <- SX %*% Si %*% y
Sigmap <- SXX - SX %*% Si %*% t(SX)


q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))



## 2d (beta, gamma)~I when gamma not equal to zero
par(mfcol=c(2,1))
persp(X1, X2, matrix(mup, ncol = 20), theta = -60, phi = 20, xlab = "beta", ylab = "gamma", zlab = "Infected", col = "orange", shade = 0.6)
persp(X1, X2, matrix(mup, ncol = 20), theta = 30, phi = 20, xlab = "beta", ylab = "gamma", zlab = "Infected", col = "orange", shade = 0.6)


image(X1, X2, matrix(mup, ncol = 20), xlab = "beta", ylab = "gamma", main = "Infected (mean)")
points(generate_points, pch=20)
#image(X1, X2, matrix(sqrt(diag(Sigmap)), ncol = 20), xlab = "beta", ylab = "gamma", main = "Infected (standard deviation)")
#points(generate_points, pch=20)

