#### Basic SIR model code

## Load "deSolve" package
## If not installed: use install.packages('deSolve')
library(deSolve)

## Create SIR function
## input parameters (time, state, parameters); return (dS, dI, dR)
sir <- function(time, state, parameters){
  
  with(as.list(c(state, parameters)),{
    
    dS <- -beta * S * I
    dI <- beta * S * I -gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

## Example:
## Set input parameters (time, state, parameters)
## Time frame
times <- seq(0, 70, by = 1)
## Initial state
init <- c(S = 1-1e-6, I = 1e-6, R = 0)
## Two parameters: beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 1.42, gamma = 0.14)

## Solve using ode
output <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
output <- as.data.frame(output)
## Delete time variable
output$time <- NULL

## Plot
matplot(x = times, y = output, type = 'l', xlab = 'Time', 
        ylab = 'Susceptible, Infected and Recovered', main = 'SIR Model',
        lwd = 1, lty = 1, bty = 'l', col = 2:4)
## Add legend
legend(50, 0.6, c('Susceptible', 'Infected', 'Recovered'), pch = 1, col = 2:4,
       bty = 'n')

