######################################################################
## TRACHOMA MODEL
######################################################################

##Assumptions: see table - https://docs.google.com/document/d/1sX1YBmMvrjmM7WO1WN-K5C8d-ZzGGwCd/edit

##State variables
## S - Susceptible individuals (PCR negative and TF negative)
## I - Infectious individuals with clinical signs (PCR positive and TF positive)
## D - Chronic Injury 
## R - Recoveries

#Parameters
## N - the total population size; N = S + I + D + R
## lambda - Force of infection 
## gamma - Progression of infection
## rho - Progression to chronic infection
## alpha - Reinfection rate

# Note: For this model the population size is constant

library(deSolve)

## Model
sidr <- function(t,y,parms){
  with(c(as.list(y),parms),{
    dSdt <- -lambda*S
    dIdt <- lambda*S - gamma*I + alpha*lambda*R
    dDdt <- gamma*I - rho*D
    dRdt <- rho*D - alpha*lambda*R
    return(list(c(dSdt, dIdt, dDdt, dRdt)))
  })
}

## Initial values
N0 <- 100
init <- c(
  S =  0.7 * N0,  
  I =  0.1 * N0,  
  D =  0,  
  R =  0   
)

## Parameters
params <- c(
  lambda = 0.5,        
  gamma = 2, 
  rho = 1,
  alpha = 0.6,
  N = N0
)


## Time vector
times <- seq(0, 20, by = 1)


## Solve ODE
out <- as.data.frame(lsoda(
  y = init,
  times = times,
  func = sidr,
  parms = params
))

head(out)
subset(out, time == 20)


## Plot 
plot(out$time, out$I,
     type = "l",
     xlab = "Time in years",
     ylab = "Number Infected",
     main = "Trachoma Transmission Model",
     col = "red", lwd = 2,
     ylim = c(0, max(out$I, out$S, out$D, out$R))
)
lines(out$time, out$S, col = "blue", lwd = 2)
lines(out$time, out$D, col = "black", lwd = 2)
lines(out$time, out$R, col = "green4", lwd = 2)

legend("topright",
       legend = c("Susceptible (S)", "Infected (I)", "Diseased (D)", "Recovered (R)"),
       col = c("blue", "red", "black", "green4"),
       lwd = 2,
       bty = "n",     # No border
       inset = 0.02,
       cex = 0.9,
       y.intersp = 1.2,
       bg = "white"   # Background box for contrast
)
