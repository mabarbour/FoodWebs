######### These are my attempts to incorporate McCann's models into R to help me interpret these food web models ######################

setwd("~/Documents/FoodWebs")

### Section 2.2.5
# load libraries
library("deSolve")

# parameters and state variables for the R-M model
r <- 1.0 # per capita rate of increase in resource
K <- 2.0 # indicated as 1.0 in Fig 2.1 description, but I think this was a typo in the book.
e <- 0.5
Ro <- 0.5
m <- 0.5
a <- 1.3 # initial attack rate of 1.3 instead of 1.2 more closely replicates the figures in the book
#R <- seq(0.1,2,0.1) # necessary for creating the x-axis for Fig 2.1
  
# rmcr = Rosenzweig-MacArthur (R-M) consumer-resource (C-R) model that assumes logistic resource growth and a type 2 functional response by the consumer (i.e. consumption rate saturates with resource density). This function is necessary for the ordinary differential equation solver, ode() 
rmcr <- function(t,y,p) {
  # t,y,p is the necessary format for solving using the ode() function in R
  R <- y[1]
  C <- y[2]
  with(as.list(p), {
    dR.dt <- r * R * (1 - R / K) - a * C * R / (R + Ro)
    dC.dt <- e * a * C * R / (R + Ro) - m * C
    return(list(c(dR.dt,dC.dt)))
  })
}

#### Modified Figures 2.1 and 2.6 to show how increasing attack rate influences consumer and resource isoclines as well as consumer and resource densities over time. 

pdf("AttackRate.C-R.pdf")
par(mfrow=c(1,2)) # plot graphs side-by-side

# initial state and parameter values for experiment
i.state <- c(R=0.6,C=0.1) # initial state variables 
p.rm1 <- c(r = r, e = e, a = a, K = K, Ro = Ro, m=m) # initial parameters 
Time <- 300 # set time scale
rm1 <- ode(i.state,1:Time, rmcr, p.rm1) # run the experiment

# resource isocline
Rx <- seq(0.1,2,0.1) # manipulating different Resource densities to solve R isocline.
Riso <- expression(r/a * (Rx + Ro) * (1 - Rx/K)) # set R = 0 and solved algebraically
RisoStable <- eval(Riso)

# consumer isocline
Ciso <- expression(m * Ro / (e * a - m)) # set C = 0, and solved algebraically.
CisoStable <- eval(Ciso)

## Experiment 1: a = 1.3
# plot densities
matplot(rm1[,1],rm1[,c(2,3)], type = "l", ylab="Density", xlab="Time", ylim=c(0,2), main = "a = 1.3")
legend("right", c("R","C"), lty=1:2, col=1:2, bty="n")

# plot stability around consumer and resource isoclines
plot(Rx,RisoStable, type = "l", ylab = "C", ylim=0:1, xlim = c(0,2))
abline(v=CisoStable, lty = 2, col =2) # abline for consumer isocline was determined using pen and paper. The alternative isocline is C = 0, which is uninformative.
legend("topleft", c("R-isocline","C-isocline"), lty=1:2, bty="n", cex=0.8, col=1:2)
points(i.state[1],i.state[2]) # starting point of experiment
# points(x=Riso) # need to find equilibrium
arrows(rm1[-Time,2], rm1[-Time,3], rm1[-1,2], rm1[-1,3], length=0.1, lty=1) # trace stability across different time steps. Don't know why there are so many warnings

### Experiment 2: a = 1.6
# adjust parameters
a <- 1.6
p.rm2 <- c(r = r, e = e, a = 1.6, K = K, Ro = Ro)
rm2 <- ode(i.state,1:Time, rmcr, p.rm2)

# resource isocline
Riso <- expression(r/a * (Rx + Ro) * (1 - Rx/K)) # adjusted attack rate
RisoStable <- eval(Riso)

# consumer isocline
Ciso <- expression(m * Ro / (e * a - m)) # adjusted attack rate
CisoStable <- eval(Ciso)

# plot consumer and resource densities
matplot(rm2[,1],rm2[,c(2,3)], type = "l", ylab="Density", xlab="Time", ylim=c(0,2), main = "a = 1.6")
legend("topright", c("R","C"), lty=1:2, col=1:2, bty="n")

# plot stabilities around consumer and resource isoclines. Why is it below the equilibrium?
plot(Rx, RisoStable, type = "l", ylab = "C", ylim=0:1)
abline(v=m * Ro / (e * 1.6 - m), lty = 2, col =2) # attack rate adjusted to a = 1.6
legend("topright", c("R-isocline","C-isocline"), lty=1:2, bty="n", cex=0.8, col=1:2)
points(i.state[1],i.state[2]) # starting point of experiment
arrows(rm2[-Time,2], rm2[-Time,3], rm2[-1,2], rm2[-1,3], length=0.1, lty=1) 

### Experiment 3: a = 2.0
# adjust parameter values
a <- 2.0
p.rm3 <- c(r = r, e = e, a = 2.0, K = K, Ro = Ro)
Time <- 50 # same pattern at larger time values, but I adjusted it to see the pattern more clearly
rm3 <- ode(i.state,1:Time, rmcr, p.rm3)

# resource isocline
Riso <- expression(r/a * (Rx + Ro) * (1 - Rx/K)) # adjusted attack rate
RisoStable <- eval(Riso)

# consumer isocline
Ciso <- expression(m * Ro / (e * a - m)) # adjusted attack rate
CisoStable <- eval(Ciso)

# plot densities
matplot(rm3[,1],rm3[,c(2,3)], type = "l", ylab="Density", xlab="Time", ylim=c(0,2), main = "a = 2.0")
legend("topright", c("R","C"), lty=1:2, col=1:2, bty="n", cex=0.8)

# plot stabilities around resouce and consumer isoclines 
plot(Rx,RisoStable, type = "l", ylab = "C", ylim=0:1)
abline(v=CisoStable, lty = 2, col = 2) # attack rate adjusted to a = 2.0
legend("topright", c("R-isocline","C-isocline"), lty=1:2, col=1:2, bty="n", cex=0.8)
points(i.state[1],i.state[2]) # starting point of experiment
arrows(rm3[-Time,2], rm3[-Time,3], rm3[-1,2], rm3[-1,3], length=0.1, lty=1) 

## Figure 2.7a: Bifurcation plot  (All credit goes to this website: http://www.r-bloggers.com/r-tools-for-dynamical-systems-bifurcation-plot-in-r%C2%A0for%C2%A0system%C2%A0of%C2%A0odes/)

param.name <- "a" # choose parameter to perturb
param.seq <- seq(1,2.5,length=50) # choose range of parameters
p.rm4 <- c(r = r, e = e, a = 1.0, K = K, Ro = Ro, m=m) # set starting parameters.
param.index <- which(param.name == names(p.rm4)) # tells the loop which parameter in "p.rm4" to grab for manipulation.
Time <- 50 # decreased time to see outcomes of simulations

out <- list()
for (i in 1:length(param.seq))
  out[[i]] <- matrix(0, Time, length(i.state))

for (i in 1:length(param.seq)) {
  # set params
  p.rm4.loop <- p.rm4
  p.rm4.loop[param.index] <- param.seq[i] # changes the parameter value for manipulation in the "init" function below.
  # converge
  init <- ode(i.state, 1:Time, rmcr, p.rm4.loop)
  # get converged points
  out[[i]] <- ode(init[Time,-1], 1:Time, rmcr, p.rm4.loop)[,-1]
}

range.lim <- lapply(out, function(x) apply(x, 2, range)) # don't completely understand what is going on here...
#range.lim <- apply(do.call("rbind", range.lim), 2, range)

# consumer bifurcation plot
plot.variable <- "C" # choose which variable to show
plot(0, 0, pch = "", xlab = param.name, ylab = plot.variable, xlim = range(param.seq), ylim = c(0,1), main = "Bifurcation plot") # range.lim[,plot.variable]
for (i in 1:length(param.seq)) {
  points(rep(param.seq[i], Time), out[[i]][,plot.variable])
}

# resource bifurcation plot
plot.variable <- "R" # choose which variable to show
plot(0, 0, pch = "", xlab = param.name, ylab = plot.variable, xlim = range(param.seq), ylim = c(0,2), main = "Bifurcation plot") # ylim = range.lim[,plot.variable]
for (i in 1:length(param.seq)) {
  points(rep(param.seq[i], Time), out[[i]][,plot.variable])
}

dev.off()

### Figure 2.7b: Plot of eigen values. Tips taken from Enemy-Victim Interactions chapter of a Primer of Ecology in R

# !!!!!!!!!! Currently not getting the dip in eigen value like I expected
dR.dt <- expression(r * R * (1 - R / K) - a * C * R / (R + Ro))
dC.dt <- expression(e * a * C * R / (R + Ro) - m * C)
RMjac1 <- list(D(dR.dt,"R"),D(dC.dt,"R"),D(dR.dt,"C"),D(dC.dt,"C"))

# confused about what I'm doing here...
#R <- m * Ro / (e * a - m) # expression is the consumer isocline... but why?
#C <- eval(Riso)

RM.jac2 <- matrix(sapply(RMjac1, function(pd) eval(pd)), nrow = 2)
eigen(RM.jac2)[["values"]]

#dev.off() # turn off pdf function

rmcr.jacList <- list()
for (i in 1:length(param.seq)) {
  # set params
  p.rm4.loop <- p.rm4
  p.rm4.loop[param.index] <- param.seq[i] # changes the parameter value for manipulation in the "init" function below.
  # converge
  init <- ode(i.state, 1:Time, rmcr, p.rm4.loop)
  # get converged points
  rmcr.jacList[[i]] <- ode(init[Time,-1], 1:Time, rmcr, p.rm4.loop)[,-1]
}

eigen.lim <- lapply(rmcr.jacList, function(x) apply(x, 2, range))
               
L1 <- sapply(rmcr.jacList, function(J) max(Re(eigen(J)[["values"]])))

plot(param.seq,L1,type ="l",xlab = "a", ylim=c(-0.4,0.3))

