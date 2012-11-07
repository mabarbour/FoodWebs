######### These are my attempts to incorporate McCann's models into R to help me interpret these food web models ######################

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
R <- seq(0.1,2,0.1) # necessary for creating the x-axis for Fig 2.1
  
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

Riso <- expression(r/a * (R + Ro) * (1 - R/K)) # solved for isocline separately using pen and paper
RisoStable <- eval(Riso)

## Figure 2.1a !!!!! Consumer isocline not drawing !!!
plot(R,RisoStable, type = "l", ylab = "C", ylim=0:1)
abline(v=m * Ro / (e * a - m), lty = 2) # abline for consumer isocline was determined using pen and paper. The alternative isocline is C = 0, whichi is uninformative.
# arrows(RM1[-Time,2], RM1[-Time,3], RM1[-1,2], RM1[-1,3], length=0.1); not accurate now but a potential technique to trace stability (from Primer for Ecology with R)

## Figure 2.6. Examine how increasing attack rate influences stability.
# create vector of parameter values for "rmcr" function
p.rm1 <- c(r = r, e = e, a = 1.3, K = K, Ro = Ro)
p.rm2 <- c(r = r, e = e, a = 1.6, K = K, Ro = Ro)
p.rm3 <- c(r = r, e = e, a = 2.0, K = K, Ro = Ro)

Time <- 300 # set time scale

# create ordinary differential equations with increasing attack rates (all other parameter values held constant).  Initial state values (consumer and resource densities) are both set at 0.1
i.state <- c(R=0.1,C=0.1)
rm1 <- ode(i.state,1:Time, rmcr, p.rm1)
rm2 <- ode(i.state,1:Time, rmcr, p.rm2)
rm3 <- ode(i.state,1:Time, rmcr, p.rm3)
# Fig 2.1b and 2.6a, respectively
plot(rm1, main = c("Resource", "Consumer"), ylab = c("Density","Density"))
# Fig 2.6b (only Consumer shown in book)
plot(rm2, main = c("Resource", "Consumer"), ylab = c("Density","Density"))
# Fig 2.6c (only Consumer shown in book)
plot(rm3, main = c("Resource", "Consumer"), ylab = c("Density","Density"))

## Figure 2.7a.  All credit goes to this website: http://www.r-bloggers.com/r-tools-for-dynamical-systems-bifurcation-plot-in-r%C2%A0for%C2%A0system%C2%A0of%C2%A0odes/
param.name <- "a" # choose parameter to perturb
param.seq <- seq(1,2.5,length=50) # choose range of parameters
p.rm4 <- c(r = r, e = e, a = 1.0, K = K, Ro = Ro, m=m) # set starting parameters.

param.index <- which(param.name == names(p.rm4)) # tells the loop which parameter in "p.rm4" to grab for manipulation.
out <- list()
for (i in 1:length(param.seq))
  out[[i]] <- matrix(0, n, length(i.state))

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
range.lim <- apply(do.call("rbind", range.lim), 2, range)
plot.variable <- "C" # choose which variable to show
plot(0, 0, pch = "", xlab = param.name, ylab = plot.variable,
     xlim = range(param.seq), ylim = range.lim[,plot.variable])
for (i in 1:length(param.seq)) {
  points(rep(param.seq[i], Time), out[[i]][,plot.variable])
}

### Figure 2.7b. Tips taken from Enemy-Victim Interactions chapter of a Primer of Ecology in R

#### Currently not getting the dip in eigen value like I expected
dR.dt <- expression(r * R * (1 - R / K) - a * C * R / (R + Ro))
dC.dt <- expression(e * a * C * R / (R + Ro) - m * C)
RMjac1 <- list(D(dR.dt,"R"),D(dC.dt,"R"),D(dR.dt,"C"),D(dC.dt,"C"))

# confused about what I'm doing here...
R <- m * Ro / (e * a - m) # expression is the consumer isocline... but why?
C <- eval(Riso)

RM.jac2 <- matrix(sapply(RMjac1, function(pd) eval(pd)), nrow = 2)
eigen(RM.jac2)[["values"]]

rmcr.jacList <- lapply(1:length(param.seq), function(i) {
  a <- param.seq[i] # changes the parameter value for manipulation
  C <- eval(Riso)
  matrix(sapply(RMjac1, function(pd) eval(pd)), nrow=2)
})

L1 <- sapply(rmcr.jacList, function(J) max(Re(eigen(J)[["values"]])))


plot(param.seq,L1,type ="l",xlab = "a")
