######################## Functional Responses #####################################
a <- 0.1 # attack rate of predator (all functional response types)
w <- 0.1 # maximum attack rate of predator in Type 2 or 3 functional response
D <- w/a # half saturation constant. Only applicable to Type 2 or 3 functional response. Don't intuitively understand this.

par(mfrow = c(1,2)) # plot 2 graphs side-by-side

# first plot: functional responses
curve(a * x, 0, 2, xlab = "Prey Density", ylab = "Prey Killed per Predator") # Type 2
curve(w * x / (D + x), 0, 2, add = TRUE, lty = 2) # Type 2
curve(w * x^2/ (D^2 + x^2), 0, 2, add = TRUE, lty = 3) # Type 3

# second plot: functional response per prey
curve(a * x/x, 0, 2, xlab = "Prey Density", ylab = "Prey Killed per Predator per Prey", ylim = c(0,a))
curve(w * x / (D + x)/x, 0, 2, add = TRUE, lty = 2)
curve(w * x^2/ (D^2 + x^2)/x, 0, 2, add = TRUE, lty = 3)
legend("bottomright", c("Type 1", "Type 2", "Type 3"), lty = 1:3, bty = "n", cex = 0.8)

######### False target model for functional response.  This version was modified from Abrahamson & Weis 1997 pg. 305.  I believe their version is incorrect.  I tried using it but it does not replicate their graphs.  This current formulation does. ######
a <- 1
H <- 0.2 # handling time of vulnerable prey
D <- 0 # discrimination time between vulnerable and invulnerable prey
v <- 1.0 # fraction of vulnerable prey
v1 <- 0.5
v2 <- 0.25


par(mfrow=c(1,3)) # setup graphs for easy comparison

# Perfect discrimination, but decreasing the fraction of true targets
D <- 0
curve((a*v*x)/(1 + a*H*v*x + a*D*(1 - v)*x), from=0, to=50, ylim=c(0,5), ylab="Number encountered per parasitoid", xlab="Gall density", main="D = 0")
curve((a*v1*x)/(1 + a*H*v1*x + a*D*(1 - v1)*x), from=0, to=50, col=2, ylim=c(0,10), add=T)
curve((a*v2*x)/(1 + a*H*v2*x + a*D*(1 - v2)*x), from=0, to=50, col=3, add=T)
legend("bottomright", legend=c("v = 1.0", "v = 0.5", "v = 0.25"), col=1:3, lty=1, bty="n")

# weak discrimination, but decreasing the fraction of true targets
D <- H*0.5
curve((a*v*x)/(1 + a*H*v*x + a*D*(1 - v)*x), from=0, to=50, ylim=c(0,5), ylab="Number encountered per parasitoid", xlab="Gall density", main="D = H*0.5")
curve((a*v1*x)/(1 + a*H*v1*x + a*D*(1 - v1)*x), from=0, to=50, col=2, ylim=c(0,10), add=T)
curve((a*v2*x)/(1 + a*H*v2*x + a*D*(1 - v2)*x), from=0, to=50, col=3, add=T)
legend("bottomright", legend=c("v = 1.0", "v = 0.5", "v = 0.25"), col=1:3, lty=1, bty="n")

# anti-discrimination, but decreasing the fraction of true targets
D <- H*2
curve((a*v*x)/(1 + a*H*v*x + a*D*(1 - v)*x), from=0, to=50, ylim=c(0,5), ylab="Number encountered per parasitoid", xlab="Gall density", main="D = H*2")
curve((a*v1*x)/(1 + a*H*v1*x + a*D*(1 - v1)*x), from=0, to=50, col=2, ylim=c(0,10), add=T)
curve((a*v2*x)/(1 + a*H*v2*x + a*D*(1 - v2)*x), from=0, to=50, col=3, add=T)
legend("bottomright", legend=c("v = 1.0", "v = 0.5", "v = 0.25"), col=1:3, lty=1, bty="n")