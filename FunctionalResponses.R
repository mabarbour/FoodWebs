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
