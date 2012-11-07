##### Notes from Enemy-Victim Interactions Chapter of a Primer of Ecology with R

######################### Lotka-Voltera Model ###############################
# Instability is a fundamental feature.  When prey reproduce and are limited only by the predator, and the predators are limited by only the abundance of prey, these interactions are not stable.

predpreyLV <- function(t,y,params) {
  H <- y[1] # H = herbivores
  P <- y[1] # P = predators
  with(as.list(params), {
    dH.dt <- b * H - a * P * H # a = attack rate of predator; b = intrinsic growth rate
    dP.dt <- e * a * P * H - s * P 
    
    # note how the loss term of the herbivore (-a * P * H) directly translates to predator growth multiplied by e = conversion efficiency.
    # s = per captia mortality rate
    # note that predators kill herbivores at a rate = a * H.  This a a linear functional response (Type 1), indicating the predators kill a fixed proportion of prey.  In other words, the rate at which predators kill prey is INDEPENDENT of herbivore density.
    # notice how the numerical response of predator incorporates the functional response, multiplied by its conversion efficiency (i.e. number of predators per number of herbivores killed) minus the number that die from other causes (s * P).
  })
}

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
