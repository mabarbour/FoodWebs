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

############### Rosenzweig-MacArthur consumer-resource model #######################
# assumes logistic resource growth and a type 2 functional response by the consumer (i.e. consumption rate saturates with resource density). This function is necessary for the ordinary differential equation solver, ode() 
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

############# Continuous logistic growth model ##################################
# note that this is a phenomological model that tends to explain isolated single species population dynamics well
clg <- function(t,y,p) {
  R <- y[1]
  with(as.list(p), {
    dR.dt <- r * R * (1 - R / K)
    return(list(c(dR.dt)))
  })
}

############# Ricker Model - discrete time population model ###########################
# these models can account for lagged population responses, which may be common to organisms that reproduce on seasonal cycles.
# help for writing this model came from Duke University's website for Bio 292: Population Ecology. Code was written by wfmorris

#Nt1=Nt*exp(B*(1-Nt/K)
