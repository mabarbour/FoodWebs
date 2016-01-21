## This script contains a variety of functions for consumer-resource models.
# all of these functions are in the necessary format "function(t,y,p)" for solving used the ode() function in the R package "deSolve"

#### Rosenweig-MacArthur type models
# All of these models make the following assumptions:
# (1) logistic growth for the resource(s)
# (2) type 2 functional response of consumer
# (3) density-independent mortality for consumer 

## Apparent competition model: 1 consumer, 2 resources ----
# Contains a general form for a type 2 functional response for a system with 2 resource types (Murdoch and Oaten 1975). Essentially, consumers switch to alternative prey when those prey have a higher relative density.
# assumes conversion efficiency is equal on both resources.
RM_apparent <- function(t,y,p) {
  # t,y,p is the necessary format for solving using the ode() function in R
  R1 <- y[1]
  R2 <- y[2]
  C1 <- y[3]
  with(as.list(p), {
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - 
              a1 * C1 * R1 / (1 + a1 * h1 * R1 + a2 * h2 * R2)
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - 
              a2 * C1 * R2 / (1 + a2 * h2 * R2 + a1 * h1 * R1)
    dC1.dt <- e1 * (a1 * C1 * R1 / (1 + a1 * h1 * R1 + a2 * h2 * R2) + 
                   a2 * C1 * R2 / (1 + a2 * h2 * R2 + a1 * h1 * R1)) - 
             dC1* C1 
    return(list(c(dR1.dt, dR2.dt, dC1.dt)))
  })
}

## 2 consumers, 2 resources ----
# Contains a general form for a type 2 functional response for a system with 2 resource types (Murdoch and Oaten 1975)
# Followed Abrams 1980: Consumer functional response and competition consumer-resource systems
# Note that only 2 attack rates and 2 handling times are modelled. This preserves the symmetry of the system. i.e. Consumer 1's handling time and attack rate for resource 1 is equal to Consumer 2's handling time and attack rate for resource 2. Similarly, Consumer 1's handling time and attack rate for resource 2 equals Consumer 2's handling time and attack rate for resource 1. 
# Also assumes that each consumer's conversion efficiency is equal for both resources.
RM_2C_2R <- function(t,y,p) {
  R1 <- y[1]
  R2 <- y[2]
  C1 <- y[3]
  C2 <- y[4]
  with(as.list(p), {
    # resource equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - 
      a1 * C1 * R1 / (1 + a1 * h1 * R1 + a2 * h2 * R2) -
      a2 * C2 * R1 / (1 + a1 * h1 * R2 + a2 * h2 * R1)
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - 
      a2 * C1 * R2 / (1 + a1 * h1 * R1 + a2 * h2 * R2) -
      a1 * C2 * R2 / (1 + a1 * h1 * R2 + a2 * h2 * R1)
    
    # consumer equations
    dC1.dt <- e1 * (a1 * C1 * R1 / (1 + a1 * h1 * R1 + a2 * h2 * R2) + 
                    a2 * C1 * R2 / (1 + a1 * h1 * R1 + a2 * h2 * R2)) - 
              dC1 * C1 
    dC2.dt <- e2 * (a2 * C2 * R1 / (1 + a1 * h1 * R2 + a2 * h2 * R1) +
                    a1 * C2 * R2 / (1 + a1 * h1 * R2 + a2 * h2 * R1)) -
              dC2 * C2
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}

## 2 consumer, 2 resource model ----
# This model defers from the previous one in that it contains 4 attack rates, 4 handling times, and 4 conversion efficiencies.
# Contains a general form for a type 2 functional response for a system with 2 resource types (Murdoch and Oaten 1975).
RM_2C_2R_full <- function(t,y,p) {

  R1 <- y[1]
  R2 <- y[2]
  C1 <- y[3]
  C2 <- y[4]
  
  # functional responses of each pairwise consumer-resource interactions.
  C1_R1_f_response <- a11 * C1 * R1 / (1 + a11 * h11 * R1 + a12 * h12 * R2)
  C2_R1_f_response <- a21 * C2 * R1 / (1 + a21 * h21 * R1 + a22 * h22 * R2)
  C1_R2_f_response <- a12 * C1 * R2 / (1 + a12 * h12 * R2 + a11 * h11 * R1)
  C2_R2_f_response <- a22 * C2 * R2 / (1 + a22 * h22 * R2 + a21 * h21 * R1)
  
  with(as.list(p), {
    # resource equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - 
      C1_R1_f_response - C2_R1_f_response
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - 
      C1_R2_f_response - C2_R2_f_response
  
    # consumer equations
    dC1.dt <- e11 * C1_R1_f_response +
              e12 * C1_R2_f_response - 
              dC1 * C1 
    dC2.dt <- e22 * C2_R2_f_response +
              e21 * C2_R1_f_response -
              dC2 * C2
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}

## Apparent competition model (1 consumer, 2 resources) that also has a spatially implicit multispecies functional response. Details in McCann et al. 2005 and McCann 2012. ----
RM_C_2R_space <- function(t,y,p) {

  R1 <- y[1]
  R2 <- y[2]
  C1 <- y[3]
  
  # preference functions
  W11 <- w * R1 / (w * R1 + (1 - w) * R2)
  W12 <- ((1 - w) * R2) / ((1 - w) * R2 + w * R1)
  
  # functional responses of each pairwise consumer-resource interactions. 
  C1_R1_f_response <- W11 * a11 * R1 / (1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
  C1_R2_f_response <- W12 * a12 * R2 / (1 + W12 * a12 * h12 * R2 + W11 * a11 * h11 * R1)
  
  with(as.list(p), {
    # resource equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - 
      C1_R1_f_response * C1
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - 
      C1_R2_f_response * C1
    
    # consumer equations
    dC1.dt <- C1 * (e11 * C1_R1_f_response + e12 * C1_R2_f_response - dC1)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt)))
  })
}

## 2 consumer, 2 resource model that also has a spatially implicit multispecies functional response. Details in McCann et al. 2005 and McCann 2012. ----
# this model permits 4 attack rates, 4 handling times, and 4 conversion efficiencies.
RM_2C_2R_space <- function(t,y,p) {

  R1 <- y[1]
  R2 <- y[2]
  C1 <- y[3]
  C2 <- y[4]
  
  ## preference functions. For right now, I'm going to set w1 = w2. In other words, there is only 1 preference parameter (w). This maintains symmetry in the consumers. I have coded out the version where w1 and w2 are separate in case I want to return to that model.
  #W11 <- w1 * R1 / (w1 * R1 + (1 - w1) * R2)
  #W12 <- ((1 - w1) * R2) / ((1 - w1) * R2 + w1 * R1)
  #W22 <- w2 * R2 / (w2 * R2 + (1 - w2) * R1)
  #W21 <- ((1 - w2) * R1) / ((1 - w2) * R1 + w2 * R2)
  
  W11 <- w * R1 / (w * R1 + (1 - w) * R2)
  W12 <- ((1 - w) * R2) / ((1 - w) * R2 + w * R1)
  W22 <- w * R2 / (w * R2 + (1 - w) * R1)
  W21 <- ((1 - w) * R1) / ((1 - w) * R1 + w * R2)
  
  # functional responses of each pairwise consumer-resource interactions. 
  C1_R1_f_response <- W11 * a11 * R1 / (1 + W11 * a11 * h11 * R1 + W12 * a12 * h12 * R2)
  C1_R2_f_response <- W12 * a12 * R2 / (1 + W12 * a12 * h12 * R2 + W11 * a11 * h11 * R1)
  
  C2_R2_f_response <- W22 * a22 * R2 / (1 + W22 * a22 * h22 * R2 + W21 * a21 * h21 * R1)
  C2_R1_f_response <- W21 * a21 * R1 / (1 + W21 * a21 * h21 * R1 + W22 * a22 * h22 * R2)
  
  with(as.list(p), {
    # resource equations
    dR1.dt <- r1 * R1 * (1 - R1 / K1) - 
      C1_R1_f_response * C1 - C2_R1_f_response * C2
    dR2.dt <- r2 * R2 * (1 - R2 / K2) - 
      C2_R2_f_response * C2 - C1_R2_f_response * C1
    
    # consumer equations
    dC1.dt <- C1 * (e11 * C1_R1_f_response + e12 * C1_R2_f_response - dC1)
    dC2.dt <- C2 * (e22 * C2_R2_f_response + e21 * C2_R1_f_response - dC2)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}

## Simplified 2 consumer, 2 resource model with spatially implicit multispecies functional response. Details in McCann et al. 2005 and McCann 2012. ----
# Model permits 2 different attack rates, 1 conversion efficiency, and 1 handling time. 2 different resource growth rates though, and 1 resource carrying capacity.
RM_2C_2R_space_simple <- function(t,y,p) {

  R1 <- y[1]
  R2 <- y[2]
  C1 <- y[3]
  C2 <- y[4]
  
  ## preference functions. For right now, these are all based off the parameter 'w', which maintains symmetry between consumers.
  W11 <- w * R1 / (w * R1 + (1 - w) * R2)
  W12 <- ((1 - w) * R2) / ((1 - w) * R2 + w * R1)
  
  W22 <- w * R2 / (w * R2 + (1 - w) * R1)
  W21 <- ((1 - w) * R1) / ((1 - w) * R1 + w * R2)

  # functional responses of each pairwise consumer-resource interactions. a_prime = attack rate on focal resource (C1-R1, C2-R2). a_alt = attack rate on alternative resource (C1-R2, C2-R1)
  C1_R1_f_response <- W11 * a_prime * R1 / (1 + W11 * a_prime * h * R1 + W12 * a_alt * h * R2)
  C1_R2_f_response <- W12 * a_alt * R2 / (1 + W12 * a_alt * h * R2 + W11 * a_prime * h * R1)
  
  C2_R2_f_response <- W22 * a_prime * R2 / (1 + W22 * a_prime * h * R2 + W21 * a_alt * h * R1)
  C2_R1_f_response <- W21 * a_alt * R1 / (1 + W21 * a_alt * h * R1 + W22 * a_prime * h * R2)
  
  with(as.list(p), {
    # resource equations
    dR1.dt <- r1 * R1 * (1 - R1 / K) - 
      C1_R1_f_response * C1 - C2_R1_f_response * C2
    dR2.dt <- r2 * R2 * (1 - R2 / K) - 
      C2_R2_f_response * C2 - C1_R2_f_response * C1
    
    # consumer equations
    dC1.dt <- C1 * (e * C1_R1_f_response + e * C1_R2_f_response - d)
    dC2.dt <- C2 * (e * C2_R2_f_response + e * C2_R1_f_response - d)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}

## Simplified Rosenzweig-MacArthur type 2 resource, 2 consumer model
# this model is just like RM_2C_2R_space_simple, except that now the consumers perceives both resource types as well-mixed (classic assumption of multispecies function response).
RM_2C_2R_NOspace_simple <- function(t,y,p) {
  R1 <- y[1]
  R2 <- y[2]
  C1 <- y[3]
  C2 <- y[4]
  
  ## no preference function. Setting them all equal to 1 returns the model to the classic multispecies functional response (both resources are well-mixed).
  W11 <- 1
  W12 <- 1
  W22 <- 1
  W21 <- 1
  
  # functional responses of each pairwise consumer-resource interactions. a_prime = attack rate on focal resource (C1-R1, C2-R2). a_alt = attack rate on alternative resource (C1-R2, C2-R1)
  C1_R1_f_response <- W11 * a_prime * R1 / (1 + W11 * a_prime * h * R1 + W12 * a_alt * h * R2)
  C1_R2_f_response <- W12 * a_alt * R2 / (1 + W12 * a_alt * h * R2 + W11 * a_prime * h * R1)
  
  C2_R2_f_response <- W22 * a_prime * R2 / (1 + W22 * a_prime * h * R2 + W21 * a_alt * h * R1)
  C2_R1_f_response <- W21 * a_alt * R1 / (1 + W21 * a_alt * h * R1 + W22 * a_prime * h * R2)
  
  with(as.list(p), {
    # resource equations
    dR1.dt <- r1 * R1 * (1 - R1 / K) - 
      C1_R1_f_response * C1 - C2_R1_f_response * C2
    dR2.dt <- r2 * R2 * (1 - R2 / K) - 
      C2_R2_f_response * C2 - C1_R2_f_response * C1
    
    # consumer equations
    dC1.dt <- C1 * (e * C1_R1_f_response + e * C1_R2_f_response - d)
    dC2.dt <- C2 * (e * C2_R2_f_response + e * C2_R1_f_response - d)
    
    return(list(c(dR1.dt, dR2.dt, dC1.dt, dC2.dt)))
  })
}
