dynamic_matplot <- function(param.vector, # vector of parameters for the dynamical model
                            init.state, # initial state variables for the model
                            sim.length, # length of simulation
                            model, # dynamical model
                            ylim, # limits of state variables for plotting
                            #cex,
                            ...){
  
  # Run the experiment
  df <- ode(init.state, 1:sim.length, model, param.vector)
  
  # plot the results. 
  matplot(df[ ,"time"], df[ ,names(init.state)],
          type = "l", ylim=ylim, ...)
  #legend("top", paste(names(init.state)), lty=1:length(init.state), 
   #        col=1:length(init.state), bty="n", ncol = length(init.state),
    #       cex = cex)
}