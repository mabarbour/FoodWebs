# May Simulation
# based off code provided here: http://ecovirtual.ib.usp.br/lib/exe/fetch.php?media=ecovirt:roteiro:math:eq_funcoes.r

## funcao para a simulacao do May
may <- function(S,C,f,nsim=100){
  m <- diag(rep(-1,S))
  ind <- which(m==0,arr.ind=TRUE)
  n <- round((C*S^2)-S,0)
  n.estavel <- 0
  for(i in 1:nsim){
    vals <- rnorm(n,mean=0,sd=f)
    ri <- sample(1:nrow(ind),size=n)
    m2 <- m
    m2[ind[ri,]] <- vals
    autov <- eigen(m2,only.values=TRUE)$values
    n.estavel <- n.estavel+all(Re(autov)<0)
  }
  data.frame(S=S,C=C,f=f,p.estab=n.estavel/nsim)
}