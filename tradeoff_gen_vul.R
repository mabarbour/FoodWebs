
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(cowplot)
library(cheddar)

## PLOT RELATIONSHIP BETWEEN VULNERABILITY AND GENERALITY FOR EACH FOOD WEB ----

# TL84
data("TL84")
gen.vul.TL84 <- NPS(TL84, c("InDegree","OutDegree","PreyAveragedTrophicLevel","Biomass")) 
gen.vul.TL84$FoodWeb_ID <- "TL84"
ggplot(gen.vul.TL84, aes(x = OutDegree, y = InDegree, color = PreyAveragedTrophicLevel)) + geom_point(aes(size = Biomass), shape = 1) + geom_smooth(method = "glm", method.args = list(family = "poisson"))

# TL86
data("TL86")
gen.vul.TL86 <- NPS(TL86, c("InDegree","OutDegree","PreyAveragedTrophicLevel","Biomass"))
gen.vul.TL86$FoodWeb_ID <- "TL86"
ggplot(gen.vul.TL86, aes(x = OutDegree, y = InDegree, color = PreyAveragedTrophicLevel)) + geom_point(aes(size = Biomass), shape = 1) + geom_smooth(method = "glm", method.args = list(family = "poisson"))

# SkipwithPond
data("SkipwithPond")
gen.vul.SkipwithPond <- NPS(SkipwithPond, c("InDegree","OutDegree","PreyAveragedTrophicLevel","Biomass"))
gen.vul.SkipwithPond$FoodWeb_ID <- "SkipwithPond"
ggplot(gen.vul.SkipwithPond, aes(x = OutDegree, y = InDegree, color = PreyAveragedTrophicLevel)) + geom_point(shape = 1) + geom_smooth(method = "glm", method.args = list(family = "poisson"))

# Benguela
data("Benguela")
gen.vul.Benguela <- NPS(Benguela, c("InDegree","OutDegree","PreyAveragedTrophicLevel","Biomass"))
gen.vul.Benguela$FoodWeb_ID <- "Benguela"
ggplot(gen.vul.Benguela, aes(x = OutDegree, y = InDegree, color = PreyAveragedTrophicLevel)) + geom_point(shape = 1) + geom_smooth(method = "glm", method.args = list(family = "poisson"))

# BroadstoneStream
data("BroadstoneStream")
gen.vul.BroadstoneStream <- NPS(BroadstoneStream, c("InDegree","OutDegree","PreyAveragedTrophicLevel","Biomass"))
gen.vul.BroadstoneStream$FoodWeb_ID <- "BroadstoneStream"
ggplot(gen.vul.BroadstoneStream, aes(x = OutDegree, y = InDegree, color = PreyAveragedTrophicLevel)) + geom_point(aes(size = Biomass), shape = 1) + geom_smooth(method = "glm", method.args = list(family = "poisson"))

# YthanEstuary
data("YthanEstuary")
gen.vul.YthanEstuary <- NPS(YthanEstuary, c("InDegree","OutDegree","PreyAveragedTrophicLevel","Biomass"))
gen.vul.YthanEstuary$FoodWeb_ID <- "YthanEstuary"
ggplot(gen.vul.YthanEstuary, aes(x = OutDegree, y = InDegree, color = PreyAveragedTrophicLevel)) + geom_point(aes(size = Biomass), shape = 1) + geom_smooth(method = "glm", method.args = list(family = "poisson"))


## ALL FOOD WEBS TOGETHER ----
foodweb.data <- bind_rows(gen.vul.TL84, gen.vul.TL86, gen.vul.Benguela, gen.vul.BroadstoneStream, gen.vul.SkipwithPond, gen.vul.YthanEstuary)

ggplot(foodweb.data, aes(x = OutDegree, y = InDegree, color = PreyAveragedTrophicLevel, group = FoodWeb_ID)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm", method.args = list(family = "poisson")) +
  facet_wrap(~FoodWeb_ID) 


## NICHE MODEL RELATIONSHIPS ----

## from Ted Hart: https://gist.github.com/emhart/1503428

# from Williamns and Martinez: http://www.nature.com/nature/journal/v404/n6774/abs/404180a0.html
#  I have a brief undocumented function below to generate the web, it 
# requires S (# of species) and C (connectivity)

niche.model <- function(S,C){
  new.mat <- matrix(0,nrow=S,ncol=S)
  ci <- vector()
  niche <- runif(S,0,1)
  r <- rbeta(S,1,((1/(2*C))-1)) * niche
  for(i in 1:S){ci[i]<-runif(1,r[i]/2,niche[i])}
  
  #now set the smallest species niche value to have an n of 0
  r[which(niche==min(niche))] <- .00000001
  for(i in 1:S){
    for(j in 1:S){
      if(niche[j] > (ci[i]-(.5*r[i])) && niche[j]< (ci[i]+.5*r[i])){new.mat[j,i]<-1}
    }
  }
  new.mat <- new.mat[,order(apply(new.mat,2,sum))]
  return(new.mat)
}

# create parameters for simulating the food webs
connect.trials <- c(0.01, 0.05, 0.1, 0.2)
species <- 50

niche.webs <- list()
for(i in 1:length(connect.trials)){
  
  niche.webs.sameC <- list()
  for(j in 1:50){
    tmp.niche.web <- niche.model(S = species, C = connect.trials[i])
    generality <- colSums(tmp.niche.web)
    vulnerability <- rowSums(tmp.niche.web)

    niche.webs.sameC[[j]] <- data.frame(Trial = j, C = connect.trials[i], generality, vulnerability) %>% filter(generality > 0 | vulnerability > 0)
  }
  niche.webs[[i]] <- plyr::ldply(niche.webs.sameC)
}
niche.webs.df <- plyr::ldply(niche.webs) %>% mutate(Trial = as.factor(Trial))

ggplot(niche.webs.df, aes(x = vulnerability, y = generality)) + stat_smooth(aes(group = Trial), geom = "line", method = "glm", method.args = list(family = "poisson"), se = FALSE, alpha = 0.5) +
  facet_wrap(~C)

