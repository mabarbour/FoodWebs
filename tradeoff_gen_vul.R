
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
