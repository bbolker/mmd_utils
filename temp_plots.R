ecoreg <- readRDS("ecoreg.rds")
library(tidyverse)

e2 <- (ecoreg
    %>% group_by(biome)
    %>% summarise(across(c(Feat_mean,NPP_mean,mbirds,mmamm,mamph), mean))
    %>% mutate(across(biome, ~reorder(factor(.), NPP_mean)))
)

barplot(e2$Feat_mean,names=e2$biome)
library(ggplot2)
ggplot(e2, aes(NPP_mean, biome)) + geom_bar(stat="identity")
ggplot(e2, aes(NPP_mean, mmamm, label=biome)) + geom_point() + geom_label(size=16)

## log(y) = a + b*log(x)
##  ->   y = exp(a)*x^b
##  ->  b (regression coefficient of y vs x) is a power-law
## e.g. if b ~ 0.5 we could say "mammal diversity increases with
##   the square root of NPP"
