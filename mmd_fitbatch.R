L <- load("ecoreg.RData")
source("mmd_utils.R")
library(lme4)
library(gamm4)
library(parallel)

logrespvars <- paste0(c("plants","mamph","mbirds","mmamm"),"_log")
## pp <- ecoreg[c("plants_log","NPP_mean","NPP_cv_inter",
##                     "Feat_mean","Feat_cv_inter",
##                     "biome","flor_realms","area_km2","x","y")]
## pp[!complete.cases(pp),]  ## NA values in predictor
## f1 <- fit_all(logrespvars[1],use_gamm4=TRUE,verbose=TRUE)

## this was run earlier with lme4 (use_gamm4=FALSE)
## and stored in allfits_lmer.RData

allfits <- parallel::mclapply(logrespvars, fit_all,
                              platform = "gamm4",
                              verbose = TRUE,
                              mc.cores = 2  ## change as available
                              )
names(allfits) <- logrespvars
sapply(allfits, function(x) is(x,"try-error"))

save("allfits",file="allfits.RData")
