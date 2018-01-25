L <- load("ecoreg.RData")
source("mmd_utils.R")
library(lme4)
library(parallel)

logrespvars <- paste0(c("plants","mamph","mbirds","mmamm"),"_log")

allfits <- parallel::mclapply(logrespvars, fit_all,
                              ## use_gamm4 = TRUE,
                              verbose = TRUE,
                              mc.cores = 2  ## change as available
                              )
names(allfits) <- logrespvars
sapply(allfits, function(x) is(x,"try-error"))

save("allfits",file="allfits.RData")
