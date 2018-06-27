L <- load("ecoreg.RData")
source("mmd_utils.R")
library(lme4)
library(gamm4)
library(parallel)

## pick out all response variables
respvars <- c("mamph","plants","mmamm","mbirds")
logrespvars <- paste0(respvars,"_log")
tmpf <- function(r) fit_all(response=r,
                            single_fit=c(2,2),
                            platform="gamm4",
                            rterms=c("biome","flor_realms"))
allfits_restr_gamm4 <- parallel::mclapply(logrespvars,
                              tmpf,
                              mc.cores = 2  ## change as available
                              )
names(allfits_restr_gamm4) <- logrespvars

save("allfits_restr_gamm4",file="allfits_restr_gamm4.RData")
