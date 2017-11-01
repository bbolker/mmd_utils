L <- load("ecoreg2.RData")
source("mmd_utils.R")
library(lme4)
library(parallel)

## pick out all response variables
logrespvars <- paste0(respvars,"_log")
tmpf <- function(r) fit_all(response=r,
                            single_fit=c(2,2),
                            rterms=c("biome","flor_realms"))
allfits_restr <- parallel::mclapply(logrespvars,
                              tmpf,
                              mc.cores = 2  ## change as available
                              )
names(allfits_restr) <- logrespvars

save("allfits_restr",file="allfits_restr.RData")
