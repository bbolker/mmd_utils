## using brms; fit only most complex model
L <- load("ecoreg.RData")
source("mmd_utils.R")
library(brms)
library(parallel)
options(mc.cores=3)

## pick out all response variables
respvars <- c("mamph","plants","mmamm","mbirds")
logrespvars <- paste0(respvars,"_log")
tmpf <- function(r) fit_all(response=r,
                            single_fit=c(3,3,3),
                            platform="brms")
## don't fit respvars in parallel because brm is doing parallel
##  chains (maybe modify if we move to SHARCnet)
allfits_brms <- lapply(logrespvars,tmpf)
names(allfits_brms) <- logrespvars

save("allfits_brms",file="allfits_brms.RData")
