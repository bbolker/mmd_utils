## using brms; fit only most complex model
L <- load("ecoreg.rds")
source("utils.R")
library(brms)
library(parallel)
options(mc.cores=3)  ## use only 1 core for debugging; 3+ for prod

## pick out all response variables
respvars <- c("mamph","plants","mmamm","mbirds")
logrespvars <- paste0(respvars,"_log")
tmpf <- function(r) {
    res <- fit_all(response=r,
            single_fit=c(3,3,3),
            platform="brms")
    return(res)
}
## don't fit respvars in parallel because brm is doing parallel
##  chains (maybe modify if we move to SHARCnet)
## for loop rather than lapply() so we can checkpoint

allfits_brms <- setNames(vector("list",length=length(logrespvars)),logrespvars)
for (i in seq_along(allfits_brms)) {
    allfits_brms[[i]] <- tmpf(logrespvars[i])
    ## checkpoint
    save("allfits_brms",file="allfits_brms.RData")
}
save("allfits_brms",file="allfits_brms.RData")
