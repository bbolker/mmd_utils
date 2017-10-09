L <- load("ecoreg2.RData")
source("mmd_utils.R")
library(lme4)

## pick out all response variables
logrespvars <- paste0(respvars,"_log")

allfits <- list()
for (r in logrespvars) {
    cat(r,"\n")
    allfits[[r]] <- fit_all(response=r)
}
save("allfits",file="allfits.RData")
