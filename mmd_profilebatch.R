source("mmd_utils.R")
load("allfits.RData")
library(lme4)

proffun <- function(fit) {
    best_model <- fit[[get_best_name(fit)]]
    pp <- profile(best_model)
    return(pp)
}
profList <- parallel::mclapply(allfits, proffun,
                              mc.cores = 2  ## change as available
                              )
names(profList) <- names(allfits)
save("profList",file="allprofs.RData")

## profiles ugly for 1 (plants), 4 (mammals);
## plotting fails for amphibians (but confints work, falling back to lin)

