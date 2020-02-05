source("utils.R")
allfits_lme4 <- readRDS("allfits_lme4.rds")
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
saveRDS(profList,file="allprofs.rds")

## profiles ugly for 1 (plants), 4 (mammals);
## plotting fails for amphibians (but confints work, falling back to lin)

