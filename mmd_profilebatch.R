source("mmd_utils.R")
load("allfits.RData")
library(lme4)
profList <- list()
for (i in seq_along(allfits)) {
    nm <- names(allfits)[i]
    cat(nm,"\n")
    best_model <- allfits[[i]][[get_best_name(allfits[[i]])]]
    profList[[i]] <- profile(best_model)
    names(profList)[i] <- nm
    save("profList",file="allprofs.RData")
}

## profiles ugly for 1 (plants), 4 (mammals);
## plotting fails for amphibians (but confints work, falling back to lin)

