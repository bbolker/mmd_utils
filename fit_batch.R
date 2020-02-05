args <- commandArgs(trailingOnly=TRUE)
platform <- if (length(args)<1) "lme4" else args[1]
restr <- if (length(args)<2) FALSE else as.logical(args[2])

L <- load("ecoreg.RData")
source("utils.R")
library(lme4)
library(gamm4)
library(parallel)

## skip plants for now ...
logrespvars <- paste0(c("mamph","mbirds","mmamm"),"_log")

restr_fit <- function(r) fit_all(response=r,
                            single_fit=c(2,2),
                            platform=platform,
                            rterms=c("biome","flor_realms"))

ff <- if (restr) restr_fit else fit_all

allfits <- parallel::mclapply(logrespvars, ff,
                              platform = platform,
                              verbose = TRUE,
                              mc.cores = 2  ## change as available
                              )
names(allfits) <- logrespvars

if (!restr) {
    fn <- sprintf("allfits_%s.rds", platform)
} else {
    fn <- sprintf("allfits_restr_%s.rds", platform)
}
saveRDS(allfits,file=fn)
