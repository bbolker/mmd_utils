args <- commandArgs(trailingOnly=TRUE)
print(args)
## 'lme4', 'gamm4', 'brms'
platform <- if (length(args)<1) "lme4" else args[1]
if (!platform %in% c("lme4","gamm4","brms")) stop("unknown platform",platform)
## restricted fit?
restr <- if (length(args)<2) FALSE else as.logical(args[2])
## specify number of parallel cores
## don't use multicores for brms (by default) because it parallelizes chains anyway
cores <- if (length(args)<3) {
             if (platform!="brms") 2 else 1
         } else as.numeric(args[3])
## include plants?
include_plants <- if (length(args)<4) TRUE else as.logical(args[4])

ecoreg <- readRDS("ecoreg.rds")
source("utils.R")
require(platform, character.only=TRUE)
library(parallel)

respvars <- c("mamph","mbirds","mmamm")
if (include_plants) {
    respvars <- c(respvars,"plants")
}
logrespvars <- paste0(respvars,"_log")

ff <- if (restr) {
          function(r,...) fit_all(response=r,
                                  single_fit=c(2,2),
                                  rterms=c("biome","flor_realms"),
                                  ...)
      } else if (platform=="brms") {
          function(r,...) fit_all(response=r,
                                  single_fit=c(3,3,3),
                                  ...)

      } else {
          fit_all
      }

## ff(logrespvars[1], platform=platform, verbose=TRUE)
## FIXME:: test/warning on Windows?
allfits <- parallel::mclapply(logrespvars, ff,
                              platform = platform,
                              verbose = TRUE,
                              mc.cores = cores  ## change as available
                              )
## FIXME: do we need this or does mclapply assign names?
names(allfits) <- logrespvars

if (!restr) {
    fn <- sprintf("allfits_%s.rds", platform)
} else {
    fn <- sprintf("allfits_restr_%s.rds", platform)
}
saveRDS(allfits,file=fn)
