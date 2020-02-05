args <- commandArgs(trailingOnly=TRUE)
print(args)
## 'lme4', 'gamm4', 'brms'
platform <- if (length(args)<1) "lme4" else args[1]
## restricted fit?
restr <- if (length(args)<2) FALSE else as.logical(args[2])
## specify number of parallel cores
cores <- if (length(args)<3) {
             if (platform!="brms") 2 else 1
         }

ecoreg <- readRDS("ecoreg.rds")
source("utils.R")
require(platform, character.only=TRUE)
library(parallel)

## skip plants for now ...
respvars <- c("mamph","mbirds","mmamm") ## "plants"
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

## ff(logrespvars[1], platform=platform)
allfits <- parallel::mclapply(logrespvars, ff,
                              platform = platform,
                              verbose = TRUE,
                              mc.cores = cores  ## change as available
                              )
names(allfits) <- logrespvars

if (!restr) {
    fn <- sprintf("allfits_%s.rds", platform)
} else {
    fn <- sprintf("allfits_restr_%s.rds", platform)
}
saveRDS(allfits,file=fn)
