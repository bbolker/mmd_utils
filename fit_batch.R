args <- commandArgs(trailingOnly=TRUE)
## arguments: platform, restricted fit?, num cores, include plants?, exclude fir?
## args <- c("gamm4",FALSE,1,TRUE,TRUE)
print(args)
## 'lme4', 'gamm4', 'brms'
platform <- if (length(args)<1) "lme4" else args[1]
if (!platform %in% c("lme4","gamm4","brms","glmmTMB")) stop("unknown platform",platform)
## restricted fit?
restr <- if (length(args)<2) FALSE else as.logical(args[2])
## specify number of parallel cores
## don't use multicores for brms (by default) because it parallelizes chains anyway
cores <- if (length(args)<3) {
             if (platform!="brms") 2 else 1
         } else as.numeric(args[3])
## include plants?
include_plants <- if (length(args)<4) TRUE else as.logical(args[4])
exclude_fire <- if (length(args)<5) FALSE else as.logical(args[5])

cat(sprintf("platform=%s restr=%d cores=%d include_plants=%d exclude_fire=%d",
            platform,restr,cores,include_plants, exclude_fire),"\n")

ecoreg <- readRDS("ecoreg.rds")
source("fit_utils.R")
source("gamm4_utils.R")
require(platform, character.only=TRUE)
library(parallel)

respvars <- c("mamph","mbirds","mmamm")

if (include_plants) {
    respvars <- c(respvars,"plants")
}
logrespvars <- paste0(respvars,"_log")

if (!exclude_fire) {
    pars <- c("NPP_log_sc","Feat_log_sc","NPP_cv_sc","Feat_cv_sc")
} else {
    pars <- c("NPP_log_sc","NPP_cv_sc")
}

ff <- if (restr) {
          function(r,...) fit_all(response=r,
                                  single_fit=c(2,2),
                                  rterms=c("biome","flor_realms"),
                                  pars=pars,
                                  ...)
      } else if (platform=="brms") {
          function(r,...) fit_all(response=r,
                                  single_fit=c(3,3,3),
                                  pars=pars,
                                  ...)

      } else {
          function(r,...) fit_all(response=r,pars=pars, ...)
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

str1 <- "allfits"
str2 <- if (restr) "restr" else ""
str3 <- if (exclude_fire) "nofire" else ""
str4 <- platform
ss <- c(str1,str2,str3,str4)
ss <- ss[nchar(ss)>0]
fn <- paste0(paste(ss,collapse="_"),".rds")

saveRDS(allfits,file=fn)
