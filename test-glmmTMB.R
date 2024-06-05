source("utils.R")
source("fit_utils.R")

## get bug-fixed version of glmmTMB
remotes::install_github("glmmTMB/glmmTMB/glmmTMB", ref = "1044")
library(broom.mixed)
library(tidyverse)
load_all_pkgs()
library(glmmTMB)

ecoreg <- readRDS("ecoreg.rds")
system.time(f1 <- fit_all(single_fit=c(2,2,NA))) ## 0.5 seconds
system.time(f2 <- fit_all(single_fit=c(2,2,NA), platform = "glmmTMB",
                          add_sos = FALSE)) ## 6.5 seconds
system.time(f3 <- fit_all(single_fit=c(2,2,NA), platform = "glmmTMB")) ## 134 seconds
system.time(f4 <- fit_all(single_fit=c(2,2,NA), platform = "gamm4")) ## 6 seconds

fit_list <- list(lme4=f1,
                 glmmTMB_nospline = f2,
                 glmmTMB = f3,
                 gamm4 = f4)
## compare fits
(purrr::map_dfr(fit_list, glance, .id = "model")
    |> dplyr::select(model, sigma, logLik, AIC)
)

## ! is glmmTMB really that much better? or difference in loglik calc?
## coeffs don't look all that different ...
## TO DO: break out argList etc. to understand why lme4/gamm4 are so much
##  faster here?

ff <- mbirds_log ~ (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc)^2 + 
    (1 + NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc || biome) +
    (1 + NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc || flor_realms) +
    area_km2_log_sc + s(y, x, bs = "sos")
ff0 <- update(ff, . ~ . - s(y, x, bs = "sos"))

glmmTMB(ff0, data = ecoreg)

