source("utils.R")
source("fit_utils.R")


## get bug-fixed version of glmmTMB
remotes::install_github("glmmTMB/glmmTMB/glmmTMB")

library(remotes)
# oo <- options(repos = "https://cran.r-project.org/")
# utils::install.packages("Matrix")
# utils::install.packages("lme4")
# utils::install.packages("TMB")
# utils::install.packages("glmmTMB")
# options(oo)
#install.packages("lme4", type = "source")
library(lme4)
library(broom.mixed)
library(tidyverse)
# load_all_pkgs()
library(glmmTMB)
# library(Matrix)
library(gamm4)
ecoreg <- readRDS("ecoreg.rds")
system.time(f1 <- fit_all(single_fit=c(2,2,NA))) ## 0.5 seconds
system.time(f2 <- fit_all(single_fit=c(2,2,NA), platform = "glmmTMB",
                          add_sos = FALSE)) ## 6.5 seconds
system.time(f3 <- fit_all(single_fit=c(2,2,NA), platform = "glmmTMB")) ## 134 seconds  8.132 seconds
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

lmer_fit_REML=system.time(lmer(ff0,data = ecoreg,REML=TRUE)) #fitting the model using REML,  0.158 seconds
lmer_fit_REML0=system.time(lmer(ff0,data = ecoreg,REML=FALSE)) #fitting the model without using REML,  0.162 seconds

GLMTmb_fit_REML= system.time(glmmTMB(ff0,data = ecoreg, REML = TRUE) )#fitting the model using REML,  0.606
GLMTmb_fit_REML0= system.time(glmmTMB(ff0,data = ecoreg, REML = FALSE)) #fitting the model using REML,  0.518


gamm4_fit=gamm4(ff,data=ecoreg,REML = FALSE)




######fitting the glmm TMB in modular form ############


ff <- mbirds_log ~ (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc)^2 + 
  (1 + NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc || biome) +
  (1 + NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc || flor_realms) +
  area_km2_log_sc + s(y, x, bs = "sos")
ff0 <- update(ff, . ~ . - s(y, x, bs = "sos")) #the model with the smoothing splines

m0=system.time(glmmTMB(ff0, data = ecoreg)) #time taken 0.529 seconds

time_modular <- system.time({
  m1= glmmTMB(ff0,data=ecoreg,doFit = FALSE)

#names(m1)
m2= fitTMB(m1,doOptim = FALSE) #initiate the automatic differentiation 

#names(m2)
m2 <- with(m2$env,
           TMB::MakeADFun(data,
                          parameters,
                          map = map,
                          random = random,
                          silent = silent,
                          DLL = "glmmTMB"))

m3 <- with(m2, nlminb(par, objective = fn, gr = gr))
#names(m3)

m4 <- finalizeTMB(m1, m2, m3)
m4$call$doFit <- NULL ## adjust 'call' element to match
}) #time taken 0.331 seconds

all.equal(m0, m4)



##############################lets do without the smoothers #####################


ff <- mbirds_log ~ (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc)^2 + 
  (1 + NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc || biome) +
  (1 + NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc || flor_realms) +
  area_km2_log_sc + s(y, x, bs = "sos")
ff0 <- update(ff, . ~ . - s(y, x, bs = "sos")) #the model with the smoothing splines

m0=system.time(glmmTMB(ff, data = ecoreg)) #time taken 0.529 seconds

time_modular_nosmooth <- system.time({
  m1= glmmTMB(ff,data=ecoreg,doFit = FALSE)
  
  #names(m1)
  m2= fitTMB(m1,doOptim = FALSE) #initiate the automatic differentiation 
  
  #names(m2)
  m2 <- with(m2$env,
             TMB::MakeADFun(data,
                            parameters,
                            map = map,
                            random = random,
                            silent = silent,
                            DLL = "glmmTMB"))
  
  m3 <- with(m2, nlminb(par, objective = fn, gr = gr))
  #names(m3)
  
  m4 <- finalizeTMB(m1, m2, m3)
  m4$call$doFit <- NULL ## adjust 'call' element to match
}) #Warning message:
#In finalizeTMB(m1, m2, m3) :
  #Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 7.793 seconds
all.equal(m0, m4)




