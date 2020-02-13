## testing backtrans, plotting, etc.
ecoreg <- readRDS("ecoreg.rds")
best_models <- readRDS("bestmodels_gamm4.rds")
allfits_brms <- readRDS("allfits_brms.rds")

library(gamm4)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
source("utils.R")

bb <- best_models[["mmamm_log"]]
head(residuals(bb$mer))
bbr <- allfits_brms[["mmamm_log"]]


## plain
no_legend <- theme(legend.position="none")
plotfun(bb,backtrans=TRUE,log="xy",xvar="Feat_log_sc")
plotfun(bb,backtrans=TRUE,xvar="Feat_log_sc")
plotfun(bb,xvar="Feat_log_sc",backtrans=TRUE,log="xy")+ no_legend
## brms prediction
plotfun(bbr,backtrans=TRUE,log="xy",ylim=log(c(8,250)))+ no_legend
## all looks reasonable

## test remef
ecoreg$rem1 <- drop(remef_allran(best_models[["mmamm_log"]],
                                 data=ecoreg,
                                 ## set_other="median",
                                 set_other="zero",
                                 keep_intercept=TRUE,
                                 ##fixed_keep=c("(Intercept)","NPP_log")
                                 fixed_keep=c("NPP_log")
                                 ))

## raw/simple
rem2 <- residuals(bb$mer)+cc[1]+cc[2]*ecoreg$NPP_log
plot(rem1~NPP_log,data=ecoreg)
plot(rem2~NPP_log,data=ecoreg)
cc <- coef(bb$gam)[c("(Intercept)","NPP_log")]
abline(a=cc[1],b=cc[2])
## looks OK ?

## remef_allran
plotfun(bb,respvar="rem1",ylim=c(-1,0.8))+
    theme(legend.position="none")+
    geom_smooth(method="lm",aes(group=1,y=rem1),
                formula=y~x-1)
## why have regression/predicted lines disappeared now?
## x,y set to zero??

## something about intercepts: centering/ 0 values?

plotfun(bb,respvar="rem1",ylim=c(-1,0.8),auxvar=NULL,
        adjust_othervar="mean")+
    theme(legend.position="none")+
    geom_smooth(method="lm",aes(group=1,y=rem1))

plotfun(bb,respvar="rem1",ylim=c(-1,0.8))+
    theme(legend.position="none")+
    geom_smooth(method="lm",aes(group=1,y=rem1))
## doesn't match exactly, but close

plotfun(bb)+theme(legend.position="none")

plotfun(bb,auxvar=NULL,grpvar="biome")+theme(legend.position="none")
    geom_smooth(method="lm",aes(group=1))

###

dd <- data.frame(x1=rnorm(100))
dd$x2 <- -0.5*x1+rnorm(100)
dd$y <-


###
rr <- remef_allran(best_models[["mmamm_log"]],data=ecoreg,
                   fixed_keep=NA,return_components=TRUE)
pred2 <- with(rr,pp_ran+pp_fixed)
pred1 <- predict(best_models[["mmamm_log"]]) ## doesn't include smooth term
par(mfrow=c(1,2))
plot(ecoreg$NPP_log,pred1-pred2)
plot(ecoreg$NPP_log,rr$rem)
dd <- data.frame(ecoreg,rr)

       
