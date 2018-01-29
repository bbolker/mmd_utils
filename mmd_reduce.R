## extract useful/necessary information
library(purrr)
library(dplyr)
library(broom)
library(broom.mixed)

## augment needs data, but will fail if data contains NAs in columns
##  not used in the model - extract formulas from model, use these
##  to find vars we need, extract those from data ...
get_data <- function(x,data) {
    if (inherits(x,"merMod")) {
        av <- all.vars(formula(x))
    } else {  ## gamm4/list
        av <- union(all.vars(formula(x$mer)),all.vars(formula(x$gam)))
    }
    return(data[intersect(av,names(data))])
    ## includes some bogus/constructed vars (y.0,X,Xr), but harmless?
}

## FIXME: still hard-coded
aa <- function(m,data=ecoreg) {
    augment(m,data=get_data(m,data))
}

tt <- function(m) {
    tidy(m)
}
## extract:
## model table - is.singular, glance() output
gg <- function(m) {
    m0 <- if (inherits(m,"merMod")) m else m$mer
    data.frame(glance(m),
               singular=is.singular(m0),
               df=lme4:::npar.merMod(m0))
}

get_allsum <- function(allfits) {
    pred <- (map(allfits,
                     ~ (map(.,aa) %>% bind_rows(.id="model")))
        %>% bind_rows(.id="taxon"))

    sum <- (map(allfits,
                    ~ (map(.,gg) %>% bind_rows(.id="model")))
        %>% bind_rows(.id="taxon")
        %>% select(taxon,model,AIC,singular,df)
        %>% group_by(taxon)
        %>% mutate(AIC=AIC-min(AIC,na.rm=TRUE),
                   AIC_OK=ifelse(singular,NA,AIC),
                   best=!is.na(AIC_OK) & AIC==min(AIC_OK,na.rm=TRUE))
        %>% arrange(AIC)
        %>% select(-AIC_OK)
    )

    ## all_sum %>% filter(best)

    coefs <- (map(allfits,
                      ~ (map(., tt) %>% bind_rows(.id="model")))
        %>% bind_rows(.id="taxon"))

    return(lme4:::namedList(coefs,sum,pred))
}

if (FALSE) {
    ## experimenting with AIC etc for gamm4 fits
    AIC(allfits[[1]][[1]]$mer)
    lme4:::npar.merMod(allfits[[1]][[1]]$mer) ## 17
    allfits[[1]][[1]]$gam$gcv.ubre
    ## sum(allfits[[1]][[1]]$) ## 17
}

source("mmd_utils.R")
load("ecoreg.RData")

system.time(load("allfits.RData"))
gamm4_res <- get_allsum(allfits)
save(gamm4_res,file="allfits_sum_gamm4.RData")
rm(allfits)

system.time(load("allfits_lmer.RData"))
lme4_res <- get_allsum(allfits)
save(lme4_res,file="allfits_sum_lmer.RData")
rm(allfits)

## testing
#####
library(ggplot2); theme_set(theme_bw())
library(cowplot)

load("allfits_sum_lmer.RData")
load("allfits_sum_gamm4.RData")

get_best_sum <- function(x) {
    x$sum %>%   ## tidyverse-ish extraction?
    group_by(taxon) %>%
    filter(best) %>%
    select(-c(AIC,best,singular))
}

get_best_name <- function(x) {
    ## extract *name* of best model
    best_model <- x$sum %>%
        ## couldn't figure out how to use 'taxon' twice here
        filter(taxon==tt) %>%
        filter(best) %>%
        pull(model)
    return(best_model)
}

get_best_sum(lme4_res)
get_best_sum(gamm4_res)

iwh14 <- c("#ff3579","#2ace05","#b90dd3","#00ac42","#4842b6","#e5c400","#1bd4ff","#cf2700","#018e66","#ff64b0","#47591c","#ff9065","#953012","#d5ad7f")
##    c("#984e81","#5bac47","#9d5ccf","#b2b045","#cb4fb2","#54a87f","#d7447d","#5ba2d6","#cc542b","#656fc8","#d79248","#d48ec9","#7f702f","#c85b5c")
##' @param x a model-summary object (lme4_res or gamm4_res)
##' @param taxon
diag_plot <- function(x,tt) {
    pred_vals <- x$pred %>%
        filter(taxon==tt,model==get_best_name(x))
    fitres <- ggplot(pred_vals,aes(.fitted,.resid))+geom_point()+
        geom_smooth()
    qqplot <- ggplot(pred_vals,aes(sample=.resid))+
        stat_qq(aes(colour=biome))+
        stat_qq_line()+
        scale_colour_manual(values=iwh14)
    plot_grid(fitres,qqplot)
}
