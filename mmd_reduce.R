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

get_best <- function(x) {
    x$sum %>%   ## tidyverse-ish extraction?
    group_by(taxon) %>%
    filter(best) %>%
    select(-c(AIC,best,singular))
}

get_best(lme4_res)
get_best(gamm4_res)

diag_plot <- function(x,tt) {
    best_model <- x$sum %>%
        ## couldn't figure out how to use 'taxon' twice here
        filter(taxon==tt) %>%
        filter(best) %>%
        pull(model)
    pred_vals <- x$pred %>%
        filter(taxon==tt,model==best_model)
    ggplot(pred_vals,
           aes(.fitted,.resid)
