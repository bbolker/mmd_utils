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
        ## include x/y vars for future ref even though not in formula
        av <- c(all.vars(formula(x)),"x","y")
    } else if (inherits(x,"brmsfit")) {
        av <- c(all.vars(formula(x$formula)),"x","y")
    } else {  ## gamm4/list
        av <- union(all.vars(formula(x$mer)),all.vars(formula(x$gam)))
    }
    return(data[intersect(av,names(data))])
    ## includes some bogus/constructed vars (y.0,X,Xr), but harmless?
}

## FIXME: ecoreg data set still hard-coded here
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
    ## get predictions from all taxa, bind them together
    pred <- (map(allfits,
                     ~ (map(.,aa) %>% bind_rows(.id="model")))
        %>% bind_rows(.id="taxon"))

    ## get summaries from all taxa, bind them together
    sum <- (map(allfits,
                    ~ (map(.,gg) %>% bind_rows(.id="model")))
        %>% bind_rows(.id="taxon")
        %>% select(taxon,model,AIC,singular,df)
        %>% group_by(taxon)
        ## find delta-AIC and identify best non-singular model
        %>% mutate(AIC=AIC-min(AIC,na.rm=TRUE),
                   AIC_OK=ifelse(singular,NA,AIC),
                   best=!is.na(AIC_OK) & AIC==min(AIC_OK,na.rm=TRUE))
        %>% arrange(AIC)
        %>% select(-AIC_OK)
    )

    ## get coefficients from all taxa, bind them together
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

system.time(L <- load("allfits.RData"))  ## 25 seconds
## "allfits": length-4, taxa;
##                length-27, all models
gamm4_res <- get_allsum(allfits)
save(gamm4_res,file="allfits_sum_gamm4.RData")

## find best gamm4 model for each taxon,
taxa <- names(allfits)
best_names <- map(taxa,get_best_name,x=gamm4_res)
best_models <- map2(allfits,best_names, ~.x[[.y]])
best_models <- map(best_models, strip_gamm4_env)
save(best_models,file="bestmodels_gamm4.RData")

gamm4_allfits <- lapply(allfits,
            function(x) lapply(x,strip_gamm4_env))
## save(gamm4_allfits,file="allfits_strip_gamm4.RData")
## ugh. "only" 160 M

rm(allfits)
gc()

## save test fits, for queries to mailing lists/Simon Wood/etc.
gamm4_testfits <- gamm4_allfits[[1]][1:8]
save(gamm4_testfits,file="testfits.RData")  ## down to 11M

system.time(load("allfits_lmer.RData"))
lme4_res <- get_allsum(allfits)
save(lme4_res,file="allfits_sum_lmer.RData")
rm(allfits)

