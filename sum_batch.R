## extract useful/necessary information
library(purrr)
library(dplyr)
library(broom.mixed)
library(gamm4)
library(brms)

source("utils.R")

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
    tidy(m,conf.int=TRUE)
}
## extract:
## model table - is.singular, glance() output
gg <- function(m) {
    m0 <- if (inherits(m,"gamm4")) m$mer else m 
    singular <- if (inherits(m,"brmsfit")) NA else is.singular(m0)
    df <- if (inherits(m,"brmsfit")) NA else lme4:::npar.merMod(m0)
    data.frame(glance(m),singular,df)
}

get_allsum <- function(allfits, nested=TRUE, debug=FALSE) {
    ## FIXME: would like to use purrr::map but don't know enough about NSE
    if (nested) {
         mfun <- function(FUN) {
             lapply(allfits,
                 function(z) (lapply(z,FUN) %>% bind_rows(.id="model")))
         }
     ## map(allfits, ~ (map(.,FUN)) %>% bind_rows(.id="model"))
    } else {
        mfun <-function(FUN) {
            lapply(allfits, FUN)
        }
    }
    ## get predictions from all taxa, bind them together
    ## FIXME: what are ending *_log columns (all NA)?
    if (debug) cat("extracting predictions\n")
    pred <- mfun(aa) %>% bind_rows(.id="taxon")
    ## get summaries from all taxa, bind them together
    if (debug) cat("extracting summaries\n")
    sum0 <- suppressWarnings(mfun(gg) 
    	%>% bind_rows(.id="taxon")
        %>% select(one_of(c("taxon","model","AIC","singular","df")))
        ## %>% select(-c(pss,nobs))
        %>% group_by(taxon))
    if ("AIC" %in% names(sum0)) {
        sum0 <- (sum0 
            ## find delta-AIC and identify best non-singular model
            %>% mutate(AIC=AIC-min(AIC,na.rm=TRUE),
                       AIC_OK=ifelse(singular,NA,AIC),
                       best=!is.na(AIC_OK) & AIC==min(AIC_OK,na.rm=TRUE))
            %>% arrange(AIC)
            %>% select(-AIC_OK)
        )
    }

    ## get coefficients from all taxa, bind them together
    if (debug) cat("extracting coefficients\n")
    coefs <- (mfun(tt) 
        %>% bind_rows(.id="taxon")
    )

    return(lme4:::namedList(coefs,sum=sum0,pred))
}

if (FALSE) {
    ## experimenting with AIC etc for gamm4 fits
    AIC(allfits[[1]][[1]]$mer)
    lme4:::npar.merMod(allfits[[1]][[1]]$mer) ## 17
    allfits[[1]][[1]]$gam$gcv.ubre
    ## sum(allfits[[1]][[1]]$) ## 17
}

ecoreg <- readRDS("ecoreg.rds")

system.time(allfits_gamm4 <- readRDS("allfits_gamm4.rds"))
## "allfits": length-4, taxa;
##                length-27, all models
gamm4_res <- get_allsum(allfits)
saveRDS(gamm4_res,file="allfits_sum_gamm4.rds")

## find best gamm4 model for each taxon,
taxa <- names(allfits_gamm4)
best_names <- map(taxa,get_best_name,x=gamm4_res)
best_models <- map2(allfits_gamm4,best_names, ~.x[[.y]])
best_models <- map(best_models, strip_gamm4_env)
saveRDS(best_models,file="bestmodels_gamm4.rds")

## gamm4_allfits <- lapply(allfits,
##             function(x) lapply(x,strip_gamm4_env))
## save(gamm4_allfits,file="allfits_strip_gamm4.RData")
## ugh. "only" 160 M

## save test fits, for queries to mailing lists/Simon Wood/etc.
## gamm4_testfits <- allfits_gamm4[[1]][1:8]
## saveRDS(gamm4_testfits,file="testfits.rds")  ## down to 11M

rm(allfits_gamm4)
gc()

system.time(allfits_lme4 <- readRDS("allfits_lme4.rds"))
lme4_res <- get_allsum(allfits_lme4)
saveRDS(lme4_res,file="allfits_sum_lme4.rds")

system.time(allfits_brms <- readRDS("allfits_brms.rds")) ## 6 seconds
brms_res <- get_allsum(allfits_brms,nested=FALSE)
saveRDS(brms_res, file="allfits_sum_brms.rds")

