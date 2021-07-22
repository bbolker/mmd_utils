library(tidyverse)
source("utils.R")
remotes::install_github("bbolker/r2glmm")
library(r2glmm)

get_best <- . %>% .[["sum"]] %>% filter(best) %>% dplyr::select(-c(best,singular))
get_best_tm <- . %>% get_best() %>% dplyr::select(taxon, model)
lme4_res <- readRDS("allfits_sum_lme4.rds")  ## 4 taxa x 27 fits, using lmer
lme4_best_models <- lme4_res %>% get_best_tm()

allfits_lme4 <- readRDS("allfits_lme4.rds")
best_models_lme4 <- apply(lme4_best_models, 1,
                          function(r) {
                            allfits_lme4[[r[["taxon"]]]][[r[["model"]]]]
                          })
names(best_models_lme4) <- names(allfits_lme4)[1:3] ## skip plants_log

best_models_gamm4 <- readRDS("bestmodels_gamm4.rds")[1:3] ## skip plants_log


get_rsq <- function(taxon,
                    models = best_models,
                    method = "kr",
                    df = NULL) {
  model <- models[[taxon]]
  if (is(model, "gamm4")) {
      if (is.null(df)) df <- model.frame(model$mer)
      ## hack area (not sure why other vars are already present
      ## in unmunged-name form??
      if (!"area_km2_log_sc" %in% names(df))
          df$area_km2_log_sc <- df$X.area_km2_log_sc
      names(df)[1] <- taxon
      formula <- reformulate(c("(NPP_log_sc+NPP_cv_sc + Feat_log_sc+Feat_cv_sc)^2", "area_km2_log_sc"), response = taxon)
      mm <- model.matrix(formula, df)
      r2beta(models[[taxon]], formula=formula, random=NULL, method=method, data=df,
             partial.terms=colnames(mm)[-1])
  } else {
      r2beta(models[[taxon]], method = method, partial = TRUE)
  }
}

std_order <- c("Area","NPP","NPP CV","Fire","Fire CV")
## create interaction terms
oo <- outer(std_order[-1],std_order[-1],paste,sep=":")
std_order2 <- c(std_order,oo[upper.tri(oo)]) ## append interactions to main effe

mutate_rsq <- (.  
    %>% as_tibble()
    %>% filter(taxon !="plants_log")
    %>%  mutate(Effect=as.character(Effect),
                Effect=gsub("^X","",Effect),
                Effect=trans_labs(Effect),
                Effect=replace_value_chr(Effect, "Fire:NPP CV","NPP CV:Fire"),
                Effect=factor(Effect,levels=rev(c("Model",std_order2)))
                )
    %>% dplyr::select(taxon,Effect,Rsq,upper.CL,lower.CL)
)

nn <- names(best_models_gamm4)
names(nn) <- nn ## ugh, for map_dfr
all_rsq_gamm4_kr <- (nn
    %>% purrr::map_dfr(get_rsq,
                       models = best_models_gamm4,
                       method = "kr",
                       .id="taxon")
    %>% mutate_rsq()
)

all_rsq_gamm4_sgv <- (nn
    %>% purrr::map_dfr(get_rsq,
                       models = best_models_gamm4,
                       method = "sgv",
                       .id="taxon")
    %>% mutate_rsq()
)

## FIXME: lost version also does lme4 fits for comparison (easy),
## also included components for grouped fire-related, NPP-related terms
## using stuff below, *or* SGV-based using cmp_R2() ?

## special_R2 <- function(model, type="fire") {
##   fire_terms <- "(Feat_log_sc * Feat_cv_sc) + (NPP_log_sc + NPP_cv_sc):(Feat_log_sc + Feat_cv_sc)"
##   NPP_terms <- "(NPP_log_sc * NPP_cv_sc) + (NPP_log_sc + NPP_cv_sc):(Feat_log_sc + Feat_cv_sc)"
##   ## update(formula(best_models[[1]]), as.formula(sprintf(". ~ . - (%s)", fire_terms)))
##   terms <- switch(type,
##                   fire=fire_terms,
##                   NPP=NPP_terms)
##   ## ugh
##   taxon <- deparse1(formula(model$gam)[[2]])
##   df <- model.frame(model$mer)
##   names(df)[1] <- taxon
##   formula <- reformulate("(NPP_log_sc+NPP_cv_sc + Feat_log_sc+Feat_cv_sc)^2", response = taxon)
##   r2beta(model, formula = formula, partial.terms = terms, method="KR", random=NULL,
##          data=df)
## }

## get_taxon_rsq <- function(t) {
##   (purrr::map_dfr(
##               best_models, .id="taxon", ~ special_R2(., t))
##     %>% as_tibble()
##     %>% filter(taxon !="plants_log")
##     %>% mutate(across(Effect, ~ ifelse(. == "Model", "Model", t)))
##     %>% dplyr::select(taxon,Effect,Rsq,upper.CL,lower.CL)
##   )
## }

save(list = ls(pattern = "all_rsq_"), file = "calc_rsq.rda")
