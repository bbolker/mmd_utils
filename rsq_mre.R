library(tidyverse)
## required hack for gamm4 models
remotes::install_github("bbolker/r2glmm") ## maybe not needed
library(r2glmm)
library(Matrix)
library(lme4)

## from r2glmm:::r2beta.lmerMod
extract.merMod.cov <- function(model) {
    Z = lme4::getME(model, "Z")
    s2e = lme4::getME(model, "sigma")^2
    lam = lme4::getME(model, "Lambda")
    lamt = lme4::getME(model, "Lambdat")
    G = s2e * (lam %*% lamt)
    SigHat = Z %*% (G %*% Matrix::t(Z))
    Matrix::diag(SigHat) = Matrix::diag(SigHat) + s2e/model@resp$weights
    return(SigHat)
}

## stored models
best_models <- readRDS("bestmodels_gamm4.rds")

## simplest example
## try SGV with different clustering levels
taxon <- "mbirds_log"
model <- best_models[[taxon]]$mer
X <- as.matrix(getME(model,"X"))
beta <- fixef(model)
p <- length(beta)
C_model <- cbind(0,diag(p-1))
data <- model.frame(model)
SigHat <- as.matrix(extract.merMod.cov(model))
Matrix::image(Matrix(SigHat))
res <- list()
for (clust.id in c("biome","flor_realms","biome_FR")) {
    ## skip Xr (fake grouping)
    for (collapse_sighat in c("blk_sgv","det_sgv","none")) {
        cat(clust.id, collapse_sighat, "\n")
        if (clust.id=="") {
            clust.id <- tail(names(model@flist),1)  ## replicate rsq default
        }
        obsperclust <- as.numeric(table(data[[clust.id]]))
        nclusts <- length(obsperclust)
        SigHat2 <- switch(collapse_sighat,
                          blk_sgv=calc_sgv(nblocks=nclusts, blksizes=obsperclust, vmat=SigHat),
                          det_sgv=exp(determinant(SigHat,log=TRUE)$modulus/nrow(SigHat)),
                          none=SigHat)
        cc <- cmp_R2(c=C_model, x=X, SigHat=SigHat2,
               beta=beta, method = "SGV",
               obsperclust=obsperclust, nclusts=nclusts)
        res <- c(res, list(data.frame(clust.id, collapse_sighat, rbind(cc))))
    }
}

all_res <- do.call("rbind", res)
print(all_res,digits=4)
r2beta(model, partial=FALSE) ## matches flor_realms
r2beta(model, partial=FALSE, method="nsj") ## matches flor_realms

r2beta(model, partial=FALSE, method="kr")

R2_vals_0 <- purrr::cross_df(list(taxon=c("mamm","bird","amph"),
                       method=c("sgv","nsj")))

pp <- pmap_dfr(R2_vals_0,
               function(taxon, method) {
                   r2beta(model = best_models[[taxon]]$mer,
                          method = method,
                          partial=FALSE)
               })

R2_vals <- R2_vals_0 %>% bind_cols(pp) %>% select(taxon, method,Rsq,lwr=upper.CL,upr=lower.CL)

ggplot(R2_vals, aes(x=Rsq,y=taxon,xmin=lwr,xmax=upr, colour=method)) +
    geom_pointrange(position=position_dodge(width=0.5))

#####

get_fire_R2 <- function(taxon, method="sgv", clust.id=NULL) {
    ## select merMod component of gamm4 fit
    model <- best_models[[grep(taxon, names(best_models))]]$mer
    X <- as.matrix(getME(model,"X"))
    beta <- fixef(model)
    if (is.null(clust.id)) {
        clust.id <- tail(names(model@flist),1)  ## replicate rsq default
    }
    data <- model.frame(model)
    obsperclust <- as.numeric(table(data[[clust.id]]))
    nclusts <- length(obsperclust)
    SigHat <- as.matrix(extract.merMod.cov(model))
    ## (what is this step doing ???)
    SigHat <- calc_sgv(nblocks=nclusts, blksizes=obsperclust, vmat=SigHat)
    p <- length(beta)
    C_model <- cbind(0,diag(p-1))
    nb <- gsub("^X","",names(beta))  ## clean up names (gamm4 hack)
    dimnames(C_model) <- list(nb[-1],nb)
    ## construct contrast matrix for terms involving fire ("Feat" in coef name)
    ## and *not* involving fire
    C_fire <- C_model
    C_fire <- C_fire[grep("Feat",names(beta))-1,]
    C_nonfire <- C_model
    C_nonfire <- C_model[grep("Feat",names(beta),invert=TRUE)-1,]
    R2 <- (purrr::map_dfr(list(model=C_model,fire=C_fire,nonfire=C_nonfire),
                         cmp_R2, x=X, SigHat=SigHat,
                         beta=beta, method = method,
                         obsperclust=obsperclust, nclusts=nclusts,
                         .id="effect")
        %>% mutate(
            lower.CL = stats::qbeta(0.025, v1/2, v2/2, ncp),
            upper.CL = stats::qbeta(0.975, v1/2, v2/2, ncp)
        ))
    return(R2)
}

gfr1 <- get_fire_R2("amph")
stopifnot(isTRUE(all.equal(gfr1$Rsq[1],0.5813195,
                           tolerance=1e-5)))
stopifnot(identical(gfr1, get_fire_R2("amph",method="sgv")),
          identical(gfr1,
                    get_fire_R2("amph",method="sgv",clust.id="flor_realms")))

R2_vals_1 <- purrr::cross_df(list(taxon=c("mamm","bird","amph"),
                       method=c("sgv","nsj"),
                       clust.id=c("biome_FR","biome","flor_realms")))

pp_1 <- pmap_dfr(R2_vals_1, get_fire_R2, .id="row")

all_R2 <- (mutate(R2_vals_1, row=as.character(seq(n)))
    %>% full_join(pp_1, by="row")
    %>% mutate(effect=factor(effect,levels=c("fire","nonfire","model")))
    %>% select(taxon,method,clust.id,effect,ncp,Rsq,lwr=lower.CL, upr=upper.CL)
)

print(ggplot(all_R2,aes(Rsq,taxon,colour=clust.id,shape=method))
      + facet_wrap(~effect,scale="free")
      + geom_pointrange(aes(xmin=lwr,xmax=upr),
                        position=position_dodge(width=0.25))
      )



## tackle KR for gamm4 models

ff <- attr(model@frame,"formula")
null_model <-
    ## don't include 'Xr' (smooths) in null model?
    lmer(update(as.formula(glmmTMB:::reOnly(ff)), y.0 ~ . - (1|Xr)),
         data=model.frame(model))
mc = pbkrtest::KRmodcomp(model, null_model)$stats
ss = with(mc, ndf * Fstat/ddf)
R2 = data.frame(Effect = "Model", F = mc$Fstat, v1 = mc$ndf,
        v2 = mc$ddf, pval = mc$p.value, ncp = mc$ndf * mc$Fstat,
        Rsq = ss/(1 + ss))



