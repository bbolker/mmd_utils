source("gamm4_utils.R")

##' basic check for singular random effects (i.e. overfitted)
##' @param fit fitted \code{merMod} object (from \code{lme4}, or
##' \code{$mer} component of a \code{gamm4} model)
is.singular <- function(fit,tol=1e-4) {
    if (inherits(fit,"brmsfit")) return(FALSE)
    return(any(abs(lme4::getME(fit,"theta"))<tol))
}

##' pretty-print model summary
##' @keywords internal
##' @param x a model
pfun <- function(x) print(summary(x),correlation=FALSE)
##
cplot <- function(x) {
    vv <- VarCorr(x)[[1]]
    cc <- cov2cor(vv)
    corrplot.mixed(cc,upper="ellipse")
}


get_best_name_fitlist <- function(fitlist,allow_sing=FALSE) {
    aa <- AICtabx(fitlist)
    if (!allow_sing) {
        aa <- aa[!aa$singular,]
    }
    return(rownames(aa)[which.min(aa$dAIC)])
}
##' @inheritParams plotfun
##' @param aux_quantiles (0.1, 0.5, 0.9) quantiles of auxiliary variable to predict
##' @param pred_lower_lim (-3) : lower cut off values (log scale)
##' @param re.form (NA) which RE to include in *predictions* (default is none)
##' @examples
##' source("gamm4_utils.R")
##' ecoreg <- readRDS("ecoreg.rds")
##' allfits_restr_gamm4 <- readRDS("allfits_restr_gamm4.rds")
##' m1 <- allfits_restr_gamm4$mbirds_log
##' pp <- predfun(m1)
##' pp2 <- predfun(m1,auxvar=NULL,grpvar="biome",re.form=~(1+NPP_log|biome))
##' library(ggplot2)
##' library(plotly)
##' gg0 <- ggplot(pp2,aes(NPP_log,mbirds_log,colour=biome))+geom_line(aes(group=biome))+
##'     geom_point(data=ecoreg,aes(shape=flor_realms))
##' ggplotly(gg0)
##'
predfun <- function(model=best_model,
                    data = ecoreg,
                    xvar="NPP_log_sc",
                    respvar=NULL,
                    auxvar="Feat_cv_sc",
                    grpvar=NULL,
                    pred_lower_lim= -3,
                    aux_quantiles=c(0.1,0.5,0.9),
                    re.form = NA,  ## exclude REs from prediction
                    alpha=0.05,
                    npts = 51,
                    adjust_othervar="mean",
                    exclude_fire=FALSE
                    ) {
    adjust_othervar <- match.arg(adjust_othervar)
    if (inherits(model,"brmsfit")) {
        require(brms)
        ff <- lme4::nobars(formula(model)$formula)
    } else if (inherits(model,"gamm4")) {
        require(gamm4)
        ## need x/y variables
        ff <- formula(model, fixed.only=TRUE, drop.smooth=FALSE)
    } else {
        ff <- formula(model, fixed.only=TRUE)
    }
    ## get LHS of formula
    mrespvar <- deparse(ff[[2]])
    if (is.null(respvar)) respvar <- mrespvar
    ## extract variables from data set by name
    xx <- data[[xvar]]
    ## add factor equiv of auxvar to data
    av <- all.vars(ff)
    othervars <- setdiff(av,c(respvar,xvar,auxvar,mrespvar,grpvar))
    ## construct prediction frame
    if (!is.null(auxvar)) {
        aa <- drop(data[[auxvar]])
        pdata <- expand.grid(seq(min(xx),max(xx),length=npts),
                             quantile(aa,aux_quantiles))
        names(pdata) <- c(xvar,auxvar)
    } else {
        if (!is.null(grpvar)) {
            ## grouped values
            lrange <- lapply(split(xx,data[[grpvar]]),
                             function(x) if (length(x)==0) NULL else (range(x,na.rm=TRUE)))
            tmpf <- function(nm,x) if (is.null(x)) NULL else data.frame(gg=nm,xx=x)
            pdata <- do.call(rbind,
                             Map(tmpf,names(lrange),lrange))
            names(pdata) <- c(grpvar,xvar)
        } else {
            ## xvar only
            pdata <- data.frame(seq(min(xx),max(xx),length=npts))
            names(pdata) <- xvar
        }
    }
    if (exclude_fire) {
        ## set fire to a minimum value (not zero!)
        ## Fire_cv to ???  (zero?)
        firevars <- grep("Feat",othervars,value=TRUE)
        othervars <- setdiff(othervars,firevars) ## exclude fire vars from othervars
        for (i in firevars) {
            ## set firevars to min val
            ok <- !is.na(data[[respvar]])
            cur_var <- data[[i]][ok]
            pdata[[i]] <- min(cur_var)
        }
    }
    ## variables other than primary x-variable and aux (and maybe grpvar) are set to median, or zero
    mfun <- get(adjust_othervar) ## mean, or median
    for (i in othervars) {
        ok <- !is.na(data[[respvar]])
        cur_var <- data[[i]][ok]
        if (adjust_othervar=="zero") {
            ## HACK: area_km2 is used in the formula as log(...)
            ## no longer necessary?
            if (i=="area_km2") {
                pdata[[i]] <- mfun(cur_var)
            } else pdata[[i]] <- 0
        } else {
            if (is.null(grpvar)) {
                pdata[[i]] <- mfun(cur_var)
            } else {
                agg <- aggregate(cur_var, by=list(data[[grpvar]][ok]),
                                 FUN=mfun)
                names(agg) <- c("biome", i)
                pdata <- merge(pdata, agg, by = "biome")
            }
        }
    }
    if (!is.null(auxvar)) {
        fauxvar <- paste0("f",auxvar)
        pdata[[fauxvar]] <- factor(pdata[[auxvar]],labels=paste0("Q(",aux_quantiles,")"))
    }
    if (inherits(model,"brmsfit")) {
        ## use fitted(), not predict(); leave out measurement error term
        ## brms needs all parameters specified
        if (is.na(re.form)) {
            pdata$biome <- pdata$flor_realms <- pdata$biome_FR <- NA
        }
        pp <- fitted(model,newdata=pdata,re.form=re.form)
        pred <- pp[,"Estimate"]
        pdata <- transform(pdata,
                           lwr=pp[,"Q2.5"],  upr=pp[,"Q97.5"])
        pdata[[mrespvar]] <- pred
    } else {
        pred <- predict(model,newdata=pdata,re.form=re.form)
        pdata[[mrespvar]] <- pred
        ## confidence intervals (fixed-effects only) on predictions
        mm <- model.matrix(formula(model),pdata)
        pvar1 <- diag(mm %*% tcrossprod(as.matrix(vcov(model)),mm))
        pdata <- transform(pdata,
                           lwr = qnorm(alpha/2,
                                       mean=pdata[[mrespvar]],sd=sqrt(pvar1)),
                           upr = qnorm(1-alpha/2,
                                       mean=pdata[[mrespvar]],sd=sqrt(pvar1)))


    }
    if (!is.na(pred_lower_lim)) {
        ## truncate predictions below, if requested
        cut_lwr <- function(x,val=NA) {
            x[x<pred_lower_lim] <- val
            return(x)
        }
        pdata[[mrespvar]] <- cut_lwr(pdata[[mrespvar]])
        if ("lwr" %in% names(pdata)) {
            pdata <- transform(pdata,
                               lwr=cut_lwr(lwr,val=pred_lower_lim),
                               upr=cut_lwr(upr,val=pred_lower_lim))
        }
    }
    ## assign scaling/logging info to pred values
    pdata[[mrespvar]] <- copy_attributes(data[[mrespvar]],pdata[[mrespvar]])
    pdata$lwr <- copy_attributes(data[[mrespvar]],pdata$lwr)
    pdata$upr <- copy_attributes(data[[mrespvar]],pdata$upr)
    return(pdata)
}

    auto_lab_text <- c("NPP"="NPP~(g~C%.%m^{-2}%.%year^{-1})",
                       "NPP_cv"="Interannual~CV~of~NPP",
                       "Feat"="Fire~loss~(prop.~NPP~year^{-1})",
                       "Feat_cv"="Interannual~CV~of~fire~loss",
                       "mbirds"="Birds~(species/km^2)",
                       "mamph"="Amphibians~(species/km^2)",
                       "mmamm"="Mammals~(species/km^2)",
                       "plants"="Plants~(species/km^2)",
                       "area_km2"="Ecoregion~area~(km^{2})")

## without units - doesn't need to be an expression so change ~ back to space
auto_lab_text_nounits <- gsub("~"," ",
                              gsub("~[(].*$","",auto_lab_text))

auto_lab_short <- c(NPP="NPP",NPP_cv="NPP CV",
                    Feat="Fire",Feat_cv="Fire CV",
                    area_km2="Area")

test <- c("NPP_log_sc", "Feat_log_sc", "NPP_cv_sc", "Feat_cv_sc", "area_km2_log_sc",
"NPP_log_sc:Feat_log_sc", "NPP_log_sc:NPP_cv_sc", "NPP_log_sc:Feat_cv_sc",
"Feat_log_sc:NPP_cv_sc", "Feat_log_sc:Feat_cv_sc", "NPP_cv_sc:Feat_cv_sc")
trans_labs <- function(x) {
    ## match name followed by _, :, end of line
    ## (?!_) is a negative lookahead; exclude stuff followed by _,
    ## and don't count it for substitution purposes
    w <- function(y) sprintf("%s(?!_)",y)
    x <- gsub("_(log|sc)","",x)
    for (i in seq_along(auto_lab_short)) {
        x <- gsub(w(names(auto_lab_short)[i]),auto_lab_short[i],x,
                  perl=TRUE)
    }
    return(x)
}

##' @param model fitted model
##' @param data data frame containing values
##' @param xvar  x-variable
##' @param respvar (equal to model response by default): response variable
##' @param auxvar ("Feat_cv_sv"): auxiliary variable (e.g. for examining interactions)
##' @param backtrans back-transform x and y variables?
##' @param log character string, as in \code{\link{plot.default}}, that specifies
##' whether to log-scale the x and y axes ("xy"), one or the other ("x" or "y"),
##' or neither ("")
##' @param xlim x-limits: if missing, use range of data (not predictions/CI)
##' @param ylim y-limits: ditto
##' @param ... parameters passed through to predfun
##' @examples
##' source("gamm4_utils.R")
##' ecoreg <- readRDS("ecoreg.rds")
##' allfits_restr_gamm4 <- readRDS("allfits_restr_gamm4.rds")
##' m1 <- allfits_restr_gamm4$mbirds_log
##' plotfun(m1)
##' plotfun(m1,auxvar=NULL)
plotfun <- function(model=best_model,
                    data = ecoreg,
                    xvar="NPP_log_sc",
                    respvar=NULL,
                    auxvar="Feat_cv_sc",
                    grpvar=NULL,
                    xlim = NULL,
                    ylim = NULL,
                    backtrans=FALSE,
                    lty=c(2,1,3),
                    log="",
                    adjust_othervar="mean",
                    auto_label=(backtrans==TRUE),
                    sci10_axis=NULL,
                    return_data = FALSE,
                    ...
                    ) {
    require("ggplot2")
    ## response variable from model (may not be the same as respvar,
    ##  e.g. in partial residuals plots)
    if (!is.null(model)) {
        form <- formula(model)
        if (inherits(model,"brmsfit")) {
            form <- form$formula
        }
        mrespvar <- deparse(form[[2]])
        if (is.null(respvar)) respvar <- mrespvar
        if (respvar=="partial_res") {
            data$partial_res <- remef_allran(model,
                                             data=data,
                                             fixed_keep=c("(Intercept)",xvar))
        }
        pdata <- predfun(model,data=data,xvar=xvar,respvar=respvar,
                         auxvar=auxvar,grpvar=grpvar,
                         adjust_othervar=adjust_othervar,...)
    } else {
        ## no model, must specify respvar
        mrespvar <- respvar
    }
    if (backtrans) {
        ## back-transform response variable from model,
        ##  with CIs
        if (!is.null(model)) {
            pdata <- backtrans_var(pdata, mrespvar, data,
                                   othervars=c("lwr","upr"), log=TRUE)
            ## change name to transformed name (e.g. minus "_log")
            mrespvar <- attr(pdata,"transvar")
            ## back-transform x-variable
            pdata <- backtrans_var(pdata,xvar,data)
            ## back-tranform respvar in data
        }
        data <- backtrans_var(data,respvar)
        respvar <- attr(data,"transvar")
        if (is.null(model)) mrespvar <- respvar
        ## back-transform xvar in data
        data <- backtrans_var(data,xvar)
        xvar <- attr(data,"transvar")
    }
    if (return_data) {
      data[["focal"]] <- data[[xvar]]
      data[["richness"]] <- data[[respvar]]
      pdata[["focal"]] <- pdata[[xvar]]
      pdata[["richness"]] <- pdata[[mrespvar]]
      return(list(data = data, pdata = pdata, respvar = respvar, xvar = xvar))
    }
    gg0 <- ggplot(data,aes_(x=as.name(xvar)))
    if (!is.null(model)) {
        if (!is.null(auxvar)) {
            fauxvar <- paste0("f",auxvar)
            ## geom_encircle(aes(group=biome_FR),expand=0)+  ## ugly ...
            ## use model respvar for predicted values
            gg0 <- gg0 + geom_line(data=pdata,aes_(y=as.name(mrespvar),
                                                   linetype=as.name(fauxvar))) +
                geom_ribbon(data=pdata,aes_(ymin=~lwr,ymax=~upr,
                                            group=as.name(fauxvar)),
                            colour=NA,fill="black",alpha=0.1)+
                scale_linetype_manual(values=lty)
        } else {
            if (!is.null(grpvar)) {
                gg0 <- gg0 + geom_line(data=pdata,aes_(y=as.name(mrespvar),
                                                       colour=as.name(grpvar))) +
                    geom_ribbon(data=pdata,
                                aes_(ymin=~lwr,ymax=~upr,
                                     fill=as.name(grpvar)),
                                alpha=0.1,colour=NA)
            } else {
                gg0 <- gg0 + geom_line(data=pdata,aes_(y=as.name(mrespvar))) +
                    geom_ribbon(data=pdata,aes_(ymin=~lwr,ymax=~upr),
                                colour=NA,fill="black",alpha=0.1)
            }
        }
    }  ## if model
    ## finish
    labels_x <- waiver()
    ## fanciness for appropriate default
    if ((missing(sci10_axis) &&
         xvar %in% c("Feat","area_km2")) || isTRUE(sci10_axis)) {
        labels_x <- scientific_10
    }
    if (grepl("x",log)) {
            gg0 <- gg0 + scale_x_log10(labels=labels_x)
    }
    if (grepl("y",log)) {
            gg0 <- gg0 + scale_y_log10()
    }
    ## set x, y to range of *data* (not prediction/CIs) if
    ##   missing (explicit NULL recovers entire range)
    if (missing(xlim)) {
        xlim <- range(data[[xvar]],na.rm=TRUE)
    }
    if (missing(ylim)) {
        ylim <- range(data[[respvar]],na.rm=TRUE)
    }
    gg0 <- (gg0
        + theme(legend.box="horizontal")
        + coord_cartesian(xlim=xlim, ylim=ylim)
    )
    ## https://stackoverflow.com/questions/18237134/line-break-in-expression
    ## (but doesn't work)
    pfun <- function(x) {
        p <- parse(text=x)
        if (grepl("\n",x)) {
            p <- as.list(p)
        }
        return(p)
    }
    if (auto_label) {
        gg0 <- gg0 + labs(x=pfun(auto_lab_text[[xvar]]),
                          y=pfun(auto_lab_text[[mrespvar]]))
    }
    return(gg0 + geom_point(aes_(y=as.name(respvar),
                                colour=~biome,shape=~flor_realms)))
}

## test

pkgList <- c(
    'plyr'       ## manipulation -- LOAD THIS FIRST
  , 'Cairo'
  ,'TMB'
  ,'bbmle'       ## AICtab cosmetic
  ,'brms'
  ,'broom.mixed' ## better coef tables
  ,'colorspace'
  ,'cowplot'
  ,'data.table' ##
  ,'dplyr'      ## data manipulation
  ,'fields'
  ,'gamm4'      ## spatial fits
  ,'ggplot2'
  ,'ggstance'   ## horizontal geoms
  ,'gridExtra'  ## arrange plots
  ,'lattice'    ## diagnostic plots
  ,'lme4'       ## lmer etc.
  ,'plotly'
  ,'plotrix'
  ,'r2glmm'
  ,'raster'
  ,'remef'      ## remotes::install_github('https://github.com/hohenstein/remef')
  ,'remotes'
  ,'rgdal'
  ,'sp'
  ,'tibble'     ## ditto: rownames_to_column
  ,'tidyr'      ## ditto
  ,'viridis'
  , 'patchwork'
  , 'ggrepel'
  , 'directlabels'
  , 'hues'
  , 'pander'
  , 'ggplotify'
)

install_all_pkgs <- function() {
    i1 <- installed.packages()
    git_pkgs <- "remef" ## FIXME: complete
    to_install <- setdiff(pkgList, c(rownames(i1),git_pkgs))
    if (length(to_install)>0) {
        install.packages(to_install)
    }
    remotes::install_github("bbolker/r2glmm")
    if (! "remef" %in% rownames(i1)) {
        remotes::install_github('https://github.com/hohenstein/remef')
    }
}

load_all_pkgs <- function() {
    if (!require("remef", quietly=TRUE)) {
        stop("install remef via remotes::install_github('hohenstein/remef')")
    }
    sapply(pkgList,library,character.only=TRUE)
    stopifnot(packageVersion("lme4")>="1.1-14")
    ## devtools::install_github("hohenstein/remef")
    ## FIXME: print results
}

shorten_modelname <- function(mname) {
    return(gsub("diag","d",
                gsub("full","f",
                     gsub("int","i",
                          gsub("flor_realms","fr",
                               gsub("biome","b",mname))))))
}

## extraction from summaries
get_best_sum <- function(x) {
    x$sum %>%   ## tidyverse-ish extraction?
    group_by(taxon) %>%
    filter(best) %>%
    dplyr::select(-c(AIC,best,singular))
}

get_best_name <- function(x,tt) {
    ## extract *name* of best model
    best_model <- x$sum %>%
        ## couldn't figure out how to use 'taxon' twice here
        filter(taxon==tt) %>%
        filter(best) %>%
        pull(model)
    return(best_model)
}

get_best_pred <- function(x,tt,best=get_best_name(x,tt)) {
    p <- x$pred %>%
        filter(taxon==tt,model==best)
    return(p)
}

## tried to get 14 distinguishable colours from IWantHue
## (still not great)
iwh14 <- c("#ff3579","#2ace05","#b90dd3","#00ac42","#4842b6","#e5c400","#1bd4ff","#cf2700","#018e66","#ff64b0","#47591c","#ff9065","#953012","#d5ad7f")
##' @param x a model-summary object (lme4_res or gamm4_res)
##' @param taxon
diag_plot <- function(x,tt="plants_log") {
    require(cowplot)
    pred_vals <- get_best_pred(x,tt)
    fitres <- ggplot(pred_vals,aes(.fitted,.resid))+geom_point()+
        geom_smooth()
    qqplot <- ggplot(pred_vals,aes(sample=.resid))+
        stat_qq(aes(colour=biome))+
        stat_qq_line()+
        scale_colour_manual(values=iwh14)
    cowplot::plot_grid(fitres,qqplot)
}

get_allcoefs <- function(data,model,focal_taxon="plants_log") {
    (data$coef
        ## add summary info (singular, AIC)
        %>% full_join(dplyr::select(data$sum,c(taxon,model,singular,AIC)),
                      by=c("taxon","model"))
        ## pick one taxon; only fixed-effect parameters; drop intercept
        %>% filter(taxon==focal_taxon,
                   effect=="fixed",
                   term!="(Intercept)")
        %>% mutate(model=shorten_modelname(model),     ## abbreviate model name
                   model=reorder(model,-AIC),          ## reverse-order by AIC
                   ## reverse-order by magnitude of coefficient for best model
                   term=reorder(term,estimate,function(x) -x[1]))
    )
}

## utilities for post-processing coefficient tabs
add_wald_ci <- function(data) {
    return(data %>%
           mutate(lwr=estimate-1.96*std.error,upr=estimate+1.96*std.error) %>%
           dplyr::select(-c(std.error,statistic)))
}
drop_intercept <- function(data) {
    data %>% filter(term != "(Intercept)")
}

##' generate coefficients relevant to each observation
##' @examples
##' load("ecoreg.RData")
##' load("bestmodels_gamm4.RData")
##' merge_coefs(ecoreg,best_models[[1]])
merge_coefs <- function(data,model,id_vars=c("biome","biome_FR","flor_realms"),
                        extra_vars=c("x","y","eco_id")) {
    rr <- ranef(model)
    rr <- rr[names(rr)!="Xr"] ## spatial term from gamm4
    ff <- fixef(model)
    ## combine ID vars with replicated fixed effects
    coefs <- data.frame(data[c(id_vars,extra_vars)],
                        matrix(ff,
                               byrow=TRUE,
                               ncol=length(ff),nrow=nrow(data),
                               dimnames=list(rownames(data),names(ff))),
                        check.names=FALSE)
    for (i in seq_along(rr)) {
        rrc <- rr[[i]] ## current ranef
        vars <- names(rrc)
        rrc[[names(rr)[[i]]]] <- rownames(rrc)
        dd <- merge(data[id_vars],rrc)
        for (j in vars) {
            coefs[,j] <- coefs[,j] + dd[,j]
        }
    }
    return(coefs)
}


##' Partial residuals, dropping all random variables
##' @param x model
##' @param data original data
##' @param fixed_keep which fixed effects should be \strong{retained}
##' (\code{NULL} means to keep all fixed effects
##' @param na.action what to do with NA values
##' @param return_components
##' @details see gamm4_predict.{rmd,html} for more info about which
##'  components (fixed, random, smooth) are included in which bits
##'  of a gamm4 model, how to access them via predict/residuals
##
## I think I could use that information to get partial residuals
##   for a subcomponent (i.e. using newdata=), but for now I'm
##   going to use residuals (which subtracts everything (fixed/random/smooth)
##   and add specified fixed effects back in ...
remef_allran <- function(x, data,
                         fixed_keep="(Intercept)",
                         set_other=c("mean","median","zero"),
                         na.action=na.exclude,
                         return_components=FALSE) {
    set_other <- match.arg(set_other)
    require(mgcv) ## for s()
    if (!inherits(x,"gamm4")) stop("only implemented for gamm4")
    ## construct fixed-effect model frame
    ds <- glmmTMB:::drop.special
    dh <- glmmTMB:::dropHead
    ## random-effects formula, minus artificial parts
    ff <- ds(formula(x$mer,random.only=TRUE),quote((1|Xr)))
    ## RHS only
    ff <- ff[-2]
    ## fixed term, minus smooth terms
    ff_fixed <- dh(formula(x$gam),quote(s))
    mf_fixed <- model.frame(ff_fixed,data,na.action=na.action)
    na.act <- attr(mf_fixed,"na.action")
    rr <- residuals(x$mer)
    ## predict *all* ranef (not really necessary except for return_components)
    pp_ran <- predict(x$mer,random.only=TRUE,re.form=ff)
    ## predict only specified FE
    pp_fixed <- rep(0,length(rr))
    mm_fixed <- model.matrix(ff_fixed[-2],data)
    ## hack/remove response var (delete.response() ?)
    mm_fixed <- mm_fixed[,colnames(mm_fixed)!="y"]
    cc <- coef(x$gam)
    cc <- cc[intersect(names(cc),colnames(mm_fixed))]
    if (length(fixed_keep)>0) {
        if (!all(is.na(fixed_keep))) {
            ## if NA, don't zero out any variables (="keep none")
            ## NULL = "keep all"
            ## DO want to predict
            ## mm_fixed[,!colnames(mm_fixed) %in% fixed_keep] <- 0
            ## set to median or mean, not zero ... ?
            for (col in setdiff(colnames(mm_fixed),fixed_keep)) {
                mm_fixed[,col] <- switch(set_other,
                         mean=mean(mm_fixed[,col],na.rm=TRUE),
                         median=median(mm_fixed[,col],na.rm=TRUE),
                         zero=0)
            }
        }
        pp_fixed <- mm_fixed %*% cc
        if (length(na.act)>0)  {
            ## drop NA values (even those in response), in order for preds to match residuals
            pp_fixed <- pp_fixed[-na.act]
        }
    }
    ## rem <- napredict(na.act,resp-pp_ran-pp_fixed)
    rem <- napredict(na.act,rr + pp_fixed)
    if (return_components) {
        return(data.frame(rem,pp_ran,pp_fixed))
    }
    ## transfer scaling/transformation info from response var to partial resids
    ## DRY?  did we extract the formula already above?
    respvar <- deparse(formula(x, fixed.only=TRUE)[[2]])
    rr <- data[[respvar]]
    rem <- copy_attributes(rr,rem)
    return(rem)
}


## 1. compute residuals
## 2. compute prediction with *only* 'keep' variables
##   (predict(.,terms=...) if we had it)
##  remove intercept???
## 3. add pred back to residuals

## or:

## compute prediction with all *but* 'keep' variables
## subtract from observed values
remef_2 <- function(x,data,response,
                       fixed_keep="NPP_log") {
    ## predict everything *but*
    data[,fixed_keep] <- 0
    pp <- predict(x,newdata=data) ## or fitted()
    return(response-pp)
}


##' undo the effects of scale()
##' @param x primary variable to unscale
##' @param y optionally: variable containing unscaling information
##' @param warn.noscale warn if there is no scaling info?
##' @examples
##' x <- 1:3
##' x <- drop(scale(1:3))
##' y <- 2:4
##' backtrans(x)
##' backtrans(y,x)
##' head(backtrans(ecoreg$mmamm_log))
backtrans <- function(x,y=NULL,warn.noscale=FALSE) {
    ## if only x is specified, x is both the target and has unscaling info
    ## if x and y are specified, x is the target and y has unscaling info
    scvar <- if (is.null(y)) x else y
    ctr <- attr(scvar,"scaled:center")
    sc <- attr(scvar,"scaled:scale")
    ## no scaling information: (silently) return original value
    if (is.null(ctr) || is.null(sc)) {
        if (warn.noscale) warning("no scaling info, returning original var")
        return(x)
    }
    ## remove scaling information from result
    ret <- x*sc+ctr
    attr(ret,"scaled:scale") <- NULL
    attr(ret,"scaled:center") <- NULL
    return(ret)
}

##' back-transform a variable by exponentiating, based on name
##' attempt to unscale (only effective if unscaling info
##'  is stored in x or y)
##' @note data are logged before scaling; to invert this
##' we need to unscale \strong{before} exponentiating ...
##' @examples
##' head(ecoreg$mmmamm_log)
##' head(ecoreg$mmamm)
##' backtrans_magic(ecoreg$mmamm_log,"mmamm_log")
backtrans_magic <- function(x,xname,y=NULL,log=NULL) {
    r <- backtrans(x,y)  ## unscale
    re <- "(_log)?(_sc)?$" ## _log OR _sc OR both (in that order)
    ## exponentiate if var is logged or if forced
    if (isTRUE(log) || grepl(re,xname) || isTRUE(attr(r,"logged"))) {
        r <- exp(r)
        attr(r,"logged") <- NULL
    }
    ## attempt
    attr(r,"name") <- gsub(re,"",xname)
    return(r)
}

copy_attributes <- function(x,y,attrs=c("logged","scaled:scale","scaled:center")) {
    for (a  in attrs) {
        attr(y,a) <- attr(x,a)
    }
    return(y)
}

##' back-transform one or more variables
##' @param data  data frame to back-transform
##' @param xname response variable
##' @param otherdata original data set containing variables with scaling info
##' @param othervars other variables to transform with the same scaling
## what happens if we pass "rem1"?  How are we supposed to back-transform it?
## should be back-transformed according to the same rules as
## the original response variable ...
backtrans_var <- function(data, xname,otherdata=NULL,othervars=NULL,...) {
    if (xname %in% names(data)) {
        r1 <- backtrans_magic(data[[xname]],xname,
                              otherdata[[xname]],...)
    }
    ## other variables (e.g. lower/upper CIs);
    ##  transform in place, using same scaling info
    for (v in othervars) {
        data[[v]] <- backtrans_magic(data[[v]],v,
                                     otherdata[[xname]],...)
    }
    ## extract new name, store in data
    new_var <- attr(r1,"name")
    data[[new_var]] <- r1
    ## attach info about transformed var
    attr(data,"transvar") <- new_var
    return(data)
}

## fix_NAs needed for remef() results: remef_allran does it automatically
fix_NAs <- function(rem,model) {
    if (!is.null(nastuff <- attr(model.frame(model),"na.action"))) {
        return(napredict(nastuff,rem))
    } else return(rem)
}

scientific_10 <- function(x, suppress_ones=TRUE, log10_thresh = 2) {
    if (all(abs(log10(na.omit(x))) < log10_thresh)) return(parse(text = x))
    s <- scales::scientific_format()(x)
    ## substitute for exact zeros
    s[s=="0e+00"] <- "0"
    ## regex: [+]?  = "zero or one occurrences of '+'"
    s2 <- gsub("e[+]?", " %*% 10^", s )
    ## suppress 1 x
    if (suppress_ones) s2 <- gsub("1 %\\*% +","",s2)
    parse(text=s2)
}

## mfun <- function(x) {
##     s <- scales::scientific_format()(x)
##     ss <- strsplit(s,"e")
##     base <- sapply(ss,"[[",1)
##     pow <- sapply(ss,"[[",2)
##     substitute(expression(x %*% 10^y),list(x=base,y=pow))
## }

graphics_setup <- function() {
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(colorspace))
    theme_set(theme_bw())
    ## utilities for cleaning up plots
    assign("zmargin",theme(panel.spacing=grid::unit(0,"lines")),envir=.GlobalEnv)
    assign("zmarginx",theme(panel.spacing.x=grid::unit(0,"lines")),envir=.GlobalEnv)
    assign("nolegend",theme(legend.position="none"),envir=.GlobalEnv)
    ## override default colors
    assign("scale_colour_discrete",
           function(...) scale_colour_discrete_qualitative(...),envir=.GlobalEnv)
}

## partial residuals plots
remef_plot <- function(taxon="mbirds_log",xvar="NPP_log",
                       auxvar=NULL,title=NULL) {
    m <- best_models[[taxon]]
    if (is.null(title)) {
        title <- if (is.null(auxvar)) xvar else {
             paste(xvar,auxvar,sep=":")
                                              }
    }
    pp <- (plotfun(m,xvar=xvar,respvar="partial_res",
                   auxvar=auxvar,data=ecoreg,
                   backtrans=TRUE,log="xy")
      + theme(legend.position="none")
      + ggtitle(title)
    )
    return(pp)
}

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

## https://www.labnol.org/internet/direct-links-for-google-drive/28356/
get_googledrive <- function(gid,dest) {
    download.file(sprintf("https://drive.google.com/uc?export=download&id=%s",gid),
                  destfile=dest)
}


read.dot <- function (f) {
    lines <- readLines(f)
    body <- lines[grep("->", lines, fixed = TRUE)]
    nodePairs <- sub("^[[:space:]]+\"", "\"", sub("\"[;[:space:]]+$",
        "\"", unlist(strsplit(body, "->"))))
    nodeLists <- split(nodePairs, 1:length(nodePairs)%%2)
    nodes <- unique(nodePairs)
    edges <- data.frame(orig = nodeLists[[2]], dest = nodeLists[[1]])
    n <- length(nodes)
    graph <- matrix(0, n, n, dimnames = list(nodes, nodes))
    for (node in nodes) {
        graph[node, nodes %in% edges$dest[edges$orig == node]] <- 1
    }
    ## retrieve names
    labels <- lines[grep("label=", lines, fixed=TRUE)]
    return(graph)
}

replace_value_chr <- function(x, from, to) {
    replace(as.character(x), which(x==from), to)
}


figzip <- function() {
    system(sprintf("zip figs_%s.zip figures.html fig*.pdf fig*.png",
                   format(Sys.time(), "%Y%m%d")))
}
