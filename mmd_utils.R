source("gamm4_utils.R")

## basic check for singular random effects (i.e. overfitted)
is.singular <- function(fit,tol=1e-4) {
    any(abs(lme4::getME(fit,"theta"))<tol)
}
## pretty-print model summary
pfun <- function(x) print(summary(x),correlation=FALSE)
##
cplot <- function(x) {
    vv <- VarCorr(x)[[1]]
    cc <- cov2cor(vv)            
    corrplot.mixed(cc,upper="ellipse")
}

## fits all of the combinations of random effects for a particular
## response variable
##
##' @param response (character) response variable
##' @param pars (character) vector of input variables
##' @param forms list of allowed random-effect terms
##' @param rterms grouping variables
##' @param interax include (2-way) interactions among fixed effects?
##' @param data data frame
##' @param single_fit fit a single specified model instead of all combinations?
##' @param use_gamm4 fit a spatial model?
##' 
##' @examples
##' # two equivalent ways to fit a single model with diagonal terms
##' #  for biome and flor_realms and no effect for their interaction
##' f1 <- fit_all(single_fit=c(2,2,NA))
##' f2 <- fit_all(single_fit=c(2,2),rterms=c("biome","flor_realms"))
## TODO
## - select terms by name ("int", "diag", "full") rather than number?
## - repeat less (if lmer() is encapsulated may need to play with
##   environment stuff some more, ugh)
fit_all <- function(response="mbirds_log",
                    pars=c("NPP_log","Feat_log","NPP_cv_sc","Feat_cv_sc"),
                    extra_pred_vars="log(area_km2)",
                    ## possible random-effect models
                    forms=c(int="1|",  ## intercept-only
                            diag=paste("1+",paste(pars,collapse="+"),"||"),
                            full=paste("1+",paste(pars,collapse="+"),"|")
                            ),
                    ## grouping variables
                    rterms=c("biome","flor_realms","biome_FR"),
                    interax=TRUE,
                    data=ecoreg,
                    single_fit=NULL,
                    use_gamm4=FALSE,
                    verbose=FALSE) {
    if (use_gamm4) {
        ## add spherical smoothing term
        ## ?smooth.construct.sos.smooth.spec:
        ##  s(latitude,longitude,bs='sos') [lat first!]
        extra_pred_vars <- c(extra_pred_vars,"s(y,x,bs='sos')")
    }
    ## n.b. need to pass data to function so lme4 internals can find it ...
    ## set up formula with specified combination of random effects
    mform <- function(which_term,interax=TRUE,extra_pred_vars=NULL) {
        ## select non-NA terms
        dd <- na.omit(data.frame(rterms,forms=forms[which_term]))
        ## random-effect terms
        rterms2 <- paste("(",dd$forms,dd$rterms,")")
        ## construct fixed-effect terms
        if (!interax) {
            pp <- pars
        } else {
            pp <- paste0("(",paste(pars,collapse="+"),")^2")
        }
        pred_vars <- c(pp,rterms2)
        if (!is.null(extra_pred_vars)) {
            test_pv <- function(x) all(all.vars(parse(text=x)) %in% names(data))
            for (v in extra_pred_vars) {
                if (!test_pv(v)) stop("extra vars not found in data:",v)
            }
            pred_vars <- c(pred_vars,extra_pred_vars)
        }
        ff <- reformulate(pred_vars,response=response)
        environment(ff) <- parent.frame() ## ugh
        return(ff)
    }
    ## g4fit() is our wrapper that calls gamm4 and assigns class "gamm4"
    ##  to the result
    fitfun <- if (use_gamm4) g4fit else lmer
    ctrl <- lmerControl(optimizer=nloptwrap,
                        optCtrl=list(ftol_rel=1e-12,ftol_abs=1e-12))
    ## run just one model
    if (!is.null(single_fit)) {
        ff <- mform(single_fit,extra_pred_vars=extra_pred_vars)
        return(try(suppressWarnings(
            fitfun(ff,data=data,
                 na.action=na.exclude,
                 control=ctrl))))
    }
    dd <- expand.grid(seq_along(forms),seq_along(forms),seq_along(forms))
    ## For example, use RE #1 (intercept) for biome, RE #2 (diag) for realm, RE #3 (full) for biome $\times$ realm interaction ...
    ## mform(c(1,2,3))
    results <- list()
    for (i in seq(nrow(dd))) {
        w <- unlist(dd[i,])
        nm <- paste(rterms,"=",names(forms)[w],sep="",collapse="/")
        if (verbose) cat(i,w,nm,"\n")
        ff <- mform(w,extra_pred_vars=extra_pred_vars)
        results[[i]] <- try(suppressWarnings(
            fitfun(ff,data=data,
                   na.action=na.exclude,
                   control=ctrl)))
        names(results)[i] <- nm
    }
    return(results)
}

AICtabx <- function(fitlist) {
    force(fitlist)
    aa <- bbmle::AICtab(fitlist,mnames=names(fitlist))
    sing <- sapply(fitlist,is.singular)[attr(aa,"row.names")]
    aa2 <- data.frame(dAIC=round(aa$dAIC,1),df=aa$df,singular=sing)
    return(aa2)
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
##' load("ecoreg.RData")
##' load("allfits_restr_gamm4.RData")
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
                    xvar="NPP_log",
                    respvar=NULL,
                    auxvar="Feat_cv_sc",
                    grpvar=NULL,
                    pred_lower_lim= -3,
                    aux_quantiles=c(0.1,0.5,0.9),
                    re.form = NA,  ## exclude REs from prediction
                    alpha=0.05,
                    npts = 51,
                    lty=c(2,1,3)
                    backtrans=FALSE
                    ) {
    
    if (inherits(model,"gamm4")) {
        ## need x/y variables
        ff <- formula(model,fixed.only=TRUE,drop.smooth=FALSE)
    } else {
        ff <- formula(model,fixed.only=TRUE)
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
            pdata <- do.call(rbind,Map(function(nm,x) if (is.null(x)) NULL else data.frame(gg=nm,xx=x),
                                       names(lrange),lrange))
            names(pdata) <- c(grpvar,xvar)
        } else {
            ## xvar only
            pdata <- data.frame(seq(min(xx),max(xx),length=npts))
            names(pdata) <- xvar
        }
    }
    ## variables other than primary x-variable and aux (and maybe grpvar) are set to median
    ##  (more consistent to set to mean=0)?
    for (i in othervars) {
        if (is.null(grpvar)) {
            pdata[[i]] <- median(data[[i]])
        } else {
            pdata[[i]] <- aggregate(data[[i]],by=list(data[[grpvar]]),FUN=median)$x
        }
    }
    pdata[[mrespvar]] <- predict(model,newdata=pdata,re.form=re.form)
    if (!is.null(auxvar)) {
        fauxvar <- paste0("f",auxvar)
        pdata[[fauxvar]] <- factor(pdata[[auxvar]],labels=paste0("Q(",aux_quantiles,")"))
    }
    ## confidence intervals (fixed-effects only) on predictions
    mm <- model.matrix(formula(model),pdata)
    pvar1 <- diag(mm %*% tcrossprod(as.matrix(vcov(model)),mm))
    pdata <- transform(pdata,
                       lwr = qnorm(alpha/2,    mean=pdata[[mrespvar]],sd=sqrt(pvar1)),
                       upr = qnorm(1-alpha/2,mean=pdata[[mrespvar]],sd=sqrt(pvar1)))
    
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
    if (backtrans) {
    }
    return(pdata)
}

##' @param model fitted model
##' @param data (ecoreg)
##' @param xvar ("NPP_log"): x-variable
##' @param respvar (equal to model response by default): response variable
##' @param auxvar ("Feat_cv_sv"): auxiliary variable (e.g. for examining interactions)
##' @param ... parameters passed through to predfun
##' @examples
##' source("gamm4_utils.R")
##' load("ecoreg.RData")
##' load("allfits_restr_gamm4.RData")
##' m1 <- allfits_restr_gamm4$mbirds_log
##' plotfun(m1)
##' plotfun(m1,auxvar=NULL)
plotfun <- function(model=best_model,
                    data = ecoreg,
                    xvar="NPP_log",
                    respvar=NULL,
                    auxvar="Feat_cv_sc",
                    grpvar=NULL,
                    ylim=c(-3,1),
                    ...
                    ) {
    pdata <- predfun(model,data,xvar,respvar,auxvar,grpvar,...)
    mrespvar <- deparse(formula(model)[[2]])
    if (is.null(respvar)) respvar <- mrespvar
    gg0 <- ggplot(data,aes_(x=as.name(xvar)))
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
        gg0 <- gg0 + geom_line(data=pdata,aes_(y=as.name(mrespvar))) +
                   geom_ribbon(data=pdata,aes_(ymin=~lwr,ymax=~upr),
                        colour=NA,fill="black",alpha=0.1)
    }
    ## finish
    if (!is.null(ylim)) {
        gg0 <- gg0 +
            scale_y_continuous(limits=ylim,oob=scales::squish)
    }
    gg0 <- gg0 + theme(legend.box="horizontal")
    return(gg0 + geom_point(aes_(y=as.name(respvar),
                                colour=~biome,shape=~flor_realms)))
}

## test

if (FALSE) {
    load("ecoreg.RData")
    library(lme4)
    debug(lFormula)
    debug(fit_all)
    fit_all("plants_log")
}

pkgList <- c('lme4',      ## lmer etc.
             'gamm4',     ## spatial fits
             'bbmle',     ## AICtab, cosmetic
             'broom',     ## coef tables
             'broom.mixed', ## better coef tables
             'lattice',   ## diagnostic plots
             'gridExtra', ## arrange plots
             'ggplot2',
             'viridis',
             'plotly',
             'cowplot',
             'ggstance',  ## horizontal geoms
             'dplyr',     ## data manipulation
             'tidyr',     ## ditto
             'tibble',    ## ditto: rownames_to_column
             'remef',
             'r2glmm')

load_all_pkgs <- function() {
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
    select(-c(AIC,best,singular))
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
    data$coef %>%
    ## add summary info (singular, AIC)
    full_join(select(data$sum,c(taxon,model,singular,AIC)),
              by=c("taxon","model")) %>%
    ## pick one taxon; only fixed-effect parameters; drop intercept
    filter(taxon==focal_taxon,
           effect=="fixed",
           term!="(Intercept)") %>%
    mutate(model=shorten_modelname(model),     ## abbreviate model name
           model=reorder(model,-AIC),          ## reverse-order by AIC
           ## reverse-order by magnitude of coefficient for best model
           term=reorder(term,estimate,function(x) -x[1]))
}

## utilities for post-processing coefficient tabs
add_wald_ci <- function(data) {
    return(data %>%
           mutate(lwr=estimate-1.96*std.error,upr=estimate+1.96*std.error) %>%
           select(-c(std.error,statistic)))
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

remef_allran <- function(x,data,na.action=na.exclude) {
    require(mgcv) ## for s()
    if (inherits(x,"gamm4")) {
        ds <- glmmTMB:::drop.special2
        dh <- glmmTMB:::dropHead
        ff <- ds(formula(x$mer,random.only=TRUE),quote((1|Xr)))
        ff <- ff[-2]
        resp <- model.frame(x$mer)$y.0
    } else {
        stop("not implemented for lme4 yet")
    }
    ## gamm4 doesn't keep na.action, ugh
    na.act <- attr(model.frame(dh(formula(x$gam),quote(s)),
                     data,na.action=na.action),"na.action")
    pp <- predict(x$mer,random.only=TRUE,re.form=ff)
    rem <- napredict(na.act,resp-pp)
    return(rem)
}

##' undo the effects of scale()
##' @param x primary variable to unscale
##' @param y optionally: variable containing unscaling information
##' @examples
##' x <- 1:3
##' x <- drop(scale(1:3))
##' y <- 2:4
##' backtrans(x)
##' backtrans(y,x)
backtrans <- function(x,y=NULL) {
    
## if only x is specified, x is both the target and has unscaling info
    ## if x and y are specified, x is the target and y has unscaling info
    scvar <- if (is.null(y)) x else y
    ctr <- attr(scvar,"scaled:center")
    sc <- attr(scvar,"scaled:scale")
    if (is.null(sc) || is.null(ctr)) {
        stop("can't retrieve scaling information")
    }
    ## remove scaling information from result
    ret <- x*sc+ctr
    attr(ret,"scaled:scale") <- NULL
    attr(ret,"scaled:center") <- NULL
    return(ret)
}

## back-transform a variable by exponentiating
## *or* unscaling, depending on its name
backtrans_magic <- function(x,xname,y=NULL) {
    if (grepl("_log$",xname)) {
        return(exp(x))
    }
    if (grepl("_sc$",xname)) {
        return(backtrans(x,y))
    }
}
        
