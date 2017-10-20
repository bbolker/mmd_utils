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
fit_all <- function(response="mbirds_log",
                    pars=c("NPP_log","Feat_log","NPP_cv_sc","Feat_cv_sc"),
                    ## possible random-effect models
                    forms=c(int="1|",  ## intercept-only
                            diag=paste("1+",paste(pars,collapse="+"),"||"),
                            full=paste("1+",paste(pars,collapse="+"),"|")
                            ),
                    ## grouping variables
                    rterms=c("biome","flor_realms","biome_FR"),
                    interax=TRUE,
                    data=ecoreg) {
    ## n.b. need to pass data to function so lme4 internals can find it ...
    ## set up formula with specified combination of random effects
    mform <- function(which_term,interax=TRUE) {
        rterms <- paste("(",forms[which_term],rterms,")")
        if (!interax) {
            pp <- pars
        } else {
            pp <- paste0("(",paste(pars,collapse="+"),")^2")
        }
        ff <- reformulate(c(pp,rterms),response=response)
        environment(ff) <- parent.frame() ## ugh
        return(ff)
    }
    dd <- expand.grid(1:3,1:3,1:3)
    ## For example, use RE #1 (intercept) for biome, RE #2 (diag) for realm, RE #3 (full) for biome $\times$ realm interaction ...
    ## mform(c(1,2,3))
    results <- list()
    for (i in seq(nrow(dd))) {
        w <- unlist(dd[i,])
        nm <- paste(rterms,"=",names(forms)[w],sep="",collapse="/")
        ## cat(i,w,nm,"\n")
        ff <- mform(w)
        results[[i]] <- try(suppressWarnings(lmer(ff,data=data,
                            na.action=na.exclude,
                            control=lmerControl(optimizer=nloptwrap,
                            optCtrl=list(ftol_rel=1e-12,ftol_abs=1e-12)))))
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

get_best_name <- function(fitlist,allow_sing=FALSE) {
    aa <- AICtabx(fitlist)
    if (!allow_sing) {
        aa <- aa[!aa$singular,]
    }
    return(rownames(aa)[which.min(aa$dAIC)])
}

plotfun <- function(model=best_model,
                    xvar="NPP_log",
                    auxvar="Feat_cv_sc",
                    respvar=NULL,
                    aux_quantiles=c(0.1,0.5,0.9),
                    pred_lower_lim= -3,
                    data = ecoreg,
                    re.form = NA  ## exclude REs from prediction
                    ) {
    ## get LHS of formula
    mrespvar <- deparse(formula(best_model)[[2]])
    if (is.null(respvar)) respvar <- mrespvar
    ## extract variables from data set by name
    xx <- data[[xvar]]
    aa <- drop(data[[auxvar]])
    ## add factor equiv of auxvar to data
    av <- all.vars(formula(model,fixed.only=TRUE))
    othervars <- setdiff(av,c(respvar,xvar,auxvar,mrespvar))
    ## construct prediction frame
    pdata <- expand.grid(seq(min(xx),max(xx),length=51),
                         quantile(aa,aux_quantiles))
    names(pdata) <- c(xvar,auxvar)
    ## variables other than primary x-variable and aux are set to median
    ##  (more consistent to set to mean=0)?
    for (i in othervars) {
        pdata[[i]] <- median(data[[i]])
    }
    fauxvar <- paste0("f",auxvar)
    pdata[[fauxvar]] <- factor(pdata[[auxvar]],labels=paste0("(",aux_quantiles,")"))
    pdata[[focal]] <- predict(best_model,newdata=pdata,re.form=re.form)
    pdata[[focal]][pdata[[focal]] < pred_lower_lim] <- NA
    gg0 <- ggplot(data,aes_(x=as.name(xvar),y=as.name(respvar)))+
        ## geom_encircle(aes(group=biome_FR),expand=0)+  ## ugly ...
        ## use model respvar for predicted values
        geom_line(data=pdata,aes_(y=as.name(mrespvar),linetype=as.name(fauxvar))) +
        ## FIXME (don't hard-code line types)
        scale_linetype_manual(values=c(2,1,3))+
        ## scale_y_continuous(limits=c(-3,1),oob=scales::squish)+
        theme(legend.box="horizontal")
    return(gg0 + geom_point(aes(colour=biome,shape=flor_realms)))
}

## test

if (FALSE) {
    load("ecoreg2.RData")
    library(lme4)
    debug(lFormula)
    debug(fit_all)
    fit_all("plants_log")
}

pkgList <- c('lme4',      ## lmer etc.
             'bbmle',     ## AICtab, cosmetic
             'broom',     ## coef tables
             'lattice',   ## diagnostic plots
             'gridExtra', ## arrange plots
             'ggplot2',
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



    
    
