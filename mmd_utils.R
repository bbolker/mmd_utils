## basic check for singular random effects (i.e. overfitted)
is.singular <- function(fit,tol=1e-4) {
    any(abs(getME(fit,"theta"))<tol)
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
                    pars=c("NPP_log","Feat_log","NPP_cv_ctr","Feat_cv_ctr"),
                    ## possible random-effect models
                    forms=c(int="1|",  ## intercept-only
                            diag=paste("1+",paste(pars,collapse="+"),"||"),
                            full=paste("1+",paste(pars,collapse="+"),"|")
                            ),
                    ## grouping variables
                    rterms=c("biome","flor_realms","biome:flor_realms"),
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
                            control=lmerControl(optimizer=nloptwrap,
                            na.action=na.exclude,
                            optCtrl=list(ftol_rel=1e-12,ftol_abs=1e-12)))))
        names(results)[i] <- nm
    }
    return(results)
}

AICtabx <- function(fitlist) {
    aa <- bbmle::AICtab(results)
    sing <- sapply(results,is.singular)[attr(aa,"row.names")]
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

## test

if (FALSE) {
    load("ecoreg2.RData")
    library(lme4)
    debug(lFormula)
    debug(fit_all)
    fit_all("plants_log")
}

