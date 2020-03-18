##' fits all of the combinations of random effects for a particular
##' response variable
##' 
##' @param response (character) response variable
##' @param pars (character) vector of input variables
##' @param forms list of allowed random-effect terms
##' @param rterms grouping variables
##' @param interax include (2-way) interactions among fixed effects?
##' @param data data frame
##' @param single_fit fit a single specified model instead of all combinations?
##' Argument is a vector passed to mform() specifying which of forms to use
##' (1=intercept; 2=diag; 3=full)
##' @param platform which fitting package to use?
##' @param add_sos add spherical smooth to model?
##' @return what does it return?
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
                    pars=c("NPP_log_sc","Feat_log_sc","NPP_cv_sc","Feat_cv_sc"),
                    extra_pred_vars="area_km2_log_sc",
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
                    platform=c("lme4","gamm4","brms"),
                    add_sos=TRUE,
                    verbose=FALSE)
{
    platform <- match.arg(platform)
    if (add_sos && platform!="lme4") {
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
    ## g4fit() just calls gamm4 and assigns class "gamm4" to the result
    fitfun <- switch(platform,gamm4=g4fit,lme4=lmer,
                     brms=brm)
    ## set up options/arguments
    argList <- list(data=data)
    if (platform %in% c("gamm4","lme4")) {
        argList <- c(argList,
                     list(control=lmerControl(optimizer=nloptwrap,
                                              optCtrl=list(ftol_rel=1e-12,ftol_abs=1e-12),
                                              ## suppress 'singular fit' messages
                                              check.conv.singular="ignore"),
                          ## brms doesn't know about na.action ... ?
                          na.action=na.exclude))
    } else if (platform=="brms") {
        argList <- c(argList,
                     list(control=list(adapt_delta=0.99),
                          family=gaussian))
    }
    ## run just one model
    if (!is.null(single_fit)) {
        ff <- mform(single_fit,extra_pred_vars=extra_pred_vars)
        suppressWarnings(time <- system.time(res <- try(do.call(fitfun,c(list(formula=ff),argList)))))
        attr(res,"time") <- time
        return(res)
    }
    dd <- expand.grid(seq_along(forms),seq_along(forms),seq_along(forms))
    ## For example, use RE #1 (intercept) for biome, RE #2 (diag) for realm, RE #3 (full) for biome $\times$ realm interaction ...
    ## mform(c(1,2,3))
    results <- list()
    if (verbose) cat("response","index","num...","name","\n",sep=" ")
    for (i in seq(nrow(dd))) {
        w <- unlist(dd[i,])
        nm <- paste(rterms,"=",names(forms)[w],sep="",collapse="/")
        if (verbose) cat(response,i,w,nm,"\n")
        ## DRY?
        ff <- mform(w,extra_pred_vars=extra_pred_vars)
        suppressWarnings(time <- system.time(res <- try(do.call(fitfun,c(list(formula=ff),argList)))))
        attr(res,"time") <- time
        results[[i]] <- res
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
