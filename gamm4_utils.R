
## copied from glmmTMB:
reOnly <- function(f,response=FALSE,bracket=TRUE) {
    ff <- f
    if (bracket)
        ff <- lapply(findbars(ff),makeOp,quote(`(`)) ## bracket-protect terms
    ff <- Reduce(function(x,y) makeOp(x,y,op=quote(`+`)),ff)
    if (response && length(f)==3) {
        form <- makeOp(f[[2]],ff,quote(`~`))
    } else {
        form <- makeOp(ff,quote(`~`))
    }
    return(form)
}

## combine unary or binary operator + arguments (sugar for 'substitute')
## FIXME: would be nice to have multiple dispatch, so
## (arg,op) gave unary, (arg,arg,op) gave binary operator
makeOp <- function(x,y,op=NULL) {
    if (is.null(op)) {  ## unary
        substitute(OP(X),list(X=x,OP=y))
    } else substitute(OP(X,Y), list(X=x,OP=op,Y=y))
}

g4fit <- function(formula,data,na.action,...) {
    ## remove all NA-containing rows
    mf <- data[all.vars(formula)]
    na.removed <- attr(na.action(mf),"na.action")
    data <- data[complete.cases(mf),]
    ## now fit gamm4
    g <- gamm4(formula=nobars(formula),
               random=as.formula(reOnly(formula)),
               data=data,...)
    ## now restore information
    attr(g$mer@frame,"na.action") <- na.removed
    class(g) <- c("gamm4","list")
    return(g)
}

logLik.gamm4 <- function(object, ...) {
    logLik(object$mer, ...)
}

tidy.gamm4 <- function(x,...) {
    r <- tidy(x$mer,...)
    r$term <- gsub("^X","",r$term)
    return(r)
}

## ugh. result of trial-and-error with object_size(), neutralizing
##  environments where we can find them.  Wonder if extractors still
##  work?
strip_gamm4_env <- function(x) {
    gam_env_obj <- c("pred.formula","formula","pterms","terms","model")
    for (i in gam_env_obj) {
        environment(x$gam[[i]]) <- NULL
    }
    x$gam$model <- c(x$gam$model)
    mer_env_obj <- c("pp","frame")
    for (i in mer_env_obj) {
        environment(slot(x$mer,i)) <- NULL
    }
    ## should strip environment from formula; don't trash formula ...
    x$mer@frame <- as.data.frame(c(x$mer@frame))
    return(x)
}
