g4fit <- function(...) {
    g <- gamm4(...)
    class(g) <- "gamm4"
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

