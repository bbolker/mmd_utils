load("ecoreg_means.RData")
## labels: biome.code.name, short.biome.names, flor.code, zoor.code
## "bindex" ?
## "db.land", "db.islands"
## "ecoreg"
## "fd" ?

## read area definitions: this replicates some information
## found in ecoreg_means.RData ...
read_fun <- function(x) read.csv(x,quote="'",stringsAsFactors=FALSE,
                                 strip.white=TRUE,
                                 ## treat 'NA' in abbrevs as a real value!
                                 na.strings="")

biome_defs <- read_fun("biome_defs.csv")
flor_defs <- read_fun("olson_flor.csv")

## predictors (must be non-zero)
predvars <- c("NPP_mean","NPP_cv_inter","Feat_mean","Feat_cv_inter")
nm <- names(ecoreg)
## potential grouping variables
grpvars <- c("biome",grep("_(realms|regions)$",nm,value=TRUE))
## potential response variables:
##  plants, plus everything starting with m (amph, birds, mamm)
respvars <- c("plants",nm[grepl("^m",nm) & !grepl("regions$",nm)])
ecoreg <- ecoreg[c(respvars,predvars,grpvars)]

## select non-zero vals
for (v in predvars) {
    ecoreg <- ecoreg[ecoreg[[v]]>0,]
}
## eliminate rock & ice, lakes
ecoreg <- ecoreg[ecoreg$biome<98,]
## convert to factors
for (v in grpvars) {
    ecoreg[[v]] <- factor(ecoreg[[v]])
}


## set factor names
ecoreg$biome <- factor(ecoreg$biome,
                       levels=seq(nrow(biome_defs)),
                       labels=biome_defs$abbrev)
ecoreg$flor_realms <- factor(ecoreg$flor_realms,
                       levels=seq(nrow(flor_defs)),
                       labels=flor_defs$name)
## better ordering?
ecoreg$biome_FR <- factor(paste(flor_defs$abbrev[ecoreg$flor_realms],
                            ecoreg$biome,sep=":"))

## log-scale all non-CV predictors
log_vars <- c(respvars,grep("_cv_inter",predvars,invert=TRUE,value=TRUE))
for (v in log_vars) {
    scv <- gsub("(_inter|_mean)","",paste0(v,"_log"))
    ecoreg[[scv]] <- log(ecoreg[[v]]/mean(ecoreg[[v]],na.rm=TRUE))
}

## center CVs
ctr_vars <- grep("_cv_inter",predvars,value=TRUE)      
for (v in ctr_vars) {
    scv <- gsub("(_inter|_mean)","",paste0(v,"_ctr"))
    ## drop() so we can use tidyverse later if we want ...
    ecoreg[[scv]] <- drop(scale(ecoreg[[v]],scale=FALSE,center=TRUE))
}   
save("ecoreg","biome_defs","flor_defs","predvars","respvars",
     file="ecoreg2.RData")
