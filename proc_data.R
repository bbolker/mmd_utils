## Enric Batllori, 2018/05/04

## simplified version of "pipeline_full_Enric.R" that starts from a matrix (full_data) with the RAW data for each pixel
## matrix with 343429 rows that was generated/saved in the first part of "pipeline_full_Enric.R"
## 1. generate the initial 'ecoreg' dataset (including biome, ecoreg, biodiversity, NPP, Feat)
## 2. generate the working datset 'ecoreg' (log and scale transformation of original variables)
## 3. incorporate the x-y coordinates of each ecoregion

## FIXME: cull needed packages?
library(raster) ## for modal()
library(data.table)

## --------------------------------------------------------------------------------
#   compute ECOREGION MEANS -- working data set
### -------------------------------------------------------------------------------------

load("full_data.RData")
full_data <- data.table(full_data)

## !!! NOTE !!! only data from land masses is retained! 
db_land = full_data[land == 1,]
## sorting the data set by ecoregion id ('eco_id')
setkey(db_land,eco_id)
db_land = db_land[!biomes %in% c('98','99'),] ## remove rock, ice, water

## 1. ecoregion-level mean or model values
ecoreg = db_land[,list(Ncells = .N,	## total number of cells per ecoregion
                       Area = sum(AREA_cell),       # total AREA of each ecoregion
                       biome = as.numeric(modal(biomes,na.rm=TRUE, ties ='random')),
                       ## modal value for Biome (Olson et al. 2001)
                       flor_realms = as.numeric(modal(FlorRealms, na.rm=TRUE, ties ='random')),
                       ## modal value for biogeographical realms (Olson et al. 2001)
                       mamph = mean(amph,na.rm=TRUE),   ## mean AMPHIBIAN diversity - GLOBAL
                       mbirds = mean(birds,na.rm=TRUE), ## mean BIRDS diversity - GLOBAL
                       mmamm = mean(mamm,na.rm=TRUE),  	## mean MAMMALS diversity - GLOBAL                                
                       NPP_mean = mean(NPP_mean,na.rm=TRUE),	## (inter-annual) mean of yearly sum of NPP
                       NPP_cv_inter = mean(NPP_cv_inter,na.rm=TRUE),	## inter-annual)CV of yearly sum of NPP
                       Feat_mean = mean(Feat_mean,na.rm=TRUE),		## (inter-annual) mean of FRACTION of NPP eaten by fire (Emissions/NPP)
                       Feat_cv_inter = mean(Feat_cv_inter,na.rm=TRUE),	## inter-annual CV of FRACTION of NPP eaten by fire
                       area_k = mean(area,na.rm=TRUE),			## Kier et al. 2005 data - cell area
                       plants = mean(sp_wfig,na.rm=TRUE)),		## Kier et al. 2005 data - plant richness
                 by = eco_id]

## convert the data.table to a data.frame
ecoreg = as.data.frame(ecoreg)

## read realm/biome definitions
read_fun <- function(x) {
    read.csv(x,quote="'",stringsAsFactors=FALSE,
                                 strip.white=TRUE,
                                 comment="#",
                                 ## treat 'NA' in abbrevs as a real value!
                                 na.strings="")
}

biome_defs <- read_fun("biome_defs.csv")
flor_defs <- read_fun("olson_defs.csv")

## predictors (must be non-zero)
predvars <- c("NPP_mean","NPP_cv_inter","Feat_mean","Feat_cv_inter")
nm <- names(ecoreg)
## potential grouping variables
grpvars <- c("biome",grep("_(realms|regions)$",nm,value=TRUE))
## potential response variables:
##  plants, plus everything starting with m (amph, birds, mamm)
respvars <- c("plants",nm[grepl("^m",nm) & !grepl("regions$",nm)])
## Adding three additional columns that will be needed to assign the residuals to each ecoregion
ecoreg <- ecoreg[c('eco_id','Ncells','Area',respvars,predvars,grpvars)]

## select non-zero vals
for (v in predvars) {
    ecoreg <- ecoreg[ecoreg[[v]]>0,]
}
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

## log/scale/center variables

## response variables: log; don't scale (but do include mean/scale attr)
for (v in respvars) {
    scv <- paste0(v,"_log")
    tmpvar <- log(ecoreg[[v]])
    ecoreg[[scv]] <- tmpvar
    ## scale with both FALSE attaches no attributes anyway
    ## drop(scale(ecoreg[[v]],scale=FALSE, center=FALSE))
}

## non-CV predictors; log and center/scale
log_sc_vars <- grep("_cv_inter",predvars,invert=TRUE,value=TRUE)
for (v in log_sc_vars) {
    scv <- gsub("(_inter|_mean)","",paste0(v,"_log_sc"))
    tmpvar <- log(ecoreg[[v]])
    ecoreg[[scv]] <- drop(scale(tmpvar,scale=TRUE,center=TRUE))
}


## center *and* scale CVs
ctr_vars <- grep("_cv_inter",predvars,value=TRUE)      
for (v in ctr_vars) {
    scv <- gsub("(_inter|_mean)","",paste0(v,"_sc"))
    ## drop() so we can use tidyverse later if we want ...
    ecoreg[[scv]] <- drop(scale(ecoreg[[v]],scale=TRUE,center=TRUE))
}

## COMPUTING "centroids" of each Ecoregion -- XY coordinates of the largest polygon in each ecoregion (from the original wwf_terr_ecos shapefile)
# !! this needs to be done after mmd_procdata.R or the x-y coordinates are removed... 

if (FALSE) {
    ## Original code to compute centroids for each polygon in the original wwf_terr_ecos shapefile
  ## load ecoregions shapefile
    teow <- readOGR('...','wwf_terr_ecos', verbose = FALSE)
    ## computing centroids (does not differ significantly from gCentroids function from rgeos)
    xy_teow <- coordinates(teow)
    ## generating a matrix with original teow polygon-level data (teow@data) + x,y coordinates
    teow_data <- cbind(teow@data,x=xy_teow[,1],y=xy_teow[,2])
    ## data output -- this data set will be loaded below to add x,y to 'ecoreg'
    save(teow_data,file='.../teow_data.RData')
}

load("teow_data.RData")
## selecting only ecoregions within the dataset within the shapefile
## keep only eco_id, area_km2, x, y 
xy_area_teow <- data.table(subset(teow_data,ECO_ID %in% ecoreg$eco_id,
                                  select=c(ECO_ID,area_km2,Shape_Area,x,y)))
setnames(xy_area_teow,c('eco_id','area_km2','Shape_Area','x','y'))
## selecting bigger 'Shape_Area' by ecoregion -- using data.table
## rationale to select unique rows from a data table 
## in data.table .I represents the row number in the original data set
biggest <- xy_area_teow[xy_area_teow[, .I[Shape_Area == max(Shape_Area)], by=c('eco_id')]$V1]

chksc <- function(x) {
    all(c("scaled:center","scaled:scale") %in% names(attributes(x)))
}
## check that attributes are preserved
stopifnot(chksc(ecoreg$NPP_log_sc))

## https://stackoverflow.com/questions/20306853/maintain-attributes-of-data-frame-columns-after-merge
## issues with dplyr::join() and merge() dropping attributes,
## but OK with merge.data.table() ...
ecoreg <- merge(data.table(ecoreg),data.table(biggest),by="eco_id")

ecoreg <- as.data.frame(ecoreg)
## check that attributes are preserved
stopifnot(chksc(ecoreg$NPP_log))

saveRDS(ecoreg, file="ecoreg.rds")









