source("utils.R")
source("gamm4_utils.R")
L <- load("ecoreg.RData")
library(brms)
options(mc.cores=3)

## use mmd_utils model-construction code?

form1 <- mbirds_log ~ (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc)^2+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome)+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|flor_realms)+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome_FR)

t1 <- system.time(b1 <- brm(form1,
    control=list(adapt_delta=0.99),
   family=gaussian,
   data=ecoreg))

save("b1",file="brms_batch.RData")

## second half
load("brms_batch.RData")
form2 <- update(form1, . ~ . + s(y,x,bs="sos"))

t2 <- system.time(b2 <- brm(form2,
    control=list(adapt_delta=0.99),
    family=gaussian,
    data=ecoreg))

save("t1","t2","b1","b2",file="brms_batch.RData")

