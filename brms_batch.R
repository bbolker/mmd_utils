source("mmd_utils.R")
source("gamm4_utils.R")
L <- load("ecoreg.RData")

## use mmd_utils model-construction code?

library(brms)
system.time(b1 <- brm(mbirds_log ~ (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc)^2+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome)+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|flor_realms)+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome_FR),
    family=gaussian,
    data=ecoreg))

system.time(b2 <- brm(mbirds_log ~ (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc)^2+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome)+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|flor_realms)+
            (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome_FR)+
       s(y,x,bs='sos'),
    family=gaussian,
    data=ecoreg))

save("b1","b2",file="brms_batch.RData")
