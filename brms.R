source("mmd_utils.R")
source("gamm4_utils.R")
L <- load("ecoreg.RData")

library(brms)
## from ?set_prior: default prior on SDs is half-Student with 3 df,
## scale parameter based on std dev of response (>=10)
## e.g. set_prior("<prior>", class="sd", group="<group>")

## default correlation priors lkj_(eta)
## eta=1 - flat
## eta > 1 extreme correlations less likely

library(brms)
system.time(b1 <- brm(mbirds_log ~ (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc)^2+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome)+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|flor_realms)+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome_FR),
    family=gaussian,
    data=ecoreg))
## 884 secs
pairs(b1)  ## also slow!

library(broom)
library(broom.mixed)


system.time(b2 <- brm(mbirds ~ (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc)^2+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome)+
        (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|flor_realms)+
            (NPP_log+Feat_log+NPP_cv_sc+Feat_cv_sc|biome_FR)+
       s(y,x,bs='sos'),
    family=gaussian,
    data=ecoreg))
