## To do

### new

- divide document into backend/frontend (update table of contents?)
- finish cleaning up pipeline
    - test remef/backtrans/etc.
	- allfit
	   - re-add checkpointing
	   - parallelize at lower level?
	   - separate `utils.R` into model-fitting and downstream?
- prettify 
    - switch to scale_colour_discrete_qualitative
- topfour
	- drop some y-axes? facet?
- make "transformed" plots nicer? (what is the point of "transformed"? should we do one thing on non-log scales?)
- get everything working with brms?

#### high

- restricted fits (biome-by-biome estimates): order the effects so we can look at them
- add coefficient plot to "topfour" document
- can we explain the unconditional vs conditional stuff graphically?
      - what's the best example? feat_cv for amphibians? mammals?

    
```
library(corrplot)
library(tidyverse)
e2 <- (ecoreg
    %>% select(mmamm_log, NPP_log, Feat_log, NPP_cv_sc, Feat_cv_sc)
)
cc <- function(x) {
    corrplot.mixed(x,lower="ellipse",upper="number")
}
## correlations in the raw data
cc(cor(e2))
cc(cov2cor(as.matrix(vcov(best_models$mmamm_log$mer)))[-1,-1])

library(rgl)
with(e2,plot3d(Feat_cv_sc, NPP_log, mmamm_log))
library(gg3D)
## do something here?
```
- interactions as panel plots; for each panel, show
    - points in the range
    - predicted line for *overall* median aux var
	- predicted line for median of aux var in the plot (e.g. 0.2 quantile
	  for the low-range plot)

- ggpairs?
- effect of fixed_keep="NPP_log" on rem1 with mammals
     - make sure to remove smooth term involved?
- conditional coefficient plots with uncertainties
- plots with actual values: unscaling/unlogging
  - test
  - remove fixed effects as well as random effects in remef plots? (done, not tested)
  
#### medium

  - update Makefile/workflow description
	- back-transform?
	- check singularity calcs
	- re-do partial-effects plots with a single data frame/faceted? ggplotly?
	- quadratic??
	
#### old/low
	
	- ?? compare GAMM4-best with lme4-best models ...
 	     - how should I calculate AIC/compare models?
		 - gamm4-stripping; be more careful removing stuff from $mer@frame
    - understand `sos` plot theta/phi?
    - other models
        - brms (`cor_car`, `cor_sar` or ?? `mgcv::s(.,bs="mrf")`: see [here](https://github.com/paul-buerkner/brms/issues/6))
		- INLA

### Enric: 
   - x/y in correlations: right now it looks like these are being computed with lat/long; if so then the horizontal axis is a bit weird. Is there a way to do this with real distances? how did you do it the last time, when we had km on the x-axis?
   - it would be nice, for each taxon, to plot correlations of gamm4 resids, correlations of lmer residuals, and correlations of log diversity on the same plot (easier comparison, even if the scales are different for each curve -- it will emphasize how much autocorrelation decreases from log-diversity -> lmer resids -> gamm4 resids


	
