## To do

### BMB

#### high
  
  - redo among-biome plots with Feat_log instead of NPP_log (etc.)
  - plot random effects for biome/FR models, especially with CIs ???
  - try stuff with `brms` (Bayesian) ???
  - nicer version of coefplot (with only `gamm4_2`)
  
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


	
