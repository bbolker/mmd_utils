## To do

### BMB

#### high

- check `Enric/Figures_Ms20180523.html`
- figure 1: screwed up?
bird/feat_cv_sc slope
mammals top 4 effects


- ggpairs?
- stop hard-coding y limits?
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


	
