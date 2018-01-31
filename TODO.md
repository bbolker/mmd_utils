## To do

- BMB
    - incorporate the streamlined results into the main document
	- ?? compare GAMM4-best with lme4-best models ...
 	     - how should I calculate AIC/compare models?
		 - gamm4-stripping; be more careful removing stuff from $mer@frame
    - understand `sos` plot theta/phi?
    - other models
        - brms (`cor_car`, `cor_sar` or ?? `mgcv::s(.,bs="mrf")`: see [here](https://github.com/paul-buerkner/brms/issues/6))
		- INLA
	
- Enric: 
    - x/y in correlations: is there a way to do this with real distances? how did you last time, when we had km on the x-axis?
	- plotting correlations of gamm4 resids, correlations of lmer residuals, and correlations of log diversity on the same plot 
    - make some more beautiful maps
	
