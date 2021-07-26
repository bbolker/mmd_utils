## To do

- add words/description about biome-specific plots to supplement text

- understand R^2, really? (Jaeger)
- areas etc. (Enric?)

- sort out R^2: are gamm4 KR results really sensible??? should we be including the os() smooth? compare lme4 fits? SGV fits?
    - try R^2 on lme4 models where I didn't have to hack so much.  What is "null" model ?
- check Enric's FC_FCC stuff: google drive/fire_diversity/FC_FCC_20201019.Rmd

- decide on which panels, tables, etc. to use
- compare elasticities with scaled parameters/R^2, especially for mammals NPP and fire effects, explain the change in rank clearly; graphics? numerical illustrations?

- figures
    - pretty labels in Rsq/coef plots
    - faceting (i.e. compact multi-panel plots by suppressing redundant scales/axis labels)
- fire comparisons
    - combined fire R^2 with `cmp_R2`
	- zero-fire predictions
    - no-fire models
- Feat -> Fire (fraction 
- try redrawing Figure 2



### low

- figure out amphibian intercept shift???
- R-squared: does "full model" Rsq include random effects??

- all-fire-excluded: R^2
- 
## Cleanup

- finish cleaning up pipeline
	- allfit
	   - re-add checkpointing
	   - parallelize at lower level?
	   - separate `utils.R` into model-fitting and downstream?
- divide document into backend/frontend (update table of contents?)

## Long-term/wishlist

- get everything working with brms?

