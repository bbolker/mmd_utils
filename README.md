# fire/diversity analyses with mixed models

## model summaries

`allfits_sum_(lme4|gamm4).rds` are files containing
lists (based on `allfits_(lme4|gamm4)_res.rds`)). The difference between these two model sets is that `gamm4_res` models use a spherical spline to model spatial autocorrelation.


with the following components:

- `coefs`: coefficients
- `sum`: model summaries
- `pred`: fitted values and residuals

All three components are data frames. All share the columns

- `taxon`: "plants_log", "mamph_log", "mamm_log", "mbirds_log"
- `model`: e.g. "biome=full/flor_realms=full/biome_FR=full"

- `coefs` has fixed and random-effect parameters, with standard errors (no profile CIs)
- `pred` has original variables from `ecoreg` used in the model, plus `.fitted`, `.resid`, `.fixed` (fixed-effect-only predictions)
- `sum` has AIC, `singular` (logical: is model fit singular?), `df` (number of model parameters), `best` (logical: is model fit the best non-singular fit for the given taxon?)

## file structure

- `MixedEffects.Rmd` is the primary document
- there is a Makefile that describes the dependencies/build rules
