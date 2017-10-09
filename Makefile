%.html: %.Rmd
	Rscript -e "rmarkdown::render(\"$<\")"

MixedEffects.html: ecoreg2.RData mmd_utils.R allfits.RData MixedEffects.Rmd	

ecoreg2.RData: ecoreg_means.RData mmd_procdata.R biome_defs.csv olson_flor.csv
	R CMD BATCH mmd_procdata.R

allfits.RData: ecoreg2.RData mmd_utils.R mmd_fitbatch.R
	R CMD BATCH mmd_fitbatch.R

## requires 'dot': sudo apt-get install graphviz
make.png: Makefile
	~/bin/genmakegraph


