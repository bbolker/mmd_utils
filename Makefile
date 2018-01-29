%.html: %.Rmd
	Rscript -e "rmarkdown::render(\"$<\")"

MixedEffects.html: gamm4_utils.R mmd_utils.R ecoreg.RData allfits.RData allfits_restr.RData allprofs.RData MixedEffects.Rmd

## Don't do it, will overwrite 'proper' ecoreg
## ecoreg.RData: ecoreg_means.RData mmd_procdata.R biome_defs.csv olson_flor.csv
##	R CMD BATCH mmd_procdata.R

## all combinations of fits
allfits.RData: ecoreg.RData mmd_utils.R mmd_fitbatch.R
	R CMD BATCH mmd_fitbatch.R

## biome=diag, flor_realm=diag fits
allfits_restr.RData: ecoreg.RData mmd_utils.R mmd_fitbatch.R
	R CMD BATCH mmd_fitbatch2.R

allprofs.RData: allfits.RData
	R CMD BATCH mmd_profilebatch.R

## requires 'dot': sudo apt-get install graphviz
make.png: Makefile
	~/bin/genmakegraph

clean:
	rm -f *~ .#* .RData

## copy archive to my google drive folder
to_gd:
	make clean
	rsync -auv --exclude='.git/' --exclude='*cache*' --exclude="allfits.RData" . /media/sf_Google_Drive/Enric

