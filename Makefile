%.html: %.Rmd
	Rscript -e "rmarkdown::render(\"$<\")"

## primary output
MixedEffects.html: gamm4_utils.R mmd_utils.R ecoreg.RData allfits.RData allfits_restr.RData allprofs.RData MixedEffects.Rmd make.png

## process 'input' data to useful version
ecoreg.RData: full_data.RData teow_data.RData mmd_procdata.R biome_defs.csv olson_defs.csv
	R CMD BATCH mmd_procdata.R

## all combinations of fits
allfits.RData: ecoreg.RData mmd_utils.R mmd_fitbatch.R
	R CMD BATCH mmd_fitbatch.R

## biome=diag, flor_realm=diag fits
allfits_restr.RData: ecoreg.RData mmd_utils.R mmd_fitbatch.R
	R CMD BATCH mmd_fitbatch2.R

allprofs.RData: allfits.RData
	R CMD BATCH mmd_profilebatch.R

## requires 'dot': sudo apt-get install graphviz
## may need to edit /etc/ImageMagick-6/policy.xml
make.png: Makefile
	. ./genmakegraph

clean:
	rm -f *~ .#* .RData

## copy archive to my google drive folder
to_gd:
	make clean
	rsync -auv --exclude='.git/' --exclude='*cache*' --exclude="allfits.RData" . /media/sf_Google_Drive/fire_diversity


