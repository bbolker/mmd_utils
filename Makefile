%.Rout: %.R
	Rscript --vanilla $< >@

## primary output
MixedEffects.html: ecoreg.RData allfits_sum_lmer.RData allfits_sum_gamm4.RData bestmodels_gamm4.RData allfits_restr_gamm4.RData mmd_utils.R gamm4_utils.R make.png

sumfiles = mmd_reduce.R mmd_utils.R ecoreg.RData allfits.RData allfits_lmer.RData allfits_brms.RData

## DRY ???
allfits_sum_lmer.RData: $(sumfiles) mmd_reduce.Rout

allfits_sum_brms.RData: $(sumfiles) mmd_reduce.Rout

bestmodels_gamm4.RData: $(sumfiles) mmd_reduce.Rout
testfits.RData:         $(sumfiles) mmd_reduce.Rout

## process 'input' data to useful version
ecoreg.RData: full_data.RData teow_data.RData mmd_procdata.R biome_defs.csv olson_defs.csv
	$(RBATCH) mmd_procdata.R

## all combinations of fits
allfits.RData: ecoreg.RData mmd_utils.R mmd_fitbatch.R
	Rscript --vanilla mmd_fitbatch.R

## biome=diag, flor_realm=diag fits
allfits_restr.RData: ecoreg.RData mmd_utils.R mmd_fitbatch.R
	$(RBATCH) mmd_fitbatch2.R

allprofs.RData: allfits.RData
	$(RBATCH) mmd_profilebatch.R

## requires 'dot': sudo apt-get install graphviz
## may need to edit /etc/ImageMagick-6/policy.xml
## from https://github.com/lindenb/makefile2graph
make.png: Makefile
	. ./genmakegraph

clean:
	rm -f *~ .#* .RData

## copy archive to my google drive folder
to_gd:
	make clean
	rsync -auv --exclude='.git/' --exclude='*cache*' --exclude="allfits.RData" . /media/sf_Google_Drive/fire_diversity

%.html: %.Rmd
	Rscript -e "rmarkdown::render(\"$<\")"

