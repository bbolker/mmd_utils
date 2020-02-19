NCORES=4
## doesn't include args properly?
%.Rout: %.R
	Rscript --vanilla $< $@

## primary output
MixedEffects.html: ecoreg.RData allfits_sum_lme4.rds allfits_sum_gamm4.rds bestmodels_gamm4.rds allfits_restr_gamm4.rds utils.R gamm4_utils.R make.png

topfour.html: topfour.Rmd utils.R

## process 'input' data to useful version
ecoreg.rds: full_data.RData teow_data.RData proc_data.R biome_defs.csv olson_defs.csv
	Rscript --vanilla proc_data.R >proc_data.Rout

## all combinations of fits
allfits_lme4.rds: ecoreg.rds utils.R fit_batch.R
	Rscript --vanilla fit_batch.R lme4 >allfits_lme4.Rout

allfits_gamm4.rds: ecoreg.rds utils.R fit_batch.R
	Rscript --vanilla fit_batch.R gamm4 FALSE $(NCORES) >allfits_gamm4.Rout

allfits_brms.rds: ecoreg.rds utils.R fit_batch.R
	Rscript --vanilla fit_batch.R brms >allfits_brms.Rout

## biome=diag, flor_realm=diag fits
allfits_restr_gamm4.rds: ecoreg.rds utils.R fit_batch.R
	Rscript --vanilla fit_batch.R gamm4 TRUE > allfits_restr_gamm4.Rout

sumfiles = sum_batch.R utils.R ecoreg.rds allfits_gamm4.rds allfits_lme4.rds

## summaries. DRY ???
allfits_sum_lme4.rds: allfits_lme4.rds sum_batch.R
	Rscript --vanilla sum_batch.R lme4 > allfits_sum_lme4.Rout

allfits_sum_brms.rds: allfits_brms.rds sum_batch.R
	Rscript --vanilla sum_batch.R brms > allfits_sum_brms.Rout

bestmodels_gamm4.rds allfits_sum_gamm4.rds: allfits_gamm4.rds sum_batch.R
	Rscript --vanilla sum_batch.R gamm4 > allfits_sum_gamm4.Rout

## testfits.rds:         $(sumfiles) sum_batch.Rout

allprofs.rds: allfits_lme4.rds
	Rscript -e profile_batch.R >profile_batch.Rout

## requires 'dot': sudo apt-get install graphviz
## may need to edit /etc/ImageMagick-6/policy.xml
## from https://github.com/lindenb/makefile2graph
## https://unix.stackexchange.com/questions/145402/regex-alternation-or-operator-foobar-in-gnu-or-bsd-sed
make.dot: Makefile
	exec make -nd MixedEffects.html | ./make2graph | sed -E 's/red|green/gray/' >make.dot 

make.png: make.dot
	dot -Tpng make.dot >make.png

clean:
	rm -f *~ .#* .RData

## copy archive to my google drive folder
to_gd:
	make clean
	rsync -auv --exclude='.git/' --exclude='*cache*' --exclude="allfits.RData" . /media/sf_Google_Drive/fire_diversity

from_cw:
	rsync -auv collywobbles.mcmaster.ca:~/Documents/projects/mmd_fire/*.rds .

%.html: %.Rmd
	Rscript -e "rmarkdown::render(\"$<\")"

