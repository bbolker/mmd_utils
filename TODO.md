## To do

- Enric: look up Area in metadata
    - centroids; weighted centroid of polygons?
	ecoreg_xy ...
- BMB:
    - use coefficients rather than conditional modes/deviations
    - spatial autocorrelation stuff:
        - brms: either use `cor_car`, `cor_sar` or ?? `mgcv::s(.,bs="mrf")` (see [here](https://github.com/paul-buerkner/brms/issues/6))
	    - `gamm4` ?
- parallelize `mmd_fitbatch.R` ?
