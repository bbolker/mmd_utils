## To do

- BMB
    - unbreak stuff with 'effects' (switch to broom.mixed)
    - hack correct predictions from gamm4 objects
	     - need to combine predictions (fixed+smooth only), predictions (random only)
    - understand `sos` plot theta/phi?
    - other models
        - brms (`cor_car`, `cor_sar` or ?? `mgcv::s(.,bs="mrf")`: see [here](https://github.com/paul-buerkner/brms/issues/6))
		- INLA
    - re-do fits with `gamm4`

- Enric: look up Area in metadata
    - centroids; weighted centroid of polygons?
	ecoreg_xy ...

- parallelize `mmd_fitbatch.R` ?
