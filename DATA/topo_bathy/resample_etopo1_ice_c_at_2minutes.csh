#!/bin/csh

### Dimitri Komatitsch, CNRS Marseille, France, May 2014

# =ns means that we store the output in integer format
grdsample etopo1_ice_c.grd -Getopo1_ice_c_resampled_at_2minutes.grd=ns -Rd -I2m

