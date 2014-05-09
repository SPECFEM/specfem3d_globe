#!/bin/csh

### Dimitri Komatitsch, CNRS Marseille, France, May 2014

grd2xyz -Rg etopo2v2c.grd | awk '{ print $3 }' > topo_bathy_etopo2v2c_original_unmodified_unsmoothed.dat

