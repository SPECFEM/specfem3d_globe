#!/bin/bash
#
# creates adjoint source files
#
# note: run this script in directory SEM/
#############################################################

# station
station="GAR"
network="II"

sta=$network.$station

# time window
t1="190."
t2="225."

#window out single phase arrival on vertical component between t1 to t2 :
~/SPECFEM3D_GLOBE/utils/adjoint_sources/traveltime/xcreate_adjsrc_traveltime $t1 $t2 3 ../REF_SEIS/$sta.MX*
mv ../REF_SEIS/$sta.MX*adj ./

# rename adjoint source files:
rename .sem.ascii.adj .adj $sta.MX*adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT

