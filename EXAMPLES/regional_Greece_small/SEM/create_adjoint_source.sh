#!/bin/bash
#
# creates adjoint source files
#
# note: run this script in directory SEM/
#############################################################

# station
station="BEKI"
network="YL"

# time window
t1="130.0"
t2="147.5"

#window out single phase arrival on vertical component between t1 to t2 :
~/SPECFEM3D_GLOBE/UTILS/adjoint_sources/traveltime/xcreate_adjsrc_traveltime $t1 $t2 3 ../REF_SEIS/$station.$network.MX*
mv ../REF_SEIS/$station.$network.MX*adj ./

# rename adjoint source files:
rename .sem.ascii.adj .adj $station.$network.MX*adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT

