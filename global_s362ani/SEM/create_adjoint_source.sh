#!/bin/bash
#
# creates adjoint source files
#
# note: run this script in directory SEM/
#############################################################

# station
station="PAS"
network="CI"

# time window
t1="734.5"
t2="755.0"

#window out single phase arrival on vertical component between t1 to t2 :
~/SPECFEM3D_GLOBE/UTILS/cut_velocity/xcut_velocity $t1 $t2 3 ../REF_SEIS/$station.$network.MX*
mv ../REF_SEIS/$station.$network.MX*adj ./

# rename adjoint source files:
rename .sem.ascii.adj .adj $station.$network.MX*adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT

