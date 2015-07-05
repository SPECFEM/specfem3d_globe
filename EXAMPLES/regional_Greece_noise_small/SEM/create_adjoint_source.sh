#!/bin/bash
#
# creates adjoint source files
#
# note: run this script in directory SEM/
#############################################################

# station
station="ZKR"
network="GE"

# fortran compiler
f90=gfortran  #ifort
#############################################################

sta=$network.$station

#
# use adj_traveltime_filter.f90
#

# for single ZKR station
rm -f ./$sta.*ascii
cp -v ../OUTPUT_FILES_0.2/$sta.*ascii ./

# tapers/filters station seismograms
cd ..
$f90 adj_traveltime_filter.f90
./a.out
cd SEM/

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT

# returns
echo
echo "done"
exit

#
# please add here more adjoint sources if you like...
#

# deprecated
#
# creates traveltime adjoint source using xcut_velocity
#

# time window
t1="8."
t2="14."

#window out single phase arrival on vertical component between t1 to t2 :
#~/SPECFEM3D_GLOBE/UTILS/cut_velocity/xcut_velocity $t1 $t2 3 ../REF_SEIS/$sta.MX*
~/SPECFEM3D_GLOBE/utils/adjoint_sources/traveltime/xcreate_adjsrc_traveltime $t1 $t2 3 ./$sta.MX*.sem.ascii

# rename adjoint source files:
rename .sem.ascii.adj .adj $sta.MX*.adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT

echo
echo "done"
exit
