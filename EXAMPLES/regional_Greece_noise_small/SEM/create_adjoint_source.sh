#!/bin/bash
#
# creates adjoint source files
#
# note: run this script in directory SEM/
#############################################################

# station
station="ZKR"
network="GE"

#
# use adj_traveltime_filter.f90
#

# for single ZKR station
cp -v ../OUTPUT_FILES_0.2/ZKR*ascii ./
cd ..
ifort adj_traveltime_filter.f90
./a.out
cd SEM/

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT



exit



#
# using xcut_velocity
#

# time window
t1="8."
t2="14."

#window out single phase arrival on vertical component between t1 to t2 :
~/SPECFEM3D_GLOBE/UTILS/cut_velocity/xcut_velocity $t1 $t2 3 ../REF_SEIS/$station.$network.MX*
mv ../REF_SEIS/$station.$network.MX*adj ./

# rename adjoint source files:
rename .sem.ascii.adj .adj $station.$network.MX*adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT



#
# please add here more adjoint sources if you like...
#

