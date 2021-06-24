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

# trace directory
dir=$1

#############################################################

if [ "$dir" == "" ]; then echo "./usage: ./create_adjoint_source.sh dir[e.g. ../OUTPUT_FILES_2]"; exit 1; fi

# station name
sta=$network.$station

#
# use adj_traveltime_filter.f90
#

# for single ZKR station
rm -f ./$sta.*ascii
cp -v $dir/$sta.*ascii ./
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# tapers/filters station seismograms
cd ..

# compiles
$f90 adj_traveltime_filter.f90
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs adjoint source creation
./a.out
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

cd SEM/

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT
cp -v ./STATIONS_ADJOINT ../DATA/

# returns
echo
echo "done"


#
# please add here more adjoint sources if you like...
#
#
# deprecated
#
# creates traveltime adjoint source using xcut_velocity
#
#
# time window
#t1="8."
#t2="14."
#
#window out single phase arrival on vertical component between t1 to t2 :
#~/SPECFEM3D_GLOBE/UTILS/cut_velocity/xcut_velocity $t1 $t2 3 ../REF_SEIS/$sta.MX*
#~/SPECFEM3D_GLOBE/utils/adjoint_sources/traveltime/xcreate_adjsrc_traveltime $t1 $t2 3 ./$sta.MX*.sem.ascii
#
# rename adjoint source files:
#rename .sem.ascii.adj .adj $sta.MX*.adj
#
# create STATIONS_ADJOINT file with adjoint source location
#fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT
#cp -v ./STATIONS_ADJOINT ../DATA/
#
#echo
#echo "done"
#exit

