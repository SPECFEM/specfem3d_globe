#!/bin/bash
#
# creates adjoint source files
#
# note: run this script in directory SEM/
#############################################################

dir=$1
if [ "$dir" == "" ]; then echo "usage: ./create_adjoint_sources.sh OUTPUT_FILES/"; exit 1; fi

# station
sta="BEKI"
net="YL"

name=${net}.${sta}

# time window
t1="130.0"
t2="147.5"

# setup SEM/
mkdir -p SEM/
rm -rf SEM/*.adj SEM/$name.MX*
cp -v $dir/$name.MX* SEM/

echo
echo

cd SEM/

#window out single phase arrival on vertical component between t1 to t2 :
../../../utils/adjoint_sources/traveltime/xcreate_adjsrc_traveltime $t1 $t2 3 ./$name.MX*
if [[ $? -ne 0 ]]; then echo "xcreate_adjsrc_traveltime failed."; exit 1; fi

# rename adjoint source files:
#rename .sem.ascii.adj .adj $sta.MX*adj
rename -v 's/\.sem.ascii.adj/\.adj/' *.sem.ascii.adj
if [ ! -e $tmp ]; then echo "renaming failed."; exit 1; fi

# clean up
rm -f *.sem.ascii

# create STATIONS_ADJOINT file with adjoint source location
fgrep $sta ../DATA/STATIONS > ./STATIONS_ADJOINT

# second station as fake station
# (used for checking output adjoint wavefield traces)
sta="BGIO"
net="SR"

# checks if BEKI adjoint trace exists
tmp=$name.MXZ.adj
if [ ! -e $tmp ]; then echo "no trace $name.MXZ.adj, exiting."; exit 1; fi

# zero out adjoint trace for second station
awk '{print $1,"0.0"}' $tmp > ${net}.${sta}.MXE.adj
awk '{print $1,"0.0"}' $tmp > ${net}.${sta}.MXN.adj
awk '{print $1,"0.0"}' $tmp > ${net}.${sta}.MXZ.adj

fgrep $sta ../DATA/STATIONS >> ./STATIONS_ADJOINT

# copy setup
cp -v STATIONS_ADJOINT ../DATA/


echo
echo "done"
echo

