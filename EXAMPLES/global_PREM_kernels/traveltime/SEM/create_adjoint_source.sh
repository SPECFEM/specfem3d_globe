#!/bin/bash
#
# creates traveltime adjoint source files
#
# notes: 1- run this script in directory SEM/
#        2- make sure you have SAC to use proces_syn.pl script
#        3- make sure you have compiled create_adjsrc_traveltime.f90
#          under ~/SPECFEM3D_GLOBE/UTILS/adjoint_sources/traveltime
#############################################################

# station
station="ANMO"
network="IU"

sta=$network.$station

# time window
t1=500.0
t2=600.0

# corner periods
tmin=90.0
tmax=400.0

echo "adjoint sources:"
echo "  station: $sta"
echo "  filter tmin/tmax = $tmin / $tmax"
echo ""

# copies over seismograms
rm -f ./$sta.MX*
cp -v ../OUTPUT_FILES/$sta.MX* ./

# filter synthetics within a bandpass filter with corner periods of f1 and f2
echo
echo "filtering synthetics"
echo
~/SPECFEM3D_GLOBE/utils/seis_process/process_syn.pl -t $tmin/$tmax -P 6/2 -x filt ./$sta.MX*

if [ $? -ne 0 ]; then echo "filtering the synthetics failed, please check your script and setup..."; exit 1; fi

rm -rf ./$sta.MX*sac

# convert sac files to ascii
sac2asc ./$sta.MXE*filt > ./$sta.MXE.filt.ascii
sac2asc ./$sta.MXN*filt > ./$sta.MXN.filt.ascii
sac2asc ./$sta.MXZ*filt > ./$sta.MXZ.filt.ascii

rm -rf ./$sta.MX*sac.filt

# window out single phase arrival on vertical component between t1 to t2 :
echo
echo "creating adjoint sources"
echo
~/SPECFEM3D_GLOBE/utils/adjoint_sources/traveltime/xcreate_adjsrc_traveltime $t1 $t2 3 ./$sta.MX*filt.ascii

# filter adjoint sources with the same bandpass used to filter seismograms
echo
echo "filtering adjoint sources"
echo
~/SPECFEM3D_GLOBE/utils/seis_process/process_syn.pl -t $tmin/$tmax -P 6/2 -x filt ./$sta.MX*.ascii.adj

rm -rf ./$sta.MX*adj.sac
rm -rf ./$sta.MX*filt.ascii
rm -rf ./$sta.MX*filt.ascii.adj

# convert sac files to ascii
sac2asc ./$sta.MXE*.adj.sac.filt > ./$sta.MXE.ascii.adj
sac2asc ./$sta.MXN*.adj.sac.filt > ./$sta.MXN.ascii.adj
sac2asc ./$sta.MXZ*.adj.sac.filt > ./$sta.MXZ.ascii.adj

rm -rf ./$sta.MX*sac.filt

# rename adjoint source files:
echo "renames file endings"
rename .ascii.adj .adj $sta.MX*adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT
cp -v ./STATIONS_ADJOINT ../DATA/

echo ""
echo "see adjoint sources: ./$sta.MX*.adj"
echo "done"

