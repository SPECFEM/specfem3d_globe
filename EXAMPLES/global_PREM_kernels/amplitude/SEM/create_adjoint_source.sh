#!/bin/bash
#
# creates amplitude adjoint source files
#
# notes: 1- run this script in directory SEM/
#        2- make sure you have SAC to use proces_syn.pl script
#        3- make sure you have compiled create_adjsrc_amplitude.f90
#          under ~/SPECFEM3D_GLOBE/UTILS/adjoint_sources/amplitude
#############################################################

# station
station="ANMO"
network="IU"

# time window
t1="500.0"
t2="600.0"

# corner periods
tmin="17"
tmax="400"

# filter synthetics within a bandpass filter with corner periods of f1 and f2
~/SPECFEM3D_GLOBE/UTILS/seis_process/process_syn.pl -t $tmin/$tmax -P 6/2 -x filt ../REF_SEIS/$station.$network.MX*

rm -rf ../REF_SEIS/$station.$network.MX*sac

# convert sac files to ascii
sac2asc ../REF_SEIS/$station.$network.MXE*filt > ../REF_SEIS/$station.$network.MXE.filt.ascii  
sac2asc ../REF_SEIS/$station.$network.MXN*filt > ../REF_SEIS/$station.$network.MXN.filt.ascii  
sac2asc ../REF_SEIS/$station.$network.MXZ*filt > ../REF_SEIS/$station.$network.MXZ.filt.ascii  

rm -rf ../REF_SEIS/$station.$network.MX*sac.filt

# window out single phase arrival on vertical component between t1 to t2 :
~/SPECFEM3D_GLOBE/UTILS/adjoint_sources/amplitude/xcreate_adjsrc_amplitude $t1 $t2 3 ../REF_SEIS/$station.$network.MX*filt.ascii

# filter adjoint sources with the same bandpass used to filter seismograms
~/SPECFEM3D_GLOBE/UTILS/seis_process/process_syn.pl -t $tmin/$tmax -P 6/2 -x filt ../REF_SEIS/$station.$network.MX*.ascii.adj

rm -rf ../REF_SEIS/$station.$network.MX*adj.sac
rm -rf ../REF_SEIS/$station.$network.MX*filt.ascii
rm -rf ../REF_SEIS/$station.$network.MX*filt.ascii.adj

# convert sac files to ascii
sac2asc ../REF_SEIS/$station.$network.MXE*.adj.sac.filt > ../REF_SEIS/$station.$network.MXE.ascii.adj
sac2asc ../REF_SEIS/$station.$network.MXN*.adj.sac.filt > ../REF_SEIS/$station.$network.MXN.ascii.adj
sac2asc ../REF_SEIS/$station.$network.MXZ*.adj.sac.filt > ../REF_SEIS/$station.$network.MXZ.ascii.adj

rm -rf ../REF_SEIS/$station.$network.MX*sac.filt

mv ../REF_SEIS/$station.$network.MX*adj ./

# rename adjoint source files:
rename .ascii.adj .adj $station.$network.MX*adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT

