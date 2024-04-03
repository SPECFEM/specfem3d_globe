#!/bin/bash

# this is the launch script to run a regular forward simulation
# on CITerra at Caltech with the LSF scheduler
# Qinya Liu, Caltech, May 2007

# use the normal queue unless otherwise directed
queue="-q normal"
if [ $# -eq 1 ]; then
  echo "Setting the queue to $1"
  queue="-q $1"
fi

d=`date`
echo "Starting compilation $d"
make clean
make xmeshfem3D
make xcreate_header_file
./bin/xcreate_header_file
make xspecfem3D
d=`date`
echo "Finished compilation $d"

# compute total number of nodes needed
NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2`
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep NCHUNKS DATA/Par_file | cut -d = -f 2`

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

#rm -r -f OUTPUT_FILES/*

echo "Submitting job"

# time below is given in hh:mm
bsub $queue -n $numnodes -W 48:00 -C 0 < go_mesher_solver_lsf_globe.bash

