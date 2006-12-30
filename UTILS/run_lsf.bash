#!/bin/bash

# use the normal queue unless otherwise directed
queue="-q normal"
if [ $# -eq 1 ]; then
	echo "Setting the queue to $1"
	queue="-q $1"
fi

d=`date`
echo "Starting compilation $d"
make clean
make meshfem3D
make create_header_file
xcreate_header_file
make specfem3D
d=`date`
echo "Finished compilation $d"

# compute total number of nodes needed
NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep NCHUNKS DATA/Par_file | cut -d = -f 2`


# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

echo "Submitting job"
bsub $queue -n $numnodes -W 60 -K <go_mesher_solver_lsf_globe.bash

# for kernel simulation, use the following:
#bsub $queue -n $numnodes -W 200 -K <go_mesher_solver_lsf_globe.kernel

