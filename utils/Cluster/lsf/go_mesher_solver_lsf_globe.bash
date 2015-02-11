#!/bin/bash
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a mpich_gm
#BSUB -J go_mesher_solver_lsf

if [ -z $USER ]; then
  echo "could not run go_mesher_solver_...bash as no USER env is set"
  exit 2
fi

BASEMPIDIR=/scratch/$USER/DATABASES_MPI

# script to run the mesher and the solver

# read DATA/Par_file to get information about the run

# compute total number of nodes needed
NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

rm -r -f OUTPUT_FILES
mkdir OUTPUT_FILES

# obtain lsf job information
echo "$LSB_MCPU_HOSTS" > OUTPUT_FILES/lsf_machines
echo "$LSB_JOBID" > OUTPUT_FILES/jobid

remap_lsf_machines.pl OUTPUT_FILES/lsf_machines >OUTPUT_FILES/machines

# now cleanup and make the dir (seismos are now written by the master, no more need to collect them on the nodes), this for avoiding crashes
shmux -M 50 -S all -c "rm -r -f /scratch/$USER; mkdir -p /scratch/$USER; mkdir -p $BASEMPIDIR" - < OUTPUT_FILES/machines >/dev/null

echo starting MPI mesher on $numnodes processors
echo " "
echo starting run in current directory $PWD
echo " "

sleep 20
mpirun.lsf --gm-no-shmem --gm-copy-env $PWD/bin/xmeshfem3D

sleep 20
mpirun.lsf --gm-no-shmem --gm-copy-env $PWD/bin/xspecfem3D

# cleanup after the run
sleep 10
shmux -M 50 -S all -c "rm -r -f /scratch/$USER" - < OUTPUT_FILES/machines >/dev/null

