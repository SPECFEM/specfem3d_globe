#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_mesher
#PBS -j oe
#PBS -o OUTPUT_FILES/job.o

###########################################################
# USER PARAMETERS

## 150 CPUs ( 18*8 + 6 ), walltime 1 hour
#PBS -l nodes=18:ppn=8+1:ppn=6,walltime=1:00:00

###########################################################

cd $PBS_O_WORKDIR

BASEMPIDIR=`grep LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

# cleans OUTPUT_FILES
rm -r -f OUTPUT_FILES
mkdir OUTPUT_FILES

# backup for files used for this simulation
cp DATA/Par_file OUTPUT_FILES/

rm -rf OUTPUT_FILES/src
mkdir OUTPUT_FILES/src
cp -p *.f90 OUTPUT_FILES/src/
cp -p *.c OUTPUT_FILES/src/
cp -p *.in OUTPUT_FILES/src/
cp -p *.h OUTPUT_FILES/src/

# obtain job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid


echo starting MPI mesher on $numnodes processors
echo " "

sleep 2
mpiexec -np $numnodes $PWD/bin/xmeshfem3D

# backup important files addressing.txt and list*.txt
cp OUTPUT_FILES/*.txt $BASEMPIDIR/


echo "done"
