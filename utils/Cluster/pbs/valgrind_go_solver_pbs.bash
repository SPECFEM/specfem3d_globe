#!/bin/bash
#
# Valgrind, a suite of tools for debugging and profiling
# http://valgrind.org/
#

# bash script
#PBS -S /bin/bash

# job name
#PBS -N go_solver

# joins output and error information
#PBS -j oe

# job output file
#PBS -o OUTPUT_FILES/job.o

###########################################################
# USER PARAMETERS

# Queue
#PBS -q tromp

# 150 CPUs ( 18*8+6 ), walltime 15 hour
#PBS -l nodes=18:ppn=8+1:ppn=6,walltime=15:00:00

# valgrind mpi library
PRELOAD_LIB=/my_valgrind_path/valgrind/lib/valgrind/libmpiwrap-x86-linux.so

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

mkdir -p OUTPUT_FILES

# backup for files used for this simulation
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/

rm -rf OUTPUT_FILES/src
mkdir OUTPUT_FILES/src
cp -p *.f90 OUTPUT_FILES/src/
cp -p *.c OUTPUT_FILES/src/
cp -p *.in OUTPUT_FILES/src/
cp -p *.h OUTPUT_FILES/src/


# obtain job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

echo starting run in current directory $PWD
echo " "

cp $BASEMPIDIR/addr*.txt OUTPUT_FILES/
cp $BASEMPIDIR/list*.txt OUTPUT_FILES/

sleep 2

# memory leaks
LD_PRELOAD=$PRELOAD_LIB mpiexec -np $numnodes valgrind --leak-check=full $PWD/bin/xspecfem3D

cp OUTPUT_FILES/job.o OUTPUT_FILES/job.memory-leaks.o

sleep 2

# cache misses
LD_PRELOAD=$PRELOAD_LIB mpiexec -np $numnodes valgrind --tool=cachegrind $PWD/bin/xspecfem3D

cp OUTPUT_FILES/job.o OUTPUT_FILES/job.cache-misses.o
cp cachegrind.out.* OUTPUT_FILES/


echo "finished successfully"

