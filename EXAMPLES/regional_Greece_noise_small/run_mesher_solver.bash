#!/bin/bash

###########################################################
# USER PARAMETERS

## 4 CPUs
CPUs=4

###########################################################


BASEMPIDIR=`grep LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

if [ ! "$numnodes" == "$CPUs" ]; then
  echo "error: Par_file for $numnodes CPUs"
  exit 1
fi

mkdir -p OUTPUT_FILES

# backup files used for this simulation
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/

cp DATA/Par_file DATA/Par_file.org


##
## mesh generation
##
sleep 2

echo
echo `date`
echo "starting MPI mesher on $numnodes processors"
echo

mpirun -np $numnodes $PWD/bin/xmeshfem3D
output=$?
if [ ! "$output" == "0" ]; then
  exit 1
fi

echo "  mesher done: $output `date`"
echo

# backup important files addressing.txt and list*.txt
cp OUTPUT_FILES/*.txt $BASEMPIDIR/


##
## forward simulation - ensemble forward source
##
# set up addressing
#cp $BASEMPIDIR/addr*.txt OUTPUT_FILES/
#cp $BASEMPIDIR/list*.txt OUTPUT_FILES/
sleep 1
echo
echo `date`
echo starting 0.1 run in current directory $PWD
echo

sed -i "s:SIMULATION_TYPE .*:SIMULATION_TYPE                 = 1:g" DATA/Par_file
sed -i "s:NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY                = 1:g" DATA/Par_file
sed -i "s:SAVE_FORWARD .*:SAVE_FORWARD                    = .true.:g" DATA/Par_file


mpirun -np $numnodes $PWD/bin/xspecfem3D
output=$?
if [ ! "$output" == "0" ]; then
  exit 1
fi

echo "solver done: $output `date`"
echo "moving OUTPUT_FILES..."

rm -rf OUTPUT_FILES_0.1
mv OUTPUT_FILES OUTPUT_FILES_0.1
mkdir -p OUTPUT_FILES
cp $BASEMPIDIR/addr*.txt OUTPUT_FILES/
cp $BASEMPIDIR/list*.txt OUTPUT_FILES/

##
## forward simulation - ensemble forward wavefield
##
sleep 1
echo
echo `date`
echo starting 0.2 run in current directory $PWD
echo

sed -i "s:NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY                = 2:g" DATA/Par_file

mpirun -np $numnodes $PWD/bin/xspecfem3D
output=$?
if [ ! "$output" == "0" ]; then
  exit 1
fi

echo "solver done: $output `date`"
echo "moving OUTPUT_FILES..."

rm -rf OUTPUT_FILES_0.2
mv OUTPUT_FILES OUTPUT_FILES_0.2
mkdir -p OUTPUT_FILES
cp $BASEMPIDIR/addr*.txt OUTPUT_FILES/
cp $BASEMPIDIR/list*.txt OUTPUT_FILES/


##
## prepare adjoint source
##
echo `date`
echo "prepare adjoint source..."

./SEM/create_adjoint_source.sh

echo "done"
echo `date`


