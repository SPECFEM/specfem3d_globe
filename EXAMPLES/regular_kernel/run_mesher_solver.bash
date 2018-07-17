#!/bin/bash
#
# default, optional argument to be able to run only a kernel/backward simulation
do_mesher=1
if [ "$1" == "1" ]; then
# forward simulation
echo "forward simulation..."
do_mesher=1
elif [ "$1" == "3" ]; then
# kernel/backward simulation
echo "backward/kernel simulation..."
do_mesher=0
fi

# DATABASES directory
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

mkdir -p OUTPUT_FILES
mkdir -p $BASEMPIDIR

# backup files used for this simulation
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/

if [ "$do_mesher" == "1" ]; then
##
## mesh generation
##
sleep 2

echo
echo `date`
echo "starting MPI mesher on $numnodes processors"
echo

mpirun -np $numnodes $PWD/bin/xmeshfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "  mesher done: `date`"
echo

# backup important files addressing.txt and list*.txt
cp OUTPUT_FILES/*.txt $BASEMPIDIR/

fi


##
## forward simulation
##

# set up addressing
#cp $BASEMPIDIR/addr*.txt OUTPUT_FILES/
#cp $BASEMPIDIR/list*.txt OUTPUT_FILES/

sleep 2

echo
echo `date`
echo starting run in current directory $PWD
echo

mpirun -np $numnodes $PWD/bin/xspecfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "finished successfully"
echo `date`

