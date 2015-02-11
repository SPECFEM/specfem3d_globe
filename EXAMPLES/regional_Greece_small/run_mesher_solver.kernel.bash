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
cp DATA/STATIONS_ADJOINT OUTPUT_FILES/



##
## mesh generation
##
sleep 2
echo
echo `date`
echo "starting MPI mesher on $numnodes processors"
echo

mpirun -np $numnodes $PWD/bin/xmeshfem3D

echo "  mesher done: `date`"
echo

# backup important files addressing.txt and list*.txt
cp OUTPUT_FILES/*.txt $BASEMPIDIR/

##
## forward simulation
## (saving last wavefields)
##
cp DATA/Par_file DATA/Par_file.org
./change_simulation_type.pl -F
sleep 2
# set up addressing
#cp $BASEMPIDIR/addr*.txt OUTPUT_FILES/
#cp $BASEMPIDIR/list*.txt OUTPUT_FILES/

echo
echo `date`
echo "starting forward run in current directory $PWD"
echo

mpirun -np $numnodes $PWD/bin/xspecfem3D.kernel

echo "  forward run done: `date`"
echo

# renames output files of forward run
cd OUTPUT_FILES/
mv output_solver.txt output_solver.for.txt
rename .sem. .sem.for. *.sem.*
cd ../


##
## adjoint simulation
##
./change_simulation_type.pl -b
sleep 2

echo
echo `date`
echo "starting adjoint run in current directory $PWD"
echo

mpirun -np $numnodes $PWD/bin/xspecfem3D.kernel

echo "  adjoint run done: `date`"
echo

# restore original Par_file
mv DATA/Par_file.org DATA/Par_file

echo "finished successfully"
echo `date`

