#!/bin/bash

# gets settings from Par_file

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

# backup files used for this simulation
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/

cp DATA/Par_file DATA/Par_file.org
# clean
rm -rf OUTPUT_FILES_1 OUTPUT_FILES_2 OUTPUT_FILES_3

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

# backup important files addressing.txt
cp OUTPUT_FILES/*.txt $BASEMPIDIR/


##
## forward simulation - ensemble forward source
##
sleep 1

echo
echo `date`
echo starting 1. run in current directory $PWD
echo

sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE                 = 1:g" DATA/Par_file
sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY                = 1:g" DATA/Par_file
sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD                    = .false.:g" DATA/Par_file

# main noise source: receiver id for station S001
echo "2" > NOISE_TOMOGRAPHY/irec_main_noise

mpirun -np $numnodes $PWD/bin/xspecfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "solver done: `date`"

mv -v OUTPUT_FILES OUTPUT_FILES_1
mkdir -p OUTPUT_FILES
echo


##
## forward simulation - ensemble forward wavefield
##
sleep 1

echo
echo `date`
echo starting 2. run in current directory $PWD
echo

sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE                 = 1:g" DATA/Par_file
sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY                = 2:g" DATA/Par_file
sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD                    = .true.:g" DATA/Par_file

mpirun -np $numnodes $PWD/bin/xspecfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "solver done: `date`"

mv -v OUTPUT_FILES OUTPUT_FILES_2
mkdir -p OUTPUT_FILES
echo

##
## prepare adjoint source
##
echo `date`
echo "prepare adjoint source..."
cd SEM/
./create_adjoint_source.sh ../OUTPUT_FILES_2
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
cd ../
echo

##
## adjoint simulation - ensemble kernel
##
sleep 1

echo
echo `date`
echo starting 3. run in current directory $PWD
echo

sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE                 = 3:g" DATA/Par_file
sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY                = 3:g" DATA/Par_file
sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD                    = .false.:g" DATA/Par_file

# single adjoint source
echo "1" > NOISE_TOMOGRAPHY/irec_main_noise

mpirun -np $numnodes $PWD/bin/xspecfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "solver done: `date`"

mv -v OUTPUT_FILES OUTPUT_FILES_3
#mkdir -p OUTPUT_FILES
echo

echo "finished successfully"
echo `date`

