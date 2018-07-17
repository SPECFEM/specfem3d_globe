#!/bin/bash
#
# global simulation example with a point force source
#
# script runs mesher and solver using mpirun
#
# modify accordingly for your own system specifics
##################################################

echo "running example: `date`"
currentdir=`pwd`

echo "directory: $currentdir"
echo "(will take about 7 minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES

rm -rf DATABASES_MPI/*
rm -rf OUTPUT_FILES/*

# compiles executables in root directory
# using default configuration
cd ../../

# compiles for a forward simulation
cp $currentdir/DATA/Par_file DATA/Par_file
cp $currentdir/DATA/CMTSOLUTION DATA/
make clean
make -j4 all

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup of constants setup
cp setup/* $currentdir/OUTPUT_FILES/
cp OUTPUT_FILES/values_from_mesher.h $currentdir/OUTPUT_FILES/values_from_mesher.h.compilation
cp DATA/Par_file $currentdir/OUTPUT_FILES/

cd $currentdir

# copy executables
mkdir -p bin
rm -rf bin/*
cp ../../bin/xmeshfem3D ./bin/
cp ../../bin/xspecfem3D ./bin/
cp ../../bin/xcombine_vol_data ./bin/
cp ../../bin/xcombine_vol_data_vtk ./bin/

# links data directories needed to run example in this current directory with s362ani
cd DATA/
ln -s ../../../DATA/s40rts
ln -s ../../../DATA/s20rts
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/topo_bathy
cd ../

# run mesher & solver
echo
echo "  running script..."
echo
./run_mesher_solver.bash

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo `date`

