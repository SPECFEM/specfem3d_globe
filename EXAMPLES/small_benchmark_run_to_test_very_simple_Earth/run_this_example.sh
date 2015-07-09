#!/bin/bash
#
# The script runs mesher and solver on 24 CPUs at long period, i.e. the benchmark is cheap.
# It does validate a lot of things though.
#
##################################################

echo "running example: `date`"
currentdir=`pwd`

echo "directory: $currentdir"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES

rm -rf DATABASES_MPI/*
rm -rf OUTPUT_FILES/*

# links data directories needed to run example in this current directory with s362ani
mkdir -p DATA
cd DATA/
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/s362ani
ln -s ../../../DATA/QRFSI12
ln -s ../../../DATA/topo_bathy
cd ../

# compiles executables in root directory
# using default configuration
cp DATA/Par_file ../../DATA
cd ../../
make clean
make all

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

# run mesher & solver
echo
echo "  running script..."
echo
./run_mesher_solver.bash

echo `date`

