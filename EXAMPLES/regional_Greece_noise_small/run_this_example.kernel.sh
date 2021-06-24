#!/bin/bash
#
# regional simulation example to compute an adjoint kernel
#
# script runs mesher and solver for an adjoint kernel
# using batch scripts for a PBS queueing system
# on 4 CPUs
#
# modify accordingly for your own system specifics
##################################################

echo "running example: `date`"
currentdir=`pwd`

echo "directory: $currentdir"
echo "(will take about 40 minutes)"
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
# compiles for an adjoint simulation
cp $currentdir/DATA/Par_file DATA/Par_file
sed -i "s:SAVE_FORWARD.*:SAVE_FORWARD                    = .true.:"  DATA/Par_file
make clean
make -j4 all

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup of constants setup
cp setup/* $currentdir/OUTPUT_FILES/
cp OUTPUT_FILES/values_from_mesher.h $currentdir/OUTPUT_FILES/
cp DATA/Par_file $currentdir/OUTPUT_FILES/

cd $currentdir

# copy executables
mkdir -p bin
cp ../../bin/x* ./bin/

# links data directories needed to run example in this current directory with s362ani
cd DATA/
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/s362ani
ln -s ../../../DATA/QRFSI12
ln -s ../../../DATA/topo_bathy
cd ../

cp DATA/Par_file DATA/Par_file.org

## creates noise spectrum
./run_generate_S_squared.sh 2999 0.165

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# run mesher & solver
echo
echo "  running script..."
echo
./run_mesher_solver.kernel.bash

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "done `date`"
echo
