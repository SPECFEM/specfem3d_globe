#!/bin/bash
#
# regional simulation example
#
# script runs mesher and solver
# using batch scripts for a PBS queueing system
# on 64 CPUs
#
# synthetics have an approximate shortest period ~ 17 s
#
# modify accordingly for your own system specifics
##################################################

echo "running example: `date`"
currentdir=`pwd`

echo "directory: $currentdir"
echo "(will take about 15 minutes)"
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
make clean
make all

# backup of constants setup
cp setup/* $currentdir/OUTPUT_FILES/
cp OUTPUT_FILES/values_from_mesher.h $currentdir/OUTPUT_FILES/
cp DATA/Par_file $currentdir/OUTPUT_FILES/

cd $currentdir

# copy executables
mkdir -p bin
cp ../../bin/xmeshfem3D ./bin/
cp ../../bin/xspecfem3D ./bin/

# links data directories needed to run example in this current directory with s362ani
cd DATA/
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/s362ani
ln -s ../../../DATA/QRFSI12
ln -s ../../../DATA/topo_bathy
cd ../

# submits job to run mesher & solver
echo
echo "  submitting script..."
echo
first=`qsub go_mesher_solver_pbs.bash`
echo "  submitted job: $first"

echo
echo "after job completion, see results in directory: OUTPUT_FILES/"
echo "done submission setup"
echo `date`

