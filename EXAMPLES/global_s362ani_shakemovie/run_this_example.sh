#!/bin/bash
#
# global simulation example
#
# script runs mesher and solver
# using batch scripts for a PBS queueing system
# on 24 CPUs
#
#
# modify accordingly for your own system specifics
##################################################
# USER PARAMETER
# source directory
rootdir=~/SPECFEM3D_GLOBE_GPU

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

# compiles executables in root directory
# using default configuration

cd $rootdir
# compiles for a forward simulation
cp $currentdir/DATA/Par_file DATA/Par_file
make clean
make all

# backup of constants setup
cp setup/* $currentdir/OUTPUT_FILES/
cp DATA/Par_file $currentdir/OUTPUT_FILES/

cd $currentdir

# copy executables
mkdir -p bin
cp $rootdir/bin/xmeshfem3D ./bin/
cp $rootdir/bin/xspecfem3D ./bin/

# links data directories needed to run example in this current directory with s362ani
cd DATA/
ln -s $rootdir/DATA/crust2.0
ln -s $rootdir/DATA/s362ani
ln -s $rootdir/DATA/QRFSI12
ln -s $rootdir/DATA/topo_bathy
cd ../

# submits job to run mesher & solver
echo
echo "  submitting script..."
echo

echo "please submit job now manually: "
echo "  meshing            : qsub go_mesher_pbs.bash"
echo "  forward simulation : qsub go_solver_pbs.bash"
echo
echo "after job completion, see results in directory: OUTPUT_FILES/"
echo "done: `date`"

