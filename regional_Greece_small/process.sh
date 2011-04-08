#!/bin/bash
#
# regional simulation example
#
# script runs mesher and solver using mpirun
# on 4 CPUs
#
# synthetics have an approximate shortest period ~ 15 s
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
# configures package with ifort compiler
./configure F90=ifort MPIF90=/usr/local/openmpi-ifort/bin/mpif90 FLAGS_CHECK="-O3 -assume byterecl" FLAGS_NO_CHECK="-O3 -assume byterecl" > tmp.log
# configures package with gfortran compiler
#./configure F90=gfortran MPIF90=mpif90 FLAGS_CHECK="-O3" FLAGS_NO_CHECK="-O3"

# compiles for a forward simulation
cp $currentdir/DATA/Par_file DATA/Par_file
make clean
make >& $currentdir/tmp_make_output.log
make xcombine_vol_data >> $currentdir/tmp_make_output.log

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

# links data directories needed to run example in this current directory with s362ani
cd DATA/
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/s362ani
ln -s ../../../DATA/QRFSI12
ln -s ../../../DATA/topo_bathy
cd ../

# run mesher & solver
echo
echo "  running script..."
echo
./run_mesher_solver.bash


echo `date`

