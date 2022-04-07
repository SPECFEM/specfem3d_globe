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

# checks if executables were compiled and available
if [ ! -e ../../bin/xspecfem3D ]; then
  echo "Compiling first all binaries in the root directory..."
  echo

  # compiles executables in root directory
  # using default configuration
  cd ../../

  # only in case static compilation would have been set to yes in Makefile:
  cp $currentdir/DATA/Par_file DATA/Par_file

  # compiles code
  make clean
  make -j4 all

  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  # backup of constants setup
  cp setup/* $currentdir/OUTPUT_FILES/
  if [ -e OUTPUT_FILES/values_from_mesher ]; then
    cp OUTPUT_FILES/values_from_mesher.h $currentdir/OUTPUT_FILES/values_from_mesher.h.compilation
  fi
  cp DATA/Par_file $currentdir/OUTPUT_FILES/

  cd $currentdir
fi

# copy executables
mkdir -p bin
rm -rf bin/*
cp ../../bin/xmeshfem3D ./bin/
cp ../../bin/xspecfem3D ./bin/
cp ../../bin/xcombine_vol_data* ./bin/

# links data directories needed to run example in this current directory with s362ani
cd DATA/
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/crustmap
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

