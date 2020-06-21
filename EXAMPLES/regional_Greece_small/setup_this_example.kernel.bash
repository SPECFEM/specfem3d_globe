#!/bin/bash
#
# regional simulation example to compute an adjoint kernel
#
# script runs mesher and solver for an adjoint kernel using mpirun
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

if [ ! -f SEM/STATIONS_ADJOINT  ]; then
  echo "must have adjoint source station files in directory: SEM/"
  exit 1
fi
cp SEM/STATIONS_ADJOINT DATA/

mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES

rm -rf DATABASES_MPI/*
rm -rf OUTPUT_FILES/*

# compiles executables in root directory
# using default configuration
cd ../../

# compiles for an adjoint simulation
cp $currentdir/DATA/Par_file DATA/Par_file
sed -i "s:SAVE_FORWARD.*:SAVE_FORWARD                    = .true.:g"  DATA/Par_file
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
rm -rf bin/x*
cp -v ../../bin/x* ./bin/
cp -v ../../bin/xspecfem3D ./bin/xspecfem3D.kernel

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# links data directories needed to run example in this current directory with s362ani
cd DATA/
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/s362ani
ln -s ../../../DATA/QRFSI12
ln -s ../../../DATA/topo_bathy
cd ../

# copy useful script
if [ ! -f ./change_simulation_type.pl ]; then
ln -s ../../utils/change_simulation_type.pl
fi
cp DATA/Par_file DATA/Par_file.org

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo `date`

