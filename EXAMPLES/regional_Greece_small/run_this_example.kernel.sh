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
  exit
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
cp ../../bin/xspecfem3D ./bin/xspecfem3D.kernel
cp ../../bin/xcombine_vol_data ./bin/
cp ../../bin/xcombine_vol_data_vtk ./bin/

# links data directories needed to run example in this current directory with s362ani
cd DATA/
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/s362ani
ln -s ../../../DATA/QRFSI12
ln -s ../../../DATA/topo_bathy
cd ../

# copy useful script
cp ../../utils/change_simulation_type.pl ./
cp DATA/Par_file DATA/Par_file.org

# run mesher & solver
echo
echo "  running script..."
echo
./run_mesher_solver.kernel.bash

echo `date`

