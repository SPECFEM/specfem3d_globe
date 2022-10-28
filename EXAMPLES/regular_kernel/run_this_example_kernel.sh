#!/bin/bash
#
# global simulation example
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

ln -s ../../utils/change_simulation_type.pl
./change_simulation_type.pl -F

# DATABASES directory
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# checks if executables were compiled and available
if [ ! -e ../../bin/xspecfem3D ]; then
  echo "Compiling first all binaries in the root directory..."
  echo

  # compiles executables in root directory
  # using default configuration
  cd ../../

  # only in case static compilation would have been set to yes in Makefile:
  cp $currentdir/DATA/Par_file DATA/Par_file
  sed -i "s:SAVE_FORWARD.*:SAVE_FORWARD                    = .true.:"  DATA/Par_file

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
cp ../../bin/xcombine_vol_data ./bin/
cp ../../bin/xcombine_vol_data_vtk ./bin/

# links data directories needed to run example in this current directory with s362ani
cd DATA/
ln -s ../../../DATA/s40rts
ln -s ../../../DATA/s20rts
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/topo_bathy
cd ../


echo
echo "#########################################################"
echo "forward simulation"
echo "#########################################################"
echo "(running forward simulation with saving forward wavefield)"
echo
./change_simulation_type.pl -F
# backup
cp DATA/Par_file OUTPUT_FILES/Par_file.for

# run mesher & solver
echo
echo "  running script..."
echo
./run_mesher_solver.bash 1

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "#########################################################"
echo "adjoint sources"
echo "#########################################################"
echo "setting up adjoint sources"
echo
# setup adjoint sources directory
mkdir -p SEM
rm -rf SEM/*

./create_adjoint_sources.sh

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup forward outputs
mkdir -p OUTPUT_FILES/forward_run
rm -rf OUTPUT_FILES/forward_run/*

mv OUTPUT_FILES/timestamp* OUTPUT_FILES/forward_run/
mv OUTPUT_FILES/output_* OUTPUT_FILES/forward_run/
mv OUTPUT_FILES/*.sem.* OUTPUT_FILES/forward_run/
mv OUTPUT_FILES/plot_* OUTPUT_FILES/sr.vtk OUTPUT_FILES/forward_run/

echo "#########################################################"
echo "kernel simulation"
echo "#########################################################"
echo "(running kernel simulation: SIMULATION_TYPE == 3)"
echo
./change_simulation_type.pl -b
# stores output
cp DATA/Par_file DATA/Par_file.kernel
cp DATA/CMTSOLUTION DATA/STATIONS* OUTPUT_FILES/

# only run solver
echo
echo "  running script..."
echo
./run_mesher_solver.bash 3

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "see results in directory       : OUTPUT_FILES/"
echo "    kernel outputs in directory: $BASEMPIDIR"
echo
echo "done"
echo `date`

