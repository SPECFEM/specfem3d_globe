#!/bin/bash

# gets settings from Par_file

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
BASEMPIDIR_NAME=`basename $BASEMPIDIR`

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 `

NSIM=`grep ^NUMBER_OF_SIMULTANEOUS_RUNS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))
numnodes_total=$(( $NSIM * $numnodes ))

echo
echo "setting up run directories..."
echo "  number of simultaneous runs = $NSIM"
echo "  number of processes for single simulation run  = $numnodes"
echo "  total number of processes for simultaneous run = $numnodes_total"
echo

currentdir=`pwd`

mkdir -p OUTPUT_FILES

# backup files used for this simulation
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
#cp DATA/CMTSOLUTION OUTPUT_FILES/

##
## mesh generation
##
sleep 2

echo
echo `date`
echo "starting MPI mesher on $numnodes processors"
echo

echo "running mesher in main directory: $currentdir"

# setup for mesher
# cleans run directory run0001
if [ -d "run0001" ]; then
  rm -rf run0001/*
fi
cp DATA/Par_file DATA/Par_file.org

# sets to single run for meshing only
sed -i "s:^NUMBER_OF_SIMULTANEOUS_RUNS .*:NUMBER_OF_SIMULTANEOUS_RUNS = 1:" DATA/Par_file
# will need a CMTSOLUTION file in DATA/
cp DATA/CMTSOLUTION.1 DATA/CMTSOLUTION

# mesher (on 1 MPI process)
mpirun -np $numnodes $PWD/bin/xmeshfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "  mesher done: `date`"
echo

# backup important files addressing.txt and list*.txt
cp OUTPUT_FILES/*.txt $BASEMPIDIR/

# re-sets to original Par_file
mv -v DATA/Par_file.org DATA/Par_file
rm -f DATA/CMTSOLUTION


## creates simultaneous run directories
for ((i=1; i<=$NSIM; i++)); do
  it=$(printf "%04d" $i)
  dir=run${it}

  echo
  echo "creating run directory: $dir"
  mkdir -p $dir
  rm -rf $dir/*
  cd $dir/

  # creates DATA/
  mkdir -pv DATA/
  rm -f DATA/*
  # only needs STATIONS and CMTSOLUTION for this event
  # (the Par_file for all runs will be in root directory ../DATA/)
  cp -vp ../DATA/STATIONS DATA/
  cp -vp ../DATA/CMTSOLUTION.$i DATA/CMTSOLUTION
  if [[ $? -ne 0 ]]; then exit 1; fi

  # cleans output files
  mkdir -pv OUTPUT_FILES
  rm -rf OUTPUT_FILES/*

  # sets up DATABASES_MPI/
  mkdir -pv $BASEMPIDIR
  rm -rf $BASEMPIDIR/*

  # solver needs mesh files in run0001/
  if [ "$i" == "1" ]; then
    cd $BASEMPIDIR/
    # creates symbolic links to mesh files
    # (still useful to have a local DATABASES_MPI/ folder, for example to store kernel files for this event)
    # mesh parameter file
    if [ -e ../../${BASEMPIDIR_NAME}/mesh_parameters.bin ]; then
      ln -s ../../${BASEMPIDIR_NAME}/mesh_parameters.bin
    else
      echo "error didn't find file : ../../${BASEMPIDIR_NAME}/mesh_parameters.bin"
      echo "please check if mesher ran in main directory..."
      exit 1
    fi
    # topo
    if [ -e ../../${BASEMPIDIR_NAME}/mesh_topo_bathy.bin ]; then
      ln -s ../../${BASEMPIDIR_NAME}/mesh_topo_bathy.bin
    fi
    # proc** files
    model_files=`ls ../../${BASEMPIDIR_NAME}/proc*.bin`
    for f in ${model_files}; do
      ln -s $f
    done
    # adios files
    if [ -e ../../${BASEMPIDIR_NAME}/attenuation.bp ]; then
      ln -s ../../${BASEMPIDIR_NAME}/attenuation.bp
    fi
    if [ -e ../../${BASEMPIDIR_NAME}/boundary.bp ]; then
      ln -s ../../${BASEMPIDIR_NAME}/boundary.bp
    fi
    if [ -e ../../${BASEMPIDIR_NAME}/solver_data.bp ]; then
      ln -s ../../${BASEMPIDIR_NAME}/solver_data.bp
    fi
    if [ -e ../../${BASEMPIDIR_NAME}/solver_data_mpi.bp ]; then
      ln -s ../../${BASEMPIDIR_NAME}/solver_data_mpi.bp
    fi
    if [ -e ../../${BASEMPIDIR_NAME}/stacey.bp ]; then
      ln -s ../../${BASEMPIDIR_NAME}/stacey.bp
    fi

    cd ../
    # checks exit code
    if [[ $? -ne 0 ]]; then exit 1; fi
  fi

  # stores setup
  cp DATA/CMTSOLUTION OUTPUT_FILES/
  cp DATA/STATIONS OUTPUT_FILES/
  cp ../DATA/Par_file OUTPUT_FILES/

  # back to root directory
  cd ../
done


##
## forward simulation
##
sleep 2

echo
echo `date`
echo "starting run in current directory $PWD"
echo

# must remove either DATA/Par_file or run0001/DATA/Par_file for simultaneous runs
# (assumes to have DATA/Par_file, moves other one)
if [ -e run0001/DATA/Par_file ]; then mv -v run0001/DATA/Par_file run0001/DATA/Par_file.org; fi

# solver (on all 4 MPI processes)
mpirun -np $numnodes_total $PWD/bin/xspecfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "finished successfully"
echo `date`

