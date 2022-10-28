#!/bin/bash


# CONSTANTS
# ==================
# example directory
currentdir=`pwd`

# gets settings from Par_file
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))


# SETUP ENVIROMENT
# ==================
echo "running example: `date`"

echo
echo "setting up example"
echo

# remove trash
rm -rf DATA/topo_bathy
rm -rf bin DATABASES_MPI OUTPUT_FILES OUTPUT_FILES_1 OUTPUT_FILES_2
rm -f NOISE_TOMOGRAPHY.py PetersonNoiseModel.m

# make directories
mkdir bin
mkdir DATABASES_MPI
mkdir OUTPUT_FILES

# link executables
cd bin
ln -s ../../../bin/* .
cd $currentdir

# link earth model and topography
cd DATA/
ln -s ../../../DATA/topo_bathy .
cd $currentdir

# link utilities
ln -s ../noise_examples/NOISE_TOMOGRAPHY.py
ln -s ../noise_examples/PetersonNoiseModel.m

# generate noise source time function
echo
echo " generating noise source time function..."
echo
./run_generate_S_squared.sh
if [[ $? -ne 0 ]]; then exit 1; fi

# for this example, the noise distribution was generated
# using "generate_noise_distribution_direction.py" and
# stored in NOISE_TOMOGRAPHY/noise_distribution, therefore,
# copy the files to DATABASES_MPI
cp NOISE_TOMOGRAPHY/noise_distribution/* DATABASES_MPI/

# RUN MESHER
# ==================
# run mesher
echo
echo `date`
echo "starting MPI mesher on $numnodes processors"
echo

mpirun -np $numnodes $PWD/bin/xmeshfem3D
if [[ $? -ne 0 ]]; then exit 1; fi

# backup simulation files
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp OUTPUT_FILES/*.txt $BASEMPIDIR/

echo "  mesher done: `date`"
echo


# RUN NOISE STEP 1
# ==================
# run solver
cp DATA/Par_file_step1 DATA/Par_file
echo
echo " running noise simulation step 1 on $numnodes processors..."
echo
mpirun -np $numnodes ./bin/xspecfem3D
if [[ $? -ne 0 ]]; then exit 1; fi

# save generating field movie
echo
echo " saving generating field movie..."
echo
./bin/xcreate_movie_AVS_DX < generate_movie_input.txt
if [[ $? -ne 0 ]]; then exit 1; fi
mkdir OUTPUT_FILES_1
cp OUTPUT_FILES/* OUTPUT_FILES_1


# RUN NOISE STEP 2
# ==================
# run solver
cp DATA/Par_file_step2 DATA/Par_file
echo
echo " running noise simulation step 2 on $numnodes processors..."
echo
mpirun -np $numnodes ./bin/xspecfem3D
if [[ $? -ne 0 ]]; then exit 1; fi

# save generating field movie
echo
echo " saving correlation field movie..."
echo
./bin/xcreate_movie_AVS_DX < generate_movie_input.txt
if [[ $? -ne 0 ]]; then exit 1; fi
mkdir OUTPUT_FILES_2
cp OUTPUT_FILES/* OUTPUT_FILES_2
rm -rf OUTPUT_FILES

echo
date
echo "done"
