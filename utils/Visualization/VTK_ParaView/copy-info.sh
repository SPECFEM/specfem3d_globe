#!/bin/sh

# copy-info.sh
# SPECFEM3D_GLOBE
#

PWD=`pwd`
echo "$PWD"
echo "copy cmt solution etc..."

jobid=$1

#echo "calculating slice files..."
#./global_slice_number.pl ../../DATA/CMTSOLUTION ../../DATA/STATIONS_ADJOINT ../../DATA/Par_file
#echo

echo "collecting database files..."
../collect_database/copy_m_globe_database.pl slices_all ../../OUTPUT_FILES/lsf_machines beta_kernel $jobid

echo
cd ../../
make xcombine_vol_data
cd $PWD

echo "combining mesh..."
../../xcombine_vol_data UTILS/Paraview/slices_all beta_kernel $PWD $PWD ../../OUTPUT_FILES 1 1

## echo "converting mesh to vtu..."
#./mesh2vtu.pl -i ../../OUTPUT_FILES/reg_1_beta_kernel.mesh ../../OUTPUT_FILES/reg_1_beta_kernel.vtu


echo "done"
echo
