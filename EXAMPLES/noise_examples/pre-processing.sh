#!/bin/bash -eu
################################################

# example directory
example=$1

# compiler
f90=ifort

# submit
do_submit=$2

################################################

if [ "$1" == "" ]; then echo "usage: ./pre-processing.sh directory[e.g. global_long] [1]"; exit 1; fi

dir=$PWD/$example
if [ ! -d $dir ]; then echo "Error directory $dir doesn't exist, please check..."; exit 1; fi

echo ""
echo ""
echo "you are going to run the example:"
echo "*****$example******"
echo ""
echo ""

case $example in
  global_long)
    name=zzz_global_long
    NSTEP=31599
    nodes=75
    ;;
  global_short)
    name=zzz_global_short
    NSTEP=31599
    nodes=75
    ;;
  regional)
    name=zzz_regional
    NSTEP=31599
    nodes=50
    ;;
  test_global)
    name=zzz_test_global
    NSTEP=799
    nodes=75
    ;;
  test_regional)
    name=zzz_test_regional
    NSTEP=799
    nodes=50
    ;;
  *)
    echo "error unknown directory: $example"
    echo "usage: ./pre-processing.sh directory[e.g. global_long] [1]";
    exit 1
esac


# change to root directory
cd ../../

rm -rf   SEM NOISE_TOMOGRAPHY $name
mkdir -p SEM NOISE_TOMOGRAPHY $name/job_info

cp $dir/adj_traveltime_filter.f90        ./
$f90 adj_traveltime_filter.f90 -o NOISE_adj

script=job_NOISE123_pbs.submit_atten.sh

cp EXAMPLES/noise_examples/scripts/$script ./
cp EXAMPLES/noise_examples/scripts/run_NOISE_clean_create.sh ./
cp EXAMPLES/noise_examples/scripts/run_NOISE_collect.sh      ./

sed -i "s/zzz_global_long/$name/g" ./$script
sed -i "s/^NSTEP=.*/NSTEP=$NSTEP/g" ./$script
sed -i "s/^#PBS -l nodes=.*/#PBS -l nodes=$nodes:ppn=8,walltime=100:00:00/g" ./$script

cp $dir/S_squared                        ./NOISE_TOMOGRAPHY/
cp $dir/irec_master_noise                ./NOISE_TOMOGRAPHY/
cp $dir/nu_master                        ./NOISE_TOMOGRAPHY/

cp $dir/Par_file_NOISE_*                 ./DATA/
cp $dir/CMTSOLUTION_NOISE                ./DATA/CMTSOLUTION
cp $dir/STATIONS_NOISE                   ./DATA/STATIONS
cp $dir/STATIONS_NOISE_ADJOINT           ./DATA/STATIONS_ADJOINT

cp DATA/Par_file_NOISE_1_attenuation DATA/Par_file
rm -rf OUTPUT_FILES/*
make xmeshfem3D
make xspecfem3D
make xcreate_movie_AVS_DX
mv bin/xmeshfem3D xmeshfem3D_attenuation
mv bin/xspecfem3D xspecfem3D_attenuation
mv bin/xcreate_movie_AVS_DX xcreate_movie_AVS_DX_attenuation

cp DATA/Par_file_NOISE_1_noattenuation DATA/Par_file
rm -rf OUTPUT_FILES/*
make xmeshfem3D
make xspecfem3D
make xcreate_movie_AVS_DX
make xcombine_vol_data
mv bin/xmeshfem3D xmeshfem3D_noattenuation
mv bin/xspecfem3D xspecfem3D_noattenuation
mv bin/xcreate_movie_AVS_DX xcreate_movie_AVS_DX_noattenuation
mv bin/xcombine_vol_data    xcombine_vol_data

# submit job-script
if [ "$do_submit" == "1" ]; then
  echo ""
  echo "submitting pbs-job:"
  qsub $script
fi

