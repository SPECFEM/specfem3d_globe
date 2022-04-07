#!/bin/bash
################################################

# example directory
example=$1

# compiler
f90=gfortran  #ifort

# submit
do_submit=$2

################################################

if [ "$1" == "" -o "$2" == "" ]; then
  echo "usage: ./run_this_example.sh directory[e.g. global_long] submit-job[on=1/off=0]"
  exit 1
fi

# example directory
example=`basename $example`
dir=$PWD/$example
if [ ! -d $dir ]; then echo "Error: directory $dir doesn't exist, please check..."; exit 1; fi

echo ""
echo ""
echo "you are going to run the example:"
echo "***** $example ******"
echo ""
echo ""
sleep 2

# example setup settings
case $example in
  global_long)
    name=zzz_global_long
    NSTEP=31599
    nodes=75
    #adjoint measurement
    filter=1
    windowlength=2000.d0
    ;;
  global_short)
    name=zzz_global_short
    NSTEP=31599
    nodes=75
    #adjoint measurement
    filter=1
    windowlength=2000.d0
    ;;
  regional)
    name=zzz_regional
    NSTEP=31599
    nodes=50
    #adjoint measurement
    filter=1
    windowlength=2000.d0
    ;;
  test_global)
    name=zzz_test_global
    NSTEP=799
    nodes=75
    #adjoint measurement
    filter=0
    windowlength=60.d0
    ;;
  test_regional)
    name=zzz_test_regional
    NSTEP=799
    nodes=50
    #adjoint measurement
    filter=0
    windowlength=60.d0
    ;;
  *)
    echo "error unknown directory: $example"
    echo "usage: ./run_this_example.sh directory[e.g. global_long] job-submit[1]";
    exit 1
esac


# change to example directory
cd $dir
echo
echo "directory: `pwd`"
echo
currentdir=`pwd`

rm -rf   SEM $name bin
mkdir -p SEM $name/job_info bin OUTPUT_FILES

# compiles adjoint source tool
cp ../adj_traveltime_filter.f90 ./
sed -i "s/length_time_window = .*/length_time_window = $windowlength/" ./adj_traveltime_filter.f90
sed -i "s/filter_flag = .*/filter_flag = $filter/" ./adj_traveltime_filter.f90

$f90 adj_traveltime_filter.f90 -o NOISE_adj

# copies PBS job-scripts
script=job_NOISE123_pbs.submit_atten.sh

cp -v ../scripts/$script ./
cp -v ../scripts/run_NOISE_clean_create.sh ./
cp -v ../scripts/run_NOISE_collect.sh      ./

sed -i "s/zzz_global_long/$name/g" ./$script
sed -i "s/^NSTEP=.*/NSTEP=$NSTEP/g" ./$script
sed -i "s/^#PBS -l nodes=.*/#PBS -l nodes=$nodes:ppn=8,walltime=100:00:00/g" ./$script


# attenuation simulations
echo; echo
cp -v DATA/Par_file_NOISE_1_attenuation DATA/Par_file

# checks if executables were compiled and available
if [ ! -e ../../../bin/xspecfem3D ]; then
  echo "Compiling first all binaries in the root directory..."
  echo

  # compiles executables in root directory
  # using default configuration
  cd ../../../

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

  cp bin/xmeshfem3D           $currendir/bin/xmeshfem3D_attenuation
  cp bin/xspecfem3D           $currentdir/bin/xspecfem3D_attenuation
  cp bin/xcreate_movie_AVS_DX $currentdir/bin/xcreate_movie_AVS_DX_attenuation
  cp bin/xcombine_vol_data $currentdir/bin/xcombine_vol_data
  cd $currentdir/

  # only in case of static compilation
  if [ 1 == 0 ]; then
    # no attenuation simulations
    echo; echo
    cp -v DATA/Par_file_NOISE_1_noattenuation DATA/Par_file

    # compiles in root directory
    cd ../../../
    cp -v $dir/DATA/Par_file DATA/
    rm -rf OUTPUT_FILES/*
    make clean
    make all
    cp bin/xmeshfem3D           $dir/bin/xmeshfem3D_noattenuation
    cp bin/xspecfem3D           $dir/bin/xspecfem3D_noattenuation
    cp bin/xcreate_movie_AVS_DX $dir/bin/xcreate_movie_AVS_DX_noattenuation
    cd $currentdir/
  else
    # no static compilation, can use same binaries
    cd bin/
    ln -s xmeshfem3D_attenuation xmeshfem3D_noattenuation
    ln -s xspecfem3D_attenuation xspecfem3D_noattenuation
    ln -s xcreate_movie_AVS_DX_attenuation xcreate_movie_AVS_DX_noattenuation
    cd ../
  fi
fi

echo
echo "done directory: $dir"

# submit job-script
if [ "$do_submit" == "1" ]; then
  echo ""
  echo "submitting pbs-job:"
  sleep 2
  qsub $script
else
  echo ""
  echo "to submit pbs-job use:"
  echo "> qsub $script"
fi
echo
echo
