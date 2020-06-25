#!/bin/bash

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi


###########################################################
# setup
###########################################################
# chooses example directory
case "$TESTDIR" in
  0) dir=./ ;;
  1) dir=EXAMPLES/regional_Greece_small/ ;;
  2) dir=EXAMPLES/global_small/ ;;
  3) dir=EXAMPLES/point_force/ ;;
  4) dir=EXAMPLES/regular_kernel/ ;;
  5) dir=EXAMPLES/global_sglobe/ ;;
  6) dir=EXAMPLES/global_full_sphericalharmonic_model/ ;;
  7) dir=EXAMPLES/regional_s40rts/ ;;
  8) dir=EXAMPLES/regional_small_adjoint/ ;;
  9) dir=EXAMPLES/mars_regional/ ;;
  10) dir=EXAMPLES/moon_global/ ;;
  11) dir=EXAMPLES/regional_Greece_small_LDDRK/ ;;
  *) dir=EXAMPLES/regional_Greece_small/ ;;
esac


# info
echo $TRAVIS_BUILD_DIR
echo $WORKDIR
echo
echo "**********************************************************"
echo
echo "configuration test: TESTID=${TESTID} TESTDIR=${TESTDIR} TESTCOV=${TESTCOV} TESTFLAGS=${TESTFLAGS}"
echo
echo "    test directory: $dir"
echo
echo "**********************************************************"
echo

# bash function for checking seismogram output with reference solutions
my_test(){
  echo "testing seismograms:"
  ln -s $WORKDIR/utils/compare_seismogram_correlations.py
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/
  if [[ $? -ne 0 ]]; then exit 1; fi
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.9 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
  rm -rf OUTPUT_FILES/
}

# bash function for checking seismogram output with reference solutions
my_report(){
  # report example
  if [ -f OUTPUT_FILES/output_mesher.txt ]; then cat OUTPUT_FILES/output_mesher.txt; else exit 1; fi
  if [ -f OUTPUT_FILES/output_solver.txt ]; then cat OUTPUT_FILES/output_solver.txt; else exit 1; fi
}



###########################################################
# configuration & compilation
###########################################################
# configuration
echo 'Configure...' && echo -en 'travis_fold:start:configure\\r'
echo "configuration:"

if [ "$TESTCOV" == "1" ]; then
  echo "configuration: for coverage"
  ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS} FLAGS_CHECK="-fprofile-arcs -ftest-coverage -O0" CFLAGS="-coverage -O0"
else
  if [ "$CUDA" == "true" ]; then
    if [ "$OPENCL" == "true" ]; then
      echo "configuration: for opencl" # uses libOpenCL provided from CUDA package
      ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS} OCL_CPU_FLAGS="-g -Wall -std=c99 -DWITH_MPI" OCL_GPU_FLAGS="-Werror" OCL_INC="${CUDA_HOME}/include" OCL_LIB="${CUDA_HOME}/lib64" OCL_LIBS="-lOpenCL"
    else
      echo "configuration: for cuda"
      ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS} CUDA_LIB="${CUDA_HOME}/lib64" CUDA_INC="${CUDA_HOME}/include" CUDA_FLAGS="-Xcompiler -Wall,-Wno-unused-function,-Wno-unused-const-variable,-Wfatal-errors -g -G"
    fi
  else
    echo "configuration: default"
    ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS}
  fi
fi
# we output to console
sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h

# regional w/ NGLL = 6
if [ "$TESTID" == "8" ]; then
  sed -i "s:NGLLX =.*:NGLLX = 6:" setup/constants.h
fi
echo -en 'travis_fold:end:configure\\r'

# compilation  (only cleaning)
echo 'Build...' && echo -en 'travis_fold:start:build\\r'
echo "compilation:"
make clean
echo -en 'travis_fold:end:build\\r'


###########################################################
# test examples
###########################################################
# testing internal mesher example (short & quick for all configuration)
echo 'Tests...' && echo -en 'travis_fold:start:tests\\r'
# runs test
echo "test directory: $dir"
echo
cd $dir


# testing small example (short & quick for all configurations)
if [ "$TESTID" == "3" ]; then
  # runs default tests
  make tests
else
  # limit number of time steps
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.5:" DATA/Par_file

  # regional w/ debug-checking
  if [ "$TESTID" == "7" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  fi

  # regional w/ NGLL = 6
  if [ "$TESTID" == "8" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  fi

  # global small
  if [ "$TESTID" == "9" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.1:" DATA/Par_file
  fi

  # coverage run
  if [ "$TESTCOV" == "1" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  fi

  # default script
  ./run_this_example.sh

  # seismogram comparison
  if [ "$TESTCOV" == "0" ] && [ ! "$TESTID" == "7" ] && [ ! "$TESTID" == "8" ]; then
    my_test
  fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:tests\\r'


# code coverage: https://codecov.io/gh/geodynamics/specfem3d/
# additional runs for coverage
#
# note: log becomes too long, trying to fold each test output
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.point-force\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing point_force
  ##
  cd EXAMPLES/point_force/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.point-force\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.regular-kernel\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing regular_kernel
  ##
  cd EXAMPLES/regular_kernel/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.regular-kernel\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.global-small\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing global_small
  ##
  cd EXAMPLES/global_small/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.global-small\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.regional-s40rts\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing regional s40rts
  ##
  cd EXAMPLES/regional_s40rts/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.regional-s40rts\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.mars-regional\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing mars regional
  ##
  cd EXAMPLES/mars_regional/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.mars-regional\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.moon-global\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing moon global
  ##
  cd EXAMPLES/moon_global/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.moon-global\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.regional-LDDRK\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing regional LDDRK
  ##
  cd EXAMPLES/regional_Greece_small_LDDRK/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.regional-LDDRK\\r'


# done
echo "done `pwd`"

