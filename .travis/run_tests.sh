#!/bin/bash

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi


###########################################################
# setup
###########################################################
# chooses example directory
case "$TESTMAKE" in
  0) dir=./ ;;
  1) dir=EXAMPLES/regional_Greece_small/ ;;
  2) dir=EXAMPLES/regional_Greece_small/ ;;
  3) dir=EXAMPLES/regional_Greece_small/ ;;
  4) dir=EXAMPLES/global_small/ ;;
  5) dir=EXAMPLES/point_force/ ;;
  6) dir=EXAMPLES/regular_kernel/ ;;
  *) dir=EXAMPLES/regional_Greece_small/ ;;
esac


# info
echo $TRAVIS_BUILD_DIR
echo $WORKDIR
echo
echo "**********************************************************"
echo
echo "configuration test: TESTMAKE=${TESTMAKE} TEST=${TEST} FLAGS=${TESTFLAGS}"
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
  ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TEST} FLAGS_CHECK="-fprofile-arcs -ftest-coverage -O0" CFLAGS="-coverage -O0"
else
  if [ "$CUDA" == "true" ]; then
    echo "configuration: for cuda"
    ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TEST} CUDA_LIB="${CUDA_HOME}/lib64" CUDA_INC="${CUDA_HOME}/include" CUDA_FLAGS="-Xcompiler -Wall,-Wno-unused-function,-Wno-unused-const-variable,-Wfatal-errors -g -G"
  else
    echo "configuration: default"
    ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TEST}
  fi
fi
# we output to console
sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h

# regional w/ NGLL = 6
if [ "$TESTMAKE" == "3" ]; then
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
if [ "$TESTMAKE" == "0" ]; then
  # runs default tests
  make tests
else
  # limit number of time steps
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.5:" DATA/Par_file

  # regional w/ debug-checking
  if [ "$TESTMAKE" == "2" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  fi

  # regional w/ NGLL = 6
  if [ "$TESTMAKE" == "3" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  fi

  # global small
  if [ "$TESTMAKE" == "4" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.1:" DATA/Par_file
  fi

  # coverage run
  if [ "$TESTCOV" == "1" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  fi

  # default script
  ./run_this_example.sh

  # seismogram comparison
  if [ "$TESTCOV" == "0" ] && [ ! "$TESTMAKE" == "2" ] && [ ! "$TESTMAKE" == "3" ]; then
    my_test
  fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:tests\\r'


# code coverage: https://codecov.io/gh/geodynamics/specfem3d/
# additional runs for coverage
# note: log becomes too long, trying to fold each test output

##
## global example
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.point-force\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
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
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
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
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing global_small
  ##
  cd EXAMPLES/global_small/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.global-small\\r'


# done
echo "done `pwd`"

