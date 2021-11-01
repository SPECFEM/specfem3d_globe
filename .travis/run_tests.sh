#!/bin/bash

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

# checks if anything to do
echo "run checks: $RUN_CHECKS"
if [ "$RUN_CHECKS" == "0" ]; then
  echo "  no run checks required, exiting..."
  exit 0
else
  echo "  run checks required, start testing..."
fi
echo

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
  5) dir=EXAMPLES/regional_sgloberani/ ;;
  6) dir=EXAMPLES/global_full_sphericalharmonic_model/ ;;
  7) dir=EXAMPLES/regional_s40rts/ ;;
  8) dir=EXAMPLES/regional_small_adjoint/ ;;
  9) dir=EXAMPLES/mars_regional/ ;;
  10) dir=EXAMPLES/moon_global/ ;;
  11) dir=EXAMPLES/regional_Greece_small_LDDRK/ ;;
  12) dir=EXAMPLES/regional_Greece_noise_small/ ;;
  *) dir=EXAMPLES/regional_Greece_small/ ;;
esac


# info
#echo $TRAVIS_BUILD_DIR
echo $WORKDIR
echo `date`
echo
echo "**********************************************************"
echo
echo "run test: TESTID=${TESTID} TESTDIR=${TESTDIR} TESTCOV=${TESTCOV} "
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
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.999 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
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
  # shortens output interval to avoid timeouts
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file

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

  # regional noise
  if [ "$TESTID" == "21" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.1:" DATA/Par_file
    sed -i "s:2999:199:g" run_this_example.kernel.sh
    # uses kernel script by default
    cp -v run_this_example.kernel.sh run_this_example.sh
  fi

  # coverage run
  if [ "$TESTCOV" == "1" ]; then
    sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  fi

  # default script
  ./run_this_example.sh

  # checks script return code
  if [[ $? -ne 0 ]]; then
    # cleanup
    rm -rf OUTPUT_FILES* DATABASES_MPI*
    exit 1
  fi

  # simulation done
  echo
  echo "simulation done: `pwd`"
  echo `date`
  echo

  # seismogram comparison
  if [ "$TESTCOV" == "0" ] && [ ! "$TESTID" == "7" ] && [ ! "$TESTID" == "8" ] && [ ! "$TESTID" == "21" ]; then
    my_test
  fi
fi

#checks
if [[ $? -ne 0 ]]; then exit 1; fi

# simulation done
echo
echo "test done: `pwd`"
echo `date`
echo

echo -en 'travis_fold:end:tests\\r'
echo

# code coverage: https://codecov.io/gh/geodynamics/specfem3d/
# additional runs for coverage
#
# note: log becomes too long, trying to fold each test output
cd $WORKDIR

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.regular-kernel\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing regular_kernel
  ##
  echo "##################################################################"
  echo "cd EXAMPLES/regular_kernel/"
  echo
  cd EXAMPLES/regular_kernel/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.regular-kernel\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.regional-LDDRK\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing regional LDDRK
  ##
  echo "##################################################################"
  echo "EXAMPLES/regional_Greece_small_LDDRK/"
  echo
  cd EXAMPLES/regional_Greece_small_LDDRK/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.regional-LDDRK\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.regional-noise\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing regional noise
  ##
  echo "##################################################################"
  echo "EXAMPLES/regional_Greece_noise_small/"
  echo
  cd EXAMPLES/regional_Greece_noise_small/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.1:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  sed -i "s:2999:199:g" run_this_example.kernel.sh
  ./run_this_example.kernel.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.regional-noise\\r'

## second testcov run built with vectorization
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.global-small\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing global_small
  ##
  echo "##################################################################"
  echo "cd EXAMPLES/global_small/"
  echo
  cd EXAMPLES/global_small/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.global-small\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.regional-s40rts\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing regional s40rts
  ##
  echo "##################################################################"
  echo "EXAMPLES/regional_s40rts/"
  echo
  cd EXAMPLES/regional_s40rts/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.regional-s40rts\\r'


#################################################################
##
## tested by github actions
##
#################################################################

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.regional-sgloberani\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing regional_sgloberani
  ##
  echo "##################################################################"
  echo "cd EXAMPLES/regional_sgloberani/"
  echo
  cd EXAMPLES/regional_sgloberani/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.regional-sgloberani\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.point-force\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing point_force
  ##
  echo "##################################################################"
  echo "cd EXAMPLES/point_force/"
  echo
  cd EXAMPLES/point_force/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.point-force\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.moon-global\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing moon global
  ##
  echo "##################################################################"
  echo "EXAMPLES/moon_global/"
  echo
  cd EXAMPLES/moon_global/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.moon-global\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.mars-regional\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing mars regional
  ##
  echo "##################################################################"
  echo "EXAMPLES/mars_regional/"
  echo
  cd EXAMPLES/mars_regional/
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.mars-regional\\r'


# done
echo "all done"
echo `date`
echo
