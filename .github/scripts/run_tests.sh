#!/bin/bash
#
# runs a test example case
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`
dir=${TESTDIR}

# info
echo "work directory: $WORKDIR"
echo `date`
echo
echo "**********************************************************"
echo
echo "test directory: $dir"
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
}

# test example
cd $dir

# default setup
# limit number of time steps
sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.5:" DATA/Par_file
# shortens output interval to avoid timeouts
sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file

# specific example setups
if [ "${TESTDIR}" == "EXAMPLES/global_small" ]; then
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.1:" DATA/Par_file
fi

# debug
if [ "${DEBUG}" == "true" ]; then
  # limit for debugging
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
fi

# full gravity
if [ "${FULL_GRAVITY}" == "true" ]; then
  # set NSTEP for short check
  echo "NSTEP = 2" >> DATA/Par_file
  # turns on full gravity
  sed -i "s:^FULL_GRAVITY .*:FULL_GRAVITY = .true.:" DATA/Par_file
  # switch to PETSc Poisson solver
  if [ "${PETSC}" == "true" ]; then
    sed -i "s:^POISSON_SOLVER .*:POISSON_SOLVER = 1:" DATA/Par_file
  fi
fi

# default script
./run_this_example.sh

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# simulation done
echo
echo "simulation done: `pwd`"
echo `date`
echo

# seismogram comparison
if [ "${DEBUG}" == "true" ] || [ "${FULL_GRAVITY}" == "true" ]; then
  # no comparisons
  :     # do nothing
else
  my_test
fi

# cleanup
rm -rf OUTPUT_FILES* DATABASES_MPI*

echo
echo "all good"
echo `date`
echo
