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
  rm -rf OUTPUT_FILES/
}

# test example
cd $dir

# default setup
# limit number of time steps
sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.5:" DATA/Par_file
# shortens output interval to avoid timeouts
sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file

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
my_test

echo
echo "all good"
echo `date`
echo
