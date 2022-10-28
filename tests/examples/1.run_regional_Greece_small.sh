#!/bin/bash
###################################################

# example directory name
var="regional_Greece_small"

# relative location of SPECFEM3D_GLOBE EXAMPLES/ directory for tests directories (e.g. SPECFEM3D_GLOBE/tests/examples)
EXAMPLES="../../EXAMPLES/"

###################################################
testdir=`pwd`
me=`basename "$0"`

# bash function for checking seismogram output with reference solutions
my_test(){
  echo "testing seismograms:";
  ../../utils/compare_seismogram_correlations.py OUTPUT_FILES/ REF_SEIS/;
  ../../utils/compare_seismogram_correlations.py OUTPUT_FILES/ REF_SEIS/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.9 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}';
}


# title
echo >> $testdir/results.log
echo "$me in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

# checks if example directory exists
if [ ! -e $EXAMPLES/$var/run_this_example.sh ]; then
  echo "run-script in directory $EXAMPLES/$var not found, please check..." >> $testdir/results.log
  exit 1
fi

#cleanup output
rm -rf ./DATABASES_MPI ./OUTPUT_FILES ./DATA ./run_this_example.sh ./run_mesher_solver.bash ./REF_SEIS
mkdir -p DATA OUTPUT_FILES

# setup
cp -rp $EXAMPLES/$var/DATA/Par_file ./DATA/
cp -rp $EXAMPLES/$var/DATA/CMTSOLUTION ./DATA/
cp -rp $EXAMPLES/$var/DATA/STATIONS ./DATA/
cp -p $EXAMPLES/$var/run_this_example.sh .
cp -p $EXAMPLES/$var/run_mesher_solver.bash .
ln -s $EXAMPLES/$var/REF_SEIS

# makes sure only 1 x 2 x 2 procs in total are used
#sed -i "s:^NPROC_XI .*:NPROC_XI      = 2:" DATA/Par_file
#sed -i "s:^NPROC_ETA .*:NPROC_ETA    = 2:" DATA/Par_file
# limits length
sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES    = 0.0:" DATA/Par_file

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "setup failed, please check..." >> $testdir/results.log
  exit 1
fi

./run_this_example.sh >> $testdir/results.log 2>&1

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "run failed, please check..." >> $testdir/results.log
  exit 1
fi

# test seismograms
my_test >> $testdir/results.log

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "seismograms failed, please check..." >> $testdir/results.log
  exit 1
fi

# cleanup
rm -rf ./DATABASES_MPI ./OUTPUT_FILES ./DATA ./run_this_example.sh ./run_mesher_solver.bash ./REF_SEIS

echo "successful run" >> $testdir/results.log

