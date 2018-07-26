#!/bin/bash
###################################################

# example directory name
NAME="global_s362ani_shakemovie"

# relative location of SPECFEM3D_GLOBE EXAMPLES/ directory for tests directories (e.g. SPECFEM3D_GLOBE/tests/examples)
EXAMPLES="../../EXAMPLES/"

###################################################

# bash function for checking profile
my_test(){
  echo "testing profile:";
  paste CARDS_th0042_ph0338 OUTPUT_FILES/CARDS_th0042_ph0338 > tmp.dat;
  # compare radius
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($1 - $10)**2;}}END{print "radius L2: ",val;if(val>0.01){print "radius failed",val;exit 1;}else{print "radius good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare rho
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($2 - $11)**2;}}END{print "rho L2: ",val;if(val>0.01){print "rho failed",val;exit 1;}else{print "rho good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare vpv
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($3 - $12)**2;}}END{print "vpv L2: ",val;if(val>0.01){print "vpv failed",val;exit 1;}else{print "vpv good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare vsv
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($4 - $13)**2;}}END{print "vsv L2: ",val;if(val>0.01){print "vsv failed",val;exit 1;}else{print "vsv good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
}

testdir=`pwd`

# title
echo >> $testdir/results.log
echo "$NAME in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

#cleanup output
rm -rf ./DATABASES_MPI ./OUTPUT_FILES ./DATA ./run_this_example.sh ./run_mesher_solver.bash ./REF_SEIS
mkdir -p OUTPUT_FILES DATABASES_MPI

# setup
cp -rp $EXAMPLES/$NAME/DATA .
cd ./DATA
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/s362ani
ln -s ../../../DATA/QRFSI12
ln -s ../../../DATA/topo_bathy
cd ../


make xwrite_profile >> $testdir/results.log 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "setup failed, please check..." >> $testdir/results.log
  exit 1
fi

# single profile
./bin/xwrite_profile 1 1 47.5 -22.0 >> $testdir/results.log 2>&1

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "running xwrite_profile failed, please check..." >> $testdir/results.log
  exit 1
fi

# test seismograms
my_test >> $testdir/results.log

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "profile failed, please check..." >> $testdir/results.log
  exit 1
fi

# cleanup
rm -rf ./OUTPUT_FILES ./DATABASES_MPI ./DATA

echo "successful compilation" >> $testdir/results.log

