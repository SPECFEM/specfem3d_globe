#!/bin/bash
###################################################

# test name
NAME="xwrite_profile"

# relative location of SPECFEM3D_GLOBE EXAMPLES/ directory for test
EXAMPLES="../../EXAMPLES/global_s362ani_shakemovie"

###################################################

# bash function for checking profile
my_test(){
  echo "testing profile:";
  paste REF_DATA/CARDS_th0042_ph0338 OUTPUT_FILES/CARDS_th0042_ph0338 > tmp.dat;
  # compare radius
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($1 - $10)**2;}}END{print "L2 radius = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare rho
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($2 - $11)**2;}}END{print "L2 rho = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare vpv
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($3 - $12)**2;}}END{print "L2 vpv = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare vsv
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($4 - $13)**2;}}END{print "L2 vsv = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare Qkappa
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($5 - $14)**2;}}END{print "L2 Qkappa = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare Qmu
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($6 - $15)**2;}}END{print "L2 Qmu = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare vph
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($7 - $16)**2;}}END{print "L2 vph = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare vsh
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($8 - $17)**2;}}END{print "L2 vsh = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
  # compare eta
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($9 - $18)**2;}}END{print "L2 eta = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then exit 1; fi
}

testdir=`pwd`

# title
echo >> $testdir/results.log
echo "$NAME in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

#cleanup output
rm -rf ./DATABASES_MPI ./OUTPUT_FILES ./DATA ./tmp.dat
mkdir -p OUTPUT_FILES DATABASES_MPI

# setup
cp -rp $EXAMPLES/DATA .
cd ./DATA
ln -s ../../../DATA/crust1.0
ln -s ../../../DATA/crust2.0
ln -s ../../../DATA/crustmap
ln -s ../../../DATA/epcrust
ln -s ../../../DATA/eucrust-07
ln -s ../../../DATA/full_sphericalharmonic_model

ln -s ../../../DATA/s20rts
ln -s ../../../DATA/s40rts
ln -s ../../../DATA/s362ani
ln -s ../../../DATA/sglobe
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
  echo "run failed, please check..." >> $testdir/results.log
  exit 1
fi

# cleanup
rm -rf ./OUTPUT_FILES ./DATABASES_MPI ./DATA ./tmp.dat

echo "successful run" >> $testdir/results.log

