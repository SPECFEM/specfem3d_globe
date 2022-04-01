#!/bin/bash
###################################################

# executable
var=xspecfem3D

# configuration parameters
CONF_PARAM="--enable-vectorization --enable-openmp"

###################################################
testdir=`pwd`
me=`basename "$0"`

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi

# sets source directory
cd $ROOT/
srcdir=`pwd`

cd $testdir/

# title
echo >> $testdir/results.log
echo "$me in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

#cleanup
rm -rf config.log config.status
rm -rf ./bin ./obj ./setup ./OUTPUT_FILES ./DATA/Par_file

# configuration
# (out-of-source compilation)
echo "configuration: $srcdir/configure ${CONF_PARAM}" >> $testdir/results.log
$srcdir/configure ${CONF_PARAM} >> $testdir/results.log 2>&1

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "configuration failed, please check..." >> $testdir/results.log
  exit 1
fi

# we need to output to console output, otherwise tests will fail by timeout in travis
sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h >> $testdir/results.log

# after configuration, we should have a default Par_file in DATA/
# checks exit code
if [ ! -e DATA/Par_file ]; then
  echo "Error: DATA/Par_file not found" >> $testdir/results.log
  echo "configuration failed, please check..." >> $testdir/results.log
  exit 1
fi

# setup model for tests
# overimposes a fixed model setup for testing
# (check also with test_save.f90 in ../meshfem3D test directory)
sed -i "s:^NCHUNKS .*:NCHUNKS = 1:" DATA/Par_file
sed -i "s:^NPROC_XI .*:NPROC_XI  = 2:" DATA/Par_file
sed -i "s:^NPROC_ETA .*:NPROC_ETA = 2:" DATA/Par_file

sed -i "s:^NEX_XI .*:NEX_XI    = 48:" DATA/Par_file
sed -i "s:^NEX_ETA .*:NEX_ETA   = 48:" DATA/Par_file

sed -i "s:^MODEL .*:MODEL = transversely_isotropic_prem_plus_3D_crust_1.0:" DATA/Par_file

sed -i "s:^ANGULAR_WIDTH_XI_IN_DEGREES .*:ANGULAR_WIDTH_XI_IN_DEGREES  = 90.d0:" DATA/Par_file
sed -i "s:^ANGULAR_WIDTH_ETA_IN_DEGREES .*:ANGULAR_WIDTH_ETA_IN_DEGREES = 90.d0:" DATA/Par_file
sed -i "s:^CENTER_LATITUDE_IN_DEGREES .*:CENTER_LATITUDE_IN_DEGREES   = 90.d0:" DATA/Par_file
sed -i "s:^CENTER_LONGITUDE_IN_DEGREES .*:CENTER_LONGITUDE_IN_DEGREES  = 0.d0:" DATA/Par_file
sed -i "s:^GAMMA_ROTATION_AZIMUTH .*:GAMMA_ROTATION_AZIMUTH       = 0.d0:" DATA/Par_file

sed -i "s:^OCEANS .*:OCEANS      = .true.:" DATA/Par_file
sed -i "s:^TOPOGRAPHY .*:TOPOGRAPHY  = .true.:" DATA/Par_file
sed -i "s:^ELLIPTICITY .*:ELLIPTICITY = .true.:" DATA/Par_file
sed -i "s:^ATTENUATION .*:ATTENUATION = .true.:" DATA/Par_file
sed -i "s:^GRAVITY .*:GRAVITY     = .true.:" DATA/Par_file
sed -i "s:^ROTATION .*:ROTATION    = .true.:" DATA/Par_file

sed -i "s:ABSORBING_CONDITIONS .*:ABSORBING_CONDITIONS = .false.:" DATA/Par_file

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "compilation failed, please check..." >> $testdir/results.log
  exit 1
fi

echo "" >> $testdir/results.log
echo "successful configuration" >> $testdir/results.log

# single compilation
echo "compilation: $var" >> $testdir/results.log
make clean >> $testdir/results.log 2>&1
make -j 4 $var >> $testdir/results.log 2>&1

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "compilation failed, please check..." >> $testdir/results.log
  exit 1
fi

# checks binary
if [ ! -e bin/$var ]; then
  echo "compilation of $var failed, please check..." >> $testdir/results.log
  exit 1
else
  echo "binary exists: $var" >> $testdir/results.log
fi
echo "" >> $testdir/results.log

#cleanup
rm -rf ./bin/*

echo "successful compilation" >> $testdir/results.log


