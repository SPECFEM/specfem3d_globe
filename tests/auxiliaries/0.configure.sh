#!/bin/bash
###################################################

# test name
NAME="0.configure"

# configuration parameters
CONF_PARAM="--enable-vectorization"

###################################################


testdir=`pwd`

# sets source directory
cd $ROOT/
srcdir=`pwd`

cd $testdir/

# title
echo >> $testdir/results.log
echo "$NAME in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

#cleanup
rm -rf config.log config.status
rm -rf ./bin ./obj ./setup ./OUTPUT_FILES ./DATA

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

echo "" >> $testdir/results.log
echo "successful configuration" >> $testdir/results.log

