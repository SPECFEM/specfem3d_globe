#!/bin/bash
###################################################
#
# Runs test_gf_capi: the C round trip over lib/libgf3d.so and
# include/gf3d.h, and the assertion that no path through the facade can end
# the process.
#
# Needs the HDF5 build from 5.configure.hdf5_make.sh and one of the example
# databases, so it skips cleanly on a fresh checkout and in CI, exactly as
# the other database-backed runners here do.
#
###################################################

testdir=`pwd`

# executable
var=test_gf_capi

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi
cd $ROOT/
srcdir=`pwd`
cd $testdir/

# title
echo >> $testdir/results.log
echo "test: $var" >> $testdir/results.log
echo >> $testdir/results.log
echo "directory: `pwd`" >> $testdir/results.log

# needs the shared library and the installed header
if [ ! -e ./lib/libgf3d.so ] || [ ! -e ./include/gf3d.h ]; then
  echo "skipped: no shared build of libgf3d (see 5.configure.hdf5_make.sh)" >> $testdir/results.log
  echo "skipped: no shared build of libgf3d"
  exit 0
fi

# locates an example database with a validation source
EX=""
for cand in "${GF3D_TEST_EXAMPLE}" \
            "$srcdir/EXAMPLES/green_function_database/regional" \
            "$srcdir/EXAMPLES/green_function_database/global"; do
  [ -z "$cand" ] && continue
  if [ -e "$cand/GFDB/mesh_info.h5" ] && [ -e "$cand/validation_data/CMTSOLUTION" ]; then
    EX="$cand"; break
  fi
done

if [ -z "$EX" ]; then
  echo "skipped: no built example database (set GF3D_TEST_EXAMPLE)" >> $testdir/results.log
  echo "skipped: no built example database"
  exit 0
fi

GFDB="$EX/GFDB"
CMT="$EX/validation_data/CMTSOLUTION"
echo "example: $EX" >> $testdir/results.log

# clean
mkdir -p bin obj
rm -f ./bin/$var ./obj/$var.o

# single compilation
echo "compilation: $var" >> $testdir/results.log
make -f $var.makefile $var >> $testdir/results.log 2>&1
echo "" >> $testdir/results.log

# check
if [ ! -e ./bin/$var ]; then
  echo "compilation of $var failed, please check..." >> $testdir/results.log
  exit 1
fi

# runs test
echo "run: `date`" >> $testdir/results.log
./bin/$var "$GFDB" "$CMT" >> $testdir/results.log 2>$testdir/error.log

# checks exit code
if [[ $? -ne 0 ]]; then
  echo "test failed"; echo "error log:"; cat $testdir/error.log; echo ""
  echo "results:"; tail -n 30 $testdir/results.log
  exit 1
fi

# checks error output
if [[ -s $testdir/error.log ]]; then
  echo "returned ERROR output:" >> $testdir/results.log
  cat $testdir/error.log >> $testdir/results.log
  exit 1
fi
rm -f $testdir/error.log

#cleanup
rm -f ./bin/$var ./obj/$var.o
# done
echo "successfully tested: `date`" >> $testdir/results.log
