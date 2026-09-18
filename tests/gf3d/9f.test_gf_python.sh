#!/bin/bash
###################################################
#
# Runs test_gf_python: the ctypes package against xgf3d's own output.
#
# Needs the shared library from 5.configure.hdf5_make.sh, a database, and a
# Python with numpy. All three are absent on a fresh checkout and in Test 0,
# so all three are skips rather than failures. Test 19 supplies all three
# and runs with GF3D_TEST_STRICT=1, where a skip fails.
#
# Note the Python step is judged by its exit code alone, not by whether it
# wrote to stderr: numpy and obspy warn there routinely, and the test itself
# provokes a RuntimeWarning or two on purpose. The compile-time checks in
# the other runners keep the stricter rule.
#
###################################################

testdir=`pwd`

# executable
var=test_gf_python

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

# needs the shared library and the executable to compare against
if [ ! -e ./lib/libgf3d.so ]; then
  echo "skipped: no shared build of libgf3d (see 5.configure.hdf5_make.sh)" >> $testdir/results.log
  echo "skipped: no shared build of libgf3d"
  exit 0
fi
if [ ! -e ./bin/xgf3d ]; then
  echo "skipped: no xgf3d to compare against" >> $testdir/results.log
  echo "skipped: no xgf3d"
  exit 0
fi

# a Python with numpy: GF3D_PYTHON, else whatever is on PATH
PY="${GF3D_PYTHON:-`command -v python3`}"
if [ -z "$PY" ] || ! "$PY" -c "import numpy" > /dev/null 2>&1; then
  echo "skipped: no Python with numpy (set GF3D_PYTHON)" >> $testdir/results.log
  echo "skipped: no Python with numpy"
  exit 0
fi
echo "python: $PY" >> $testdir/results.log

# resolves a database and the source files: $GF3D_TEST_GFDB, or a fixture
. ./gfdb_env.bash
if [ $? -ne 0 ]; then
  echo "skipped: no database and no fixture could be built" >> $testdir/results.log
  echo "skipped: no database and no fixture could be built"
  exit 0
fi

# runs test
echo "run: `date`" >> $testdir/results.log
PYTHONPATH="$srcdir/utils/green_function" \
GF3D_LIB="$testdir/lib/libgf3d.so" \
  "$PY" "$testdir/$var.py" "$testdir/bin/xgf3d" "$GFDB" "$CMT" "$FORCE" \
  >> $testdir/results.log 2>$testdir/error.log

# checks exit code
if [[ $? -ne 0 ]]; then
  echo "test failed"; echo "error log:"; cat $testdir/error.log; echo ""
  echo "results:"; tail -n 40 $testdir/results.log
  exit 1
fi

# stderr is recorded but does not fail the test: see the note at the top
if [[ -s $testdir/error.log ]]; then
  echo "stderr (not a failure):" >> $testdir/results.log
  cat $testdir/error.log >> $testdir/results.log
fi
rm -f $testdir/error.log

#cleanup
rm -rf $srcdir/utils/green_function/gf3d/__pycache__
# done
echo "successfully tested: `date`" >> $testdir/results.log
