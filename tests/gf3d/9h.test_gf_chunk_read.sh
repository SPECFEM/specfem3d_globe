#!/bin/bash
###################################################
#
# The direct chunk reader against h5dread_f: the fixture in both layouts
# (contiguous: the h5dread_f route; chunked as the solver writes: the raw
# route), and $GF3D_TEST_GFDB as well when one is given.
#
###################################################

testdir=`pwd`

# executable
var=test_gf_chunk_read

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

# needs the HDF5 build from 5.configure.hdf5_make.sh
if [ ! -e ./lib/libgf3d.a ]; then
  echo "skipped: no HDF5 build of libgf3d (see 5.configure.hdf5_make.sh)" >> $testdir/results.log
  echo "skipped: no HDF5 build of libgf3d"
  exit 0
fi

# resolves $GF3D_TEST_GFDB or a fixture, and defines gfdb_make_fixture
. ./gfdb_env.bash
if [ $? -ne 0 ]; then
  echo "skipped: no database and no fixture could be built" >> $testdir/results.log
  echo "skipped: no database and no fixture could be built"
  exit 0
fi

# the two fixture layouts, whatever gfdb_env.bash chose; its own fixture is
# the contiguous one unless it was asked for the chunked variant
if [ -n "$FIXTURE_DIR" ] && [ "${GF3D_FIXTURE_CHUNKED}" != "1" ]; then
  GFDB_CONTIGUOUS="$GFDB"
elif gfdb_make_fixture; then
  GFDB_CONTIGUOUS="$FIXTURE_MADE"
else
  echo "could not build the contiguous fixture" >> $testdir/results.log
  exit 1
fi
if ! gfdb_make_fixture chunked; then
  echo "could not build the chunked fixture" >> $testdir/results.log
  exit 1
fi
GFDB_CHUNKED="$FIXTURE_MADE"

# clean
mkdir -p bin
rm -f ./bin/$var

# single compilation
echo "compilation: $var" >> $testdir/results.log
make -f $var.makefile $var >> $testdir/results.log 2>&1
echo "" >> $testdir/results.log

# check
if [ ! -e ./bin/$var ]; then
  echo "compilation of $var failed, please check..." >> $testdir/results.log
  exit 1
fi

# a serial library must not have pulled MPI in
# version node stripped, see 6.test_gf_open.sh
echo "checking that $var links no MPI" >> $testdir/results.log
if nm ./bin/$var | sed 's/@.*//' | grep -i ' U .*mpi_' >> $testdir/results.log 2>&1; then
  echo "$var references MPI, please check..." >> $testdir/results.log
  exit 1
fi

# run_one <database> <raw|fallback|any>
run_one() {
  echo "run: $1 (route: $2) `date`" >> $testdir/results.log
  ./bin/$var "$1" "$2" >> $testdir/results.log 2>$testdir/error.log

  # checks exit code
  if [[ $? -ne 0 ]]; then
    echo "test failed"; echo "error log:"; cat $testdir/error.log; echo ""
    exit 1
  fi

  # checks error output (note: fortran stop returns with a zero-exit code)
  if [[ -s $testdir/error.log ]]; then
    echo "returned ERROR output:" >> $testdir/results.log
    cat $testdir/error.log >> $testdir/results.log
    exit 1
  fi
  rm -f $testdir/error.log
}

run_one "$GFDB_CONTIGUOUS" fallback
run_one "$GFDB_CHUNKED" raw
# a real database: the solver chunks it, but a CUSTOM_REAL = 8 one, or one
# written by another tool, may not be, so only the numbers are asserted
if [ -n "${GF3D_TEST_GFDB}" ]; then
  run_one "$GFDB" any
fi

#cleanup
rm -f bin/$var
rm -f *.mod ./obj/gf_manufactured.mod
# done
echo "successfully tested: `date`" >> $testdir/results.log
