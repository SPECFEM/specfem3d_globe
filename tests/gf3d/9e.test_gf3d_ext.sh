#!/bin/bash
###################################################
#
# Runs test_gf3d_ext: an external Fortran program built the way a
# downstream caller builds one, against include/ and lib/libgf3d.a only.
#
# It checks two things at once -- that `use gf3d` alone is enough, and that
# the public Fortran API and the C ABI produce the same numbers.
#
# Needs the HDF5 build from 5.configure.hdf5_make.sh and one of the example
# databases, so it skips cleanly on a fresh checkout and in CI.
#
###################################################

testdir=`pwd`

# executable
var=test_gf3d_ext

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

# needs the archive and the installed module
if [ ! -e ./lib/libgf3d.a ]; then
  echo "skipped: no HDF5 build of libgf3d (see 5.configure.hdf5_make.sh)" >> $testdir/results.log
  echo "skipped: no HDF5 build of libgf3d"
  exit 0
fi
if [ ! -e ./include/gf3d.mod ] && [ ! -e ./include/GF3D.mod ]; then
  echo "skipped: the public module is not installed in include/" >> $testdir/results.log
  echo "skipped: no include/gf3d.mod"
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
FORCE="$EX/validation_data/FORCESOLUTION"
[ -e "$FORCE" ] || FORCE=""
echo "example: $EX" >> $testdir/results.log

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

# runs test
echo "run: `date`" >> $testdir/results.log
./bin/$var "$GFDB" "$CMT" $FORCE >> $testdir/results.log 2>$testdir/error.log

# checks exit code
if [[ $? -ne 0 ]]; then
  echo "test failed"; echo "error log:"; cat $testdir/error.log; echo ""
  echo "results:"; tail -n 30 $testdir/results.log
  exit 1
fi

# checks error output (note: fortran stop returns with a zero-exit code)
if [[ -s $testdir/error.log ]]; then
  echo "returned ERROR output:" >> $testdir/results.log
  cat $testdir/error.log >> $testdir/results.log
  exit 1
fi
rm -f $testdir/error.log

#cleanup
rm -f ./bin/$var
rm -f *.mod
# done
echo "successfully tested: `date`" >> $testdir/results.log
