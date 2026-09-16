#!/bin/bash
###################################################
#
# Runs test_gf_anchors against a Green function database, and cross-checks it
# against `xgf3d --check-anchors`.
#
# The database comes from gfdb_env.bash. Note the synthetic fixture is affine,
# so its anchor residual is the float32 storage floor rather than a real
# mesh's -- the assertion is the same 1e-6 either way.
#
###################################################

testdir=`pwd`

# executable
var=test_gf_anchors

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
if [ ! -e ./lib/libgf3d.a ] || [ ! -e ./bin/xgf3d ]; then
  echo "skipped: no HDF5 build of libgf3d (see 5.configure.hdf5_make.sh)" >> $testdir/results.log
  echo "skipped: no HDF5 build of libgf3d"
  exit 0
fi

# resolves a database and the source files: $GF3D_TEST_GFDB, or a fixture
. ./gfdb_env.bash
if [ $? -ne 0 ]; then
  echo "skipped: no database and no fixture could be built" >> $testdir/results.log
  echo "skipped: no database and no fixture could be built"
  exit 0
fi

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
echo "checking that $var links no MPI" >> $testdir/results.log
if nm ./bin/$var | grep ' U .*mpi_' >> $testdir/results.log 2>&1; then
  echo "$var references MPI, please check..." >> $testdir/results.log
  exit 1
fi

# runs test
echo "run: `date`" >> $testdir/results.log
./bin/$var "$GFDB" >> $testdir/results.log 2>$testdir/error.log

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

###################################################
#
# cross-check against `xgf3d --check-anchors`
#
# The library routine and the tool must report the same worst residual: the
# tool is the form a user runs, the test is the form CI runs, and they should
# not be able to drift apart.
#
###################################################

echo "" >> $testdir/results.log
echo "cross-checking against xgf3d --check-anchors" >> $testdir/results.log

./bin/xgf3d --check-anchors "$GFDB" > $testdir/anchors.log 2>$testdir/error.log
if [[ $? -ne 0 ]] || [[ -s $testdir/error.log ]]; then
  echo "xgf3d --check-anchors failed:" >> $testdir/results.log
  cat $testdir/error.log >> $testdir/results.log
  cat $testdir/anchors.log >> $testdir/results.log
  exit 1
fi
rm -f $testdir/error.log

ref=`grep -E "^ +worst residual +=" $testdir/anchors.log | sed 's/.*= *//'`
got=`grep -E "^ +worst anchor residual +=" $testdir/results.log | tail -1 | sed 's/.*= *//'`

if [ -z "$ref" ] || [ -z "$got" ]; then
  echo "  could not read the worst residual from both sources" >> $testdir/results.log
  exit 1
fi

# Compares as numbers, not as strings. Both sides print es22.14, but the two
# sweep the database in their own loops, so this is a "the tool and the
# library agree" check at six significant digits, not a bitwise one.
same=`awk -v a="$ref" -v b="$got" 'BEGIN{ d=a-b; if (d<0) d=-d; r=(a<0?-a:a); print (d <= 1e-5*(r>0?r:1)) ? "yes" : "no" }'`
if [ "$same" != "yes" ]; then
  echo "  MISMATCH worst residual: xgf3d says '$ref', test_gf_anchors says '$got'" >> $testdir/results.log
  exit 1
fi
echo "  ok worst residual = $ref (library and tool agree)" >> $testdir/results.log

rm -f $testdir/anchors.log

#cleanup
rm -f bin/$var
rm -f *.mod ./obj/gf_manufactured.mod
# done
echo "successfully tested: `date`" >> $testdir/results.log
