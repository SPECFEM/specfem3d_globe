#!/bin/bash
###################################################
#
# Runs test_gf_locate against a Green function database from gfdb_env.bash.
#
# Two kinds of check live here. The kd-tree ownership guard and the geometry
# of whatever element is found need only a valid database, so they run
# against the synthetic fixture. The comparison of the geographic chain
# against the solver's own Cartesian position needs the forward run that
# produced that database, and there is no substitute for it -- the library is
# the only other implementation of that chain, so checking it against itself
# would prove nothing. Those values therefore come from a reference file,
# named by $GF3D_TEST_REFERENCE, and the comparison is made only when one is
# supplied. REF_DATA/reference_regional.txt is the one for the shipped
# regional example.
#
# Override the example directory with $GF3D_TEST_EXAMPLE.
#
###################################################

testdir=`pwd`

# executable
var=test_gf_locate

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

# resolves a database and the source files: $GF3D_TEST_GFDB, or a fixture
. ./gfdb_env.bash
if [ $? -ne 0 ]; then
  echo "skipped: no database and no fixture could be built" >> $testdir/results.log
  echo "skipped: no database and no fixture could be built"
  exit 0
fi

###################################################
#
# the request, and the solver's answer to it
#
# The position comes from the CMTSOLUTION in REF_DATA/, which is test data.
# The Cartesian position the solver settled on comes from a reference file,
# because it is a property of one forward run and not of the library: it is
# supplied only when GF3D_TEST_REFERENCE names one. Without it the test runs
# everything except that comparison.
#
###################################################

LAT=`grep -E '^latitude:'  "$CMT" | head -1 | awk '{print $2}'`
LON=`grep -E '^longitude:' "$CMT" | head -1 | awk '{print $2}'`
DEP=`grep -E '^depth:'     "$CMT" | head -1 | awk '{print $2}'`

if [ -z "$LAT" ] || [ -z "$LON" ] || [ -z "$DEP" ]; then
  echo "could not read the source position from $CMT" >> $testdir/results.log
  exit 1
fi

XREF=`read_reference x`
YREF=`read_reference y`
ZREF=`read_reference z`
XIREF=`read_reference xi`
ETAREF=`read_reference eta`
GAMREF=`read_reference gamma`

echo "request:  lat=$LAT lon=$LON depth=$DEP km" >> $testdir/results.log
if [ -n "$XREF" ] && [ -n "$YREF" ] && [ -n "$ZREF" ]; then
  echo "solver:   x=$XREF y=$YREF z=$ZREF" >> $testdir/results.log
else
  echo "solver:   (no reference supplied; position comparison not asserted)" >> $testdir/results.log
  XREF=""; YREF=""; ZREF=""
fi
if [ -n "$XIREF" ] && [ -n "$ETAREF" ] && [ -n "$GAMREF" ]; then
  echo "solver:   xi=$XIREF eta=$ETAREF gamma=$GAMREF" >> $testdir/results.log
else
  XIREF=""; ETAREF=""; GAMREF=""
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
# version node stripped, see 6.test_gf_open.sh
echo "checking that $var links no MPI" >> $testdir/results.log
if nm ./bin/$var | sed 's/@.*//' | grep -i ' U .*mpi_' >> $testdir/results.log 2>&1; then
  echo "$var references MPI, please check..." >> $testdir/results.log
  exit 1
fi

# runs test
echo "run: `date`" >> $testdir/results.log
./bin/$var "$GFDB" "$LAT" "$LON" "$DEP" $XREF $YREF $ZREF $XIREF $ETAREF $GAMREF \
  >> $testdir/results.log 2>$testdir/error.log

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
# cross-check `xgf3d --locate` against the same reference
#
# The tool is the form a user runs; it must not drift from the library the
# test drives.
#
###################################################

if [ -e ./bin/xgf3d ] && [ -n "$XREF" ]; then
  echo "" >> $testdir/results.log
  echo "cross-checking xgf3d --locate against the solver reference" >> $testdir/results.log

  ./bin/xgf3d --locate "$GFDB" "$LAT" "$LON" "$DEP" > $testdir/loc.log 2>$testdir/error.log
  if [[ $? -ne 0 ]] || [[ -s $testdir/error.log ]]; then
    echo "xgf3d --locate failed:" >> $testdir/results.log
    cat $testdir/error.log >> $testdir/results.log
    exit 1
  fi
  rm -f $testdir/error.log

  for pair in "x:$XREF" "y:$YREF" "z:$ZREF"; do
    key="${pair%%:*}"
    ref="${pair##*:}"
    got=`grep -E "^  $key +=" $testdir/loc.log | head -1 | sed 's/.*= *//'`
    if [ -z "$got" ]; then
      echo "  could not read $key from xgf3d --locate" >> $testdir/results.log
      exit 1
    fi
    # the solver prints through sngl(), so the reference itself only carries
    # float32 precision (~3e-8 at these magnitudes); 1e-7 is the bar
    ok=`awk -v a="$ref" -v b="$got" 'BEGIN{ d=a-b; if (d<0) d=-d; print (d <= 1e-7) ? "yes" : "no" }'`
    if [ "$ok" != "yes" ]; then
      echo "  MISMATCH $key: solver says '$ref', xgf3d says '$got'" >> $testdir/results.log
      exit 1
    fi
    echo "  ok $key = $ref (xgf3d: $got)" >> $testdir/results.log
  done

  rm -f $testdir/loc.log
else
  echo "bin/xgf3d not present, skipping the --locate cross-check" >> $testdir/results.log
fi

#cleanup
rm -f bin/$var
rm -f *.mod ./obj/gf_manufactured.mod
# done
echo "successfully tested: `date`" >> $testdir/results.log
