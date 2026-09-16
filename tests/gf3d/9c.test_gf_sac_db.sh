#!/bin/bash
###################################################
#
# The SAC output against the forward run's SAC headers and against the
# extraction's own ASCII output, on a built example database.
#
# Runs xgf3d --seis with --format all --partials 1 and hands the directory
# to utils/green_function/gf_sac_check.py (obspy): every file readable,
# the reference time / event / station / component headers equal to the
# forward run's, b within one stored sample of the forward run's -t0, and
# the SAC data equal to the ASCII columns in single precision.
#
# Most of what is asserted is the SAC output against the extraction's own
# ASCII output -- the same numbers through two writers, which is how the
# single-precision storage floor and the -ftz subnormal flush were caught --
# and that needs only a database, so it runs against the fixture. The header
# comparison against the solver's own *.sem.sac additionally needs the
# forward run that produced the database: set $GF3D_TEST_FORWARD_SAC (with a
# matching $GF3D_TEST_GFDB) and it happens too.
#
# Needs a Python with obspy, in $GF3D_PYTHON or on PATH.
#
###################################################

testdir=`pwd`

# executable
var=xgf3d

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi
cd $ROOT/
srcdir=`pwd`
cd $testdir/

# title
echo >> $testdir/results.log
echo "test: sac output ($var --seis --format all --partials 1)" >> $testdir/results.log
echo >> $testdir/results.log
echo "directory: `pwd`" >> $testdir/results.log

# needs the HDF5 build from 5.configure.hdf5_make.sh
if [ ! -e ./bin/$var ]; then
  echo "skipped: no HDF5 build of xgf3d (see 5.configure.hdf5_make.sh)" >> $testdir/results.log
  echo "skipped: no HDF5 build of xgf3d"
  exit 0
fi

# resolves a database and the source files: $GF3D_TEST_GFDB, or a fixture
. ./gfdb_env.bash
if [ $? -ne 0 ]; then
  echo "skipped: no database and no fixture could be built" >> $testdir/results.log
  echo "skipped: no database and no fixture could be built"
  exit 0
fi

# The forward run's own seismograms, if there are any. Most of what
# gf_sac_check.py asserts is the SAC output against the extraction's own
# ASCII -- same numbers, two writers -- and that needs no solver run. Only
# the header comparison does, so $GF3D_TEST_FORWARD_SAC is optional and the
# test runs either way.
FWD="${GF3D_TEST_FORWARD_SAC}"
FWDARG=""
if [ -n "$FWD" ] && ls "$FWD"/*.sem.sac > /dev/null 2>&1; then
  FWDARG="--fwd $FWD"
  echo "forward:  $FWD" >> $testdir/results.log
else
  FWD=""
  echo "forward:  (none; SAC headers not compared against a solver run)" >> $testdir/results.log
fi

# a Python with obspy
PY="${GF3D_PYTHON:-python3}"
if ! "$PY" -c "import obspy, numpy" > /dev/null 2>&1; then
  echo "skipped: no Python with obspy (set GF3D_PYTHON)" >> $testdir/results.log
  echo "skipped: no Python with obspy (set GF3D_PYTHON)"
  exit 0
fi

echo "python:   $PY" >> $testdir/results.log

OUT="$testdir/OUTPUT_FILES/sac_check"

echo "example:  $EX" >> $testdir/results.log
echo "python:   $PY" >> $testdir/results.log

rm -rf "$OUT"
mkdir -p "$OUT"

# extraction, every format, with the moment-tensor partials
echo "run: `date`" >> $testdir/results.log
./bin/$var --seis "$GFDB" "$CMT" "$OUT" --format all --partials 1 >> $testdir/results.log 2>$testdir/error.log
if [[ $? -ne 0 ]]; then
  echo "xgf3d failed"; echo "error log:"; cat $testdir/error.log; echo ""
  exit 1
fi
if [[ -s $testdir/error.log ]]; then
  echo "returned ERROR output:" >> $testdir/results.log
  cat $testdir/error.log >> $testdir/results.log
  exit 1
fi

# the check
"$PY" "$srcdir/utils/green_function/gf_sac_check.py" --dir "$OUT" $FWDARG >> $testdir/results.log 2>$testdir/error.log
if [[ $? -ne 0 ]]; then
  echo "test failed"; echo "error log:"; cat $testdir/error.log; echo ""
  tail -30 $testdir/results.log
  exit 1
fi
rm -f $testdir/error.log

#cleanup
rm -rf "$OUT"
# done
echo "successfully tested: `date`" >> $testdir/results.log
