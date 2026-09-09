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
# Needs, for one example:
#   <example>/GFDB/mesh_info.h5
#   <example>/validation_data/CMTSOLUTION
#   <example>/forward_cmt/OUTPUT_FILES/*.sem.sac
# and a Python with obspy: $GF3D_PYTHON, else the example's own venv
# (EXAMPLES/green_function_database/.venv/bin/python).
#
# forward_cmt/ is gitignored, so this skips cleanly on a fresh checkout and
# in CI, like its siblings. Override the example with $GF3D_TEST_EXAMPLE.
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

# locates an example with a database and a forward CMT run
EX=""
for cand in "${GF3D_TEST_EXAMPLE}" \
            "$srcdir/EXAMPLES/green_function_database/regional" \
            "$srcdir/EXAMPLES/green_function_database/global"; do
  [ -z "$cand" ] && continue
  if [ -e "$cand/GFDB/mesh_info.h5" ] && \
     [ -e "$cand/validation_data/CMTSOLUTION" ] && \
     ls "$cand"/forward_cmt/OUTPUT_FILES/*.sem.sac > /dev/null 2>&1; then
    EX="$cand"; break
  fi
done

if [ -z "$EX" ]; then
  echo "skipped: no example with both a database and forward CMT SAC files" >> $testdir/results.log
  echo "skipped: no example with both a database and forward CMT SAC files"
  exit 0
fi

# a Python with obspy
PY="${GF3D_PYTHON}"
if [ -z "$PY" ]; then PY="$srcdir/EXAMPLES/green_function_database/.venv/bin/python"; fi
if ! "$PY" -c "import obspy, numpy" > /dev/null 2>&1; then
  echo "skipped: no Python with obspy (set GF3D_PYTHON)" >> $testdir/results.log
  echo "skipped: no Python with obspy (set GF3D_PYTHON)"
  exit 0
fi

GFDB="$EX/GFDB"
CMT="$EX/validation_data/CMTSOLUTION"
FWD="$EX/forward_cmt/OUTPUT_FILES"
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
"$PY" "$srcdir/utils/green_function/gf_sac_check.py" --dir "$OUT" --fwd "$FWD" >> $testdir/results.log 2>$testdir/error.log
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
