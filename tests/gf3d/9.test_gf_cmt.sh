#!/bin/bash
###################################################
#
# The moment-tensor amplitude pin, against the solver's own numbers.
#
# OUTPUT_FILES/output_solver.txt records the scalar moment and moment
# magnitude specfem computed for this very CMTSOLUTION, at full double
# precision. They are extracted here and passed to the test, which
# reproduces them through get_cmt + scaleM.
#
# Needs, for one example:
#   <example>/GFDB/mesh_info.h5
#   <example>/validation_data/CMTSOLUTION
#   <example>/forward_cmt/OUTPUT_FILES/output_solver.txt
#
# forward_cmt/ is gitignored, so this skips cleanly on a fresh checkout and
# in CI, like its siblings. Override the example with $GF3D_TEST_EXAMPLE.
#
###################################################

testdir=`pwd`

# executable
var=test_gf_cmt

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

# locates an example with a database and a forward CMT run
EX=""
for cand in "${GF3D_TEST_EXAMPLE}" \
            "$srcdir/EXAMPLES/green_function_database/regional" \
            "$srcdir/EXAMPLES/green_function_database/global"; do
  [ -z "$cand" ] && continue
  if [ -e "$cand/GFDB/mesh_info.h5" ] && \
     [ -e "$cand/validation_data/CMTSOLUTION" ] && \
     [ -e "$cand/forward_cmt/OUTPUT_FILES/output_solver.txt" ]; then
    EX="$cand"; break
  fi
done

if [ -z "$EX" ]; then
  echo "skipped: no example with both a database and a forward CMT run" >> $testdir/results.log
  echo "skipped: no example with both a database and a forward CMT run"
  exit 0
fi

GFDB="$EX/GFDB"
CMT="$EX/validation_data/CMTSOLUTION"
SOLVER="$EX/forward_cmt/OUTPUT_FILES/output_solver.txt"

echo "example:  $EX" >> $testdir/results.log

# the solver's own source magnitude, at full double precision
M0=`grep 'scalar moment M0'    "$SOLVER" | head -1 | sed 's/.*= *//' | awk '{print $1}'`
MW=`grep 'moment magnitude Mw' "$SOLVER" | head -1 | sed 's/.*= *//' | awk '{print $1}'`
HDUR=`grep 'half duration:'    "$SOLVER" | head -1 | awk '{print $3}'`
TSHIFT=`grep 'time shift:'     "$SOLVER" | head -1 | awk '{print $3}'`

if [ -z "$M0" ] || [ -z "$MW" ] || [ -z "$HDUR" ] || [ -z "$TSHIFT" ]; then
  echo "could not extract the reference values:" >> $testdir/results.log
  echo "  M0='$M0' Mw='$MW' hdur='$HDUR' tshift='$TSHIFT'" >> $testdir/results.log
  exit 1
fi

echo "solver:   M0=$M0 Mw=$MW hdur=$HDUR tshift=$TSHIFT" >> $testdir/results.log

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
./bin/$var "$GFDB" "$CMT" "$M0" "$MW" "$HDUR" "$TSHIFT" \
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

#cleanup
rm -f bin/$var
rm -f *.mod ./obj/gf_manufactured.mod
# done
echo "successfully tested: `date`" >> $testdir/results.log
