#!/bin/bash
###################################################
#
# Runs test_gf_locate against a built example Green function database, using
# the solver's own OUTPUT_FILES/output_solver.txt as the reference.
#
# That file is the only oracle this step has: gf_cross_validate.py never
# converts lat/lon/depth to Cartesian itself, it reads x, y, z straight out
# of the solver log. So the reference values are extracted here with grep and
# awk and passed to the test as arguments.
#
# Needs, for one example:
#   <example>/GFDB/mesh_info.h5
#   <example>/validation_data/CMTSOLUTION
#   <example>/forward_cmt/OUTPUT_FILES/output_solver.txt
#
# The last of those is *not* committed -- forward_cmt/ is gitignored -- so
# this test skips cleanly on a fresh checkout and in CI, like its siblings.
#
# When both example databases are present the second is passed as well, which
# exercises the kd-tree ownership guard in src/gf3d/gf_locate.F90.
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

# locates an example that has a database *and* a forward CMT run to compare against
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
  echo "  needs <example>/GFDB, <example>/validation_data/CMTSOLUTION and" >> $testdir/results.log
  echo "  <example>/forward_cmt/OUTPUT_FILES/output_solver.txt" >> $testdir/results.log
  echo "  build one with EXAMPLES/green_function_database/*/Snakefile" >> $testdir/results.log
  echo "skipped: no example with both a database and a forward CMT run"
  exit 0
fi

GFDB="$EX/GFDB"
CMT="$EX/validation_data/CMTSOLUTION"
SOLVER="$EX/forward_cmt/OUTPUT_FILES/output_solver.txt"

echo "example:  $EX" >> $testdir/results.log

# a second database, for the kd-tree ownership guard
GFDB2=""
for cand in "$srcdir/EXAMPLES/green_function_database/regional/GFDB" \
            "$srcdir/EXAMPLES/green_function_database/global/GFDB"; do
  if [ -e "$cand/mesh_info.h5" ] && [ "$cand" != "$GFDB" ]; then GFDB2="$cand"; break; fi
done

###################################################
#
# the reference values
#
###################################################

# the requested source position, from the CMTSOLUTION the forward run used
LAT=`grep -E '^latitude:'  "$CMT" | head -1 | awk '{print $2}'`
LON=`grep -E '^longitude:' "$CMT" | head -1 | awk '{print $2}'`
DEP=`grep -E '^depth:'     "$CMT" | head -1 | awk '{print $2}'`

# the Cartesian position the solver settled on:
#    at (x,y,z)                  =   0.248383403  -0.944765568  -9.87713933E-02
XYZ=`grep 'at (x,y,z)' "$SOLVER" | head -1 | sed 's/.*= *//'`
XREF=`echo $XYZ | awk '{print $1}'`
YREF=`echo $XYZ | awk '{print $2}'`
ZREF=`echo $XYZ | awk '{print $3}'`

if [ -z "$LAT" ] || [ -z "$LON" ] || [ -z "$DEP" ] || \
   [ -z "$XREF" ] || [ -z "$YREF" ] || [ -z "$ZREF" ]; then
  echo "could not extract the reference values:" >> $testdir/results.log
  echo "  lat='$LAT' lon='$LON' depth='$DEP'" >> $testdir/results.log
  echo "  x='$XREF' y='$YREF' z='$ZREF'" >> $testdir/results.log
  exit 1
fi

echo "request:  lat=$LAT lon=$LON depth=$DEP km" >> $testdir/results.log
echo "solver:   x=$XREF y=$YREF z=$ZREF" >> $testdir/results.log
if [ -n "$GFDB2" ]; then
  echo "second:   $GFDB2" >> $testdir/results.log
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
./bin/$var "$GFDB" "$LAT" "$LON" "$DEP" "$XREF" "$YREF" "$ZREF" "$GFDB2" \
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

if [ -e ./bin/xgf3d ]; then
  echo "" >> $testdir/results.log
  echo "cross-checking xgf3d --locate against output_solver.txt" >> $testdir/results.log

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
