#!/bin/bash
###################################################
#
# Runs test_gf_open against a built example Green function database, and
# diffs `xgf3d --info` against the values h5dump reads out of mesh_info.h5.
#
# Skips cleanly when there is nothing to run against. Both are the normal
# case in CI: 5.configure.hdf5_make.sh skips without HDF5, and the example
# databases are gitignored (300 MB to 1.7 GB), so neither bin/xgf3d nor a
# database will exist there.
#
# Database search order:
#   1. $GF3D_TEST_GFDB
#   2. EXAMPLES/green_function_database/regional/GFDB
#   3. EXAMPLES/green_function_database/global/GFDB
#
###################################################

testdir=`pwd`

# executable
var=test_gf_open

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

# locates a database
GFDB=""
for cand in "${GF3D_TEST_GFDB}" \
            "$srcdir/EXAMPLES/green_function_database/regional/GFDB" \
            "$srcdir/EXAMPLES/green_function_database/global/GFDB"; do
  if [ -n "$cand" ] && [ -e "$cand/mesh_info.h5" ]; then GFDB="$cand"; break; fi
done

if [ -z "$GFDB" ]; then
  echo "skipped: no example Green function database found" >> $testdir/results.log
  echo "  build one with EXAMPLES/green_function_database/*/Snakefile," >> $testdir/results.log
  echo "  or point GF3D_TEST_GFDB at one" >> $testdir/results.log
  echo "skipped: no example Green function database found"
  exit 0
fi

echo "database: $GFDB" >> $testdir/results.log

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
if nm ./bin/$var | grep -i ' U .*mpi_' >> $testdir/results.log 2>&1; then
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
# cross-check `xgf3d --info` against h5dump
#
# The Stage 1 acceptance criterion is that the inspector prints what is
# actually in the file, so the oracle is the file itself, read by a tool
# that shares no code with ours.
#
###################################################

if command -v h5dump > /dev/null 2>&1; then
  echo "" >> $testdir/results.log
  echo "cross-checking xgf3d --info against h5dump" >> $testdir/results.log

  ./bin/xgf3d --info "$GFDB" > $testdir/info.log 2>$testdir/error.log
  if [[ $? -ne 0 ]] || [[ -s $testdir/error.log ]]; then
    echo "xgf3d --info failed:" >> $testdir/results.log
    cat $testdir/error.log >> $testdir/results.log
    exit 1
  fi
  rm -f $testdir/error.log

  # compares the integer attributes, which h5dump prints unambiguously
  for pair in "nstep:nstep" "nt_subsampled:nt_subsampled" \
              "subsample_step:subsample_step" "ngllx:ngll" \
              "NX_BATHY:NX_BATHY" "NY_BATHY:NY_BATHY" "nspl:nspl"; do
    attr="${pair%%:*}"
    key="${pair##*:}"

    ref=`h5dump -a "/$attr" "$GFDB/mesh_info.h5" 2>/dev/null | \
         grep -A1 '^ *DATA {' | tail -1 | tr -d ' ' | sed 's/^(0)://'`
    [ -z "$ref" ] && continue

    got=`grep -E "^ +$key +=" $testdir/info.log | head -1 | \
         sed 's/.*= *//' | awk '{print $1}'`

    if [ "$ref" != "$got" ]; then
      echo "  MISMATCH $attr: h5dump says '$ref', xgf3d says '$got'" >> $testdir/results.log
      exit 1
    fi
    echo "  ok $attr = $ref" >> $testdir/results.log
  done

  # element count must equal the manifest line count minus its header
  if [ -e "$GFDB/manifest.csv" ]; then
    ref=$(( `wc -l < "$GFDB/manifest.csv"` - 1 ))
    got=`grep -E "^ +nelem +=" $testdir/info.log | sed 's/.*= *//'`
    if [ "$ref" != "$got" ]; then
      echo "  MISMATCH nelem: manifest.csv has $ref entries, xgf3d says $got" >> $testdir/results.log
      exit 1
    fi
    echo "  ok nelem = $ref (== manifest.csv lines - 1)" >> $testdir/results.log
  fi

  # station count must equal the number of files in stations/
  ref=`ls -1 "$GFDB/stations/"*.h5 2>/dev/null | wc -l | tr -d ' '`
  got=`grep -E "^ +nstations +=" $testdir/info.log | sed 's/.*= *//'`
  if [ "$ref" != "$got" ]; then
    echo "  MISMATCH nstations: stations/ holds $ref files, xgf3d says $got" >> $testdir/results.log
    exit 1
  fi
  echo "  ok nstations = $ref (== files in stations/)" >> $testdir/results.log

  rm -f $testdir/info.log
else
  echo "h5dump not available, skipping the metadata cross-check" >> $testdir/results.log
fi

#cleanup
rm -f bin/$var
# done
echo "successfully tested: `date`" >> $testdir/results.log
