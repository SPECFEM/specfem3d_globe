#!/bin/bash
###################################################
#
# Reconfigures this test directory with --with-hdf5 and builds
# lib/libgf3d.a and bin/xgf3d, for the tests that need a real database.
#
# Skips cleanly when HDF5 is not available. That is the normal case in CI:
# `make tests` runs after a bare ./configure with nothing installed
# (.github/workflows/CI.yml, Test 0), and the example databases are
# gitignored besides -- the global one is 300 MB to 1.7 GB.
#
# HDF5 is located from, in order:
#   1. HDF5_INC and HDF5_LIBS in the environment -- the same names
#      .github/scripts/run_build.sh uses
#   2. h5fc / h5pfc on PATH, whose install prefix is queried
#
###################################################

# executable
var=xgf3d

###################################################
testdir=`pwd`
me=`basename "$0"`

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi

cd $ROOT/
srcdir=`pwd`
cd $testdir/

# title
echo >> $testdir/results.log
echo "$me in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

# locates HDF5
if [ -z "${HDF5_INC}" ] || [ -z "${HDF5_LIBS}" ]; then
  h5wrap=""
  for w in h5pfc h5fc; do
    if command -v $w > /dev/null 2>&1; then h5wrap=$w; break; fi
  done
  if [ -n "$h5wrap" ]; then
    h5prefix=`$h5wrap -show 2>/dev/null | tr ' ' '\n' | grep '^-L' | head -1 | sed 's/^-L//'`
    if [ -n "$h5prefix" ] && [ -d "$h5prefix" ]; then
      HDF5_LIBS="-L$h5prefix"
      # include directory sits next to the library directory in every layout
      # we have seen (lib/, lib64/)
      h5root=`dirname $h5prefix`
      if [ -d "$h5root/include" ]; then HDF5_INC="$h5root/include" ; fi
    fi
  fi
fi

if [ -z "${HDF5_INC}" ] || [ -z "${HDF5_LIBS}" ]; then
  echo "skipped: HDF5 not available" >> $testdir/results.log
  echo "  set HDF5_INC and HDF5_LIBS, or put h5fc on PATH, to run the" >> $testdir/results.log
  echo "  database-backed gf3d tests" >> $testdir/results.log
  echo "skipped: HDF5 not available"
  exit 0
fi

echo "HDF5_INC  = ${HDF5_INC}" >> $testdir/results.log
echo "HDF5_LIBS = ${HDF5_LIBS}" >> $testdir/results.log

#cleanup
rm -rf config.log config.status
rm -rf ./bin ./obj ./lib ./setup ./OUTPUT_FILES ./DATA/Par_file

# configuration
echo "configuration: $srcdir/configure --with-hdf5" >> $testdir/results.log
$srcdir/configure --with-hdf5 HDF5_INC="${HDF5_INC}" HDF5_LIBS="${HDF5_LIBS}" \
  >> $testdir/results.log 2>&1

if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "configuration with HDF5 failed, skipping the database-backed tests" >> $testdir/results.log
  echo "skipped: could not configure with HDF5"
  exit 0
fi

echo "" >> $testdir/results.log
echo "successful configuration" >> $testdir/results.log

# compilation
echo "compilation: $var" >> $testdir/results.log
make gf3d >> $testdir/results.log 2>&1

if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "compilation failed, please check..." >> $testdir/results.log
  exit 1
fi

if [ ! -e ./bin/$var ] || [ ! -e ./lib/libgf3d.a ]; then
  echo "compilation of $var failed, please check..." >> $testdir/results.log
  exit 1
fi
echo "binary exists: $var" >> $testdir/results.log
echo "library exists: lib/libgf3d.a" >> $testdir/results.log

# Stage 9's products: the shared object, the C header and the public module.
# A missing one is a build regression, not a reason to skip: `make gf3d`
# just succeeded, so it promised all four.
for f in ./lib/libgf3d.so ./include/gf3d.h; do
  if [ ! -e "$f" ]; then
    echo "make gf3d did not produce $f, please check..." >> $testdir/results.log
    exit 1
  fi
  echo "exists: $f" >> $testdir/results.log
done

# the module file's name follows the compiler's own case convention
if [ ! -e ./include/gf3d.mod ] && [ ! -e ./include/GF3D.mod ]; then
  echo "make gf3d did not install the public module into include/" >> $testdir/results.log
  exit 1
fi
echo "exists: include/gf3d.mod" >> $testdir/results.log

# the C entry points must actually be exported from the shared object
echo "checking the exported C symbols" >> $testdir/results.log
if command -v nm > /dev/null 2>&1; then
  for sym in gf3d_open gf3d_close gf3d_get_info gf3d_get_station gf3d_locate \
             gf3d_get_plan gf3d_seismograms gf3d_partials gf3d_last_error gf3d_sizeof; do
    if ! nm -D --defined-only ./lib/libgf3d.so 2>/dev/null | grep -q " T $sym$"; then
      echo "libgf3d.so does not export $sym, please check..." >> $testdir/results.log
      exit 1
    fi
  done
  echo "all C entry points are exported" >> $testdir/results.log
fi

# The point of baking an rpath into the shared object is that ctypes can
# load it with no LD_LIBRARY_PATH set. Check that, rather than assuming it.
if command -v ldd > /dev/null 2>&1; then
  echo "checking that libgf3d.so resolves HDF5 without LD_LIBRARY_PATH" >> $testdir/results.log
  env -u LD_LIBRARY_PATH ldd ./lib/libgf3d.so >> $testdir/results.log 2>&1
  if env -u LD_LIBRARY_PATH ldd ./lib/libgf3d.so 2>/dev/null | grep -i "hdf5.*not found" > /dev/null; then
    echo "libgf3d.so cannot find HDF5 without LD_LIBRARY_PATH, please check..." >> $testdir/results.log
    exit 1
  fi
fi

# the library itself must reference no MPI symbol. HDF5 may of course be a
# parallel build and drag libmpi in transitively; that is HDF5's dependency,
# not ours, so this checks our own objects rather than ldd of the binary or
# of the shared object.
echo "checking that libgf3d.a references no MPI" >> $testdir/results.log
if nm ./obj/gf_*.o ./obj/gf3d_main*.o ./obj/gf3d_capi*.o | grep -i ' U .*mpi_' >> $testdir/results.log 2>&1; then
  echo "a gf3d object references MPI, please check..." >> $testdir/results.log
  exit 1
fi

# runs the trivial mode, which needs no database
./bin/$var --version >> $testdir/results.log 2>$testdir/error.log
if [[ $? -ne 0 ]]; then
  echo "xgf3d --version failed" >> $testdir/results.log
  cat $testdir/error.log >> $testdir/results.log
  exit 1
fi
if [[ -s $testdir/error.log ]]; then
  echo "xgf3d --version wrote to stderr:" >> $testdir/results.log
  cat $testdir/error.log >> $testdir/results.log
  exit 1
fi
rm -f $testdir/error.log

echo "" >> $testdir/results.log
echo "successful compilation" >> $testdir/results.log
