#!/bin/bash
###################################################
#
# Stations over OpenMP threads. Configures a second, private copy of the
# library with --enable-openmp in ./omp_build -- the build that every other
# runner here tests stays the default one, without OpenMP -- and runs
# test_gf_cache against it, whose section 9 then compares the threaded
# stations with the serial ones, bit for bit.
#
###################################################

testdir=`pwd`

# executable
var=test_gf_cache

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi
cd $ROOT/
srcdir=`pwd`
cd $testdir/

# title
echo >> $testdir/results.log
echo "test: OpenMP build, $var" >> $testdir/results.log
echo >> $testdir/results.log
echo "directory: `pwd`" >> $testdir/results.log

# needs the HDF5 build from 5.configure.hdf5_make.sh, and HDF5 found the way
# that script finds it
if [ ! -e ./lib/libgf3d.a ]; then
  echo "skipped: no HDF5 build of libgf3d (see 5.configure.hdf5_make.sh)" >> $testdir/results.log
  echo "skipped: no HDF5 build of libgf3d"
  exit 0
fi
if [ -z "${HDF5_INC}" ] || [ -z "${HDF5_LIBS}" ]; then
  h5wrap=""
  for w in h5pfc h5fc; do
    if command -v $w > /dev/null 2>&1; then h5wrap=$w; break; fi
  done
  if [ -n "$h5wrap" ]; then
    h5prefix=`$h5wrap -show 2>/dev/null | tr ' ' '\n' | grep '^-L' | head -1 | sed 's/^-L//'`
    if [ -n "$h5prefix" ] && [ -d "$h5prefix" ]; then
      HDF5_LIBS="-L$h5prefix"
      h5root=`dirname $h5prefix`
      if [ -d "$h5root/include" ]; then HDF5_INC="$h5root/include" ; fi
    fi
  fi
fi
if [ -z "${HDF5_INC}" ] || [ -z "${HDF5_LIBS}" ]; then
  echo "skipped: HDF5 not available" >> $testdir/results.log
  echo "skipped: HDF5 not available"
  exit 0
fi

# resolves a database and the source files: $GF3D_TEST_GFDB, or a fixture
. ./gfdb_env.bash
if [ $? -ne 0 ]; then
  echo "skipped: no database and no fixture could be built" >> $testdir/results.log
  echo "skipped: no database and no fixture could be built"
  exit 0
fi

# the OpenMP build
bdir=$testdir/omp_build
rm -rf $bdir
mkdir -p $bdir
cd $bdir
echo "configuration: $srcdir/configure --with-hdf5 --enable-openmp" >> $testdir/results.log
$srcdir/configure --with-hdf5 --enable-openmp HDF5_INC="${HDF5_INC}" HDF5_LIBS="${HDF5_LIBS}" \
  >> $testdir/results.log 2>&1
if [[ $? -ne 0 ]]; then
  echo "could not configure the OpenMP build, please check..." >> $testdir/results.log
  exit 1
fi
echo "compilation: OpenMP libgf3d" >> $testdir/results.log
make gf3d >> $testdir/results.log 2>&1
if [[ $? -ne 0 ]] || [ ! -e ./lib/libgf3d.a ]; then
  echo "compilation of the OpenMP libgf3d failed, please check..." >> $testdir/results.log
  exit 1
fi

# the threaded loop must be in the library: a parallel region in the
# object that holds it (GOMP_* from gfortran, __kmpc_* from Intel)
echo "checking that gf_seismograms has a parallel region" >> $testdir/results.log
if ! nm ./obj/gf_seismograms.gf3d.o | grep -E -q 'GOMP_parallel|__kmpc_fork_call'; then
  echo "the OpenMP build's gf_seismograms has no parallel region, please check..." >> $testdir/results.log
  exit 1
fi

# the test, against that build
echo "compilation: $var (OpenMP)" >> $testdir/results.log
make -f $testdir/test_gf_openmp.makefile test_gf_cache_omp TESTDIR=$testdir >> $testdir/results.log 2>&1
if [ ! -e ./bin/$var ]; then
  echo "compilation of $var failed, please check..." >> $testdir/results.log
  exit 1
fi
cd $testdir

# runs test; its output on its own as well, for the check below
echo "run: `date`" >> $testdir/results.log
$bdir/bin/$var "$GFDB" "$CMT" "$FORCE" > $bdir/test.out 2>$testdir/error.log
status=$?
cat $bdir/test.out >> $testdir/results.log

# checks exit code
if [[ $status -ne 0 ]]; then
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

# section 9 must have run: this runner exists for it. Read from this run's
# output alone: results.log also holds the build log, in which a compiler
# message can quote the very source line looked for.
if ! grep -q "ok   CMT + 10 partials + onset, threads against one" $bdir/test.out; then
  echo "section 9 did not run: was the test compiled with OpenMP?" >> $testdir/results.log
  exit 1
fi

#cleanup
rm -rf $bdir
# done
echo "successfully tested: `date`" >> $testdir/results.log
