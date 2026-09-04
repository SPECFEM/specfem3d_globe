#!/bin/bash
###################################################
#
# Builds the HDF5-free part of src/gf3d/ from a plain ./configure.
#
# This is the check that keeps the "these files stay HDF5-free, MPI-free and
# specfem_par-free" contract honest. The manufactured-solution tests that
# arrive with the extraction kernels drive those kernels with synthetic
# arrays -- no database, no HDF5, no MPI -- and that only works if the
# kernels never acquire such a dependency. CI runs `make tests` after a bare
# ./configure (see .github/workflows/CI.yml, Test 0), so this runs on every
# commit with nothing installed.
#
# The object list comes from the `gf3d_kernels` target in
# src/gf3d/rules.mk rather than being spelled out here, so that adding a
# kernel is a one-line change in one place.
#
###################################################

# target
var=gf3d_kernels

# configuration parameters
CONF_PARAM=""

###################################################
testdir=`pwd`
me=`basename "$0"`

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi

# sets source directory
cd $ROOT/
srcdir=`pwd`

cd $testdir/

# title
echo >> $testdir/results.log
echo "$me in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

#cleanup
rm -rf config.log config.status
rm -rf ./bin ./obj ./lib ./setup ./OUTPUT_FILES ./DATA/Par_file

# configuration
# (out-of-source compilation, deliberately without --with-hdf5)
echo "configuration: $srcdir/configure ${CONF_PARAM}" >> $testdir/results.log
$srcdir/configure ${CONF_PARAM} >> $testdir/results.log 2>&1

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "configuration failed, please check..." >> $testdir/results.log
  exit 1
fi

echo "" >> $testdir/results.log
echo "successful configuration" >> $testdir/results.log

# `make gf3d` must skip cleanly rather than fail when HDF5 is off
echo "checking that make gf3d skips without HDF5" >> $testdir/results.log
make gf3d >> $testdir/results.log 2>&1
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "make gf3d should skip cleanly without HDF5, but it failed..." >> $testdir/results.log
  exit 1
fi

# no library or executable may have been produced
if [ -e ./lib/libgf3d.a ] || [ -e ./bin/xgf3d ]; then
  echo "make gf3d built something without HDF5, please check..." >> $testdir/results.log
  exit 1
fi

# compilation of the HDF5-free kernels
echo "compilation: $var" >> $testdir/results.log
make $var >> $testdir/results.log 2>&1

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "compilation failed, please check..." >> $testdir/results.log
  exit 1
fi

# checks objects exist
nobj=`ls -1 ./obj/gf_*.o 2>/dev/null | wc -l`
if [[ $nobj -eq 0 ]]; then
  echo "compilation of $var produced no objects, please check..." >> $testdir/results.log
  exit 1
fi
echo "built $nobj kernel objects" >> $testdir/results.log

# MPI-freeness: an accidental `use specfem_par` in a kernel shows up here long
# before Stage 9's shared object would notice
echo "checking that no kernel object references MPI" >> $testdir/results.log
if nm ./obj/gf_*.o | grep -i ' U .*mpi_' >> $testdir/results.log 2>&1; then
  echo "a gf3d kernel object references MPI, please check..." >> $testdir/results.log
  exit 1
fi

# same for HDF5
echo "checking that no kernel object references HDF5" >> $testdir/results.log
if nm ./obj/gf_*.o | grep -i ' U .*h5[a-z]*_' >> $testdir/results.log 2>&1; then
  echo "a gf3d kernel object references HDF5, please check..." >> $testdir/results.log
  exit 1
fi

echo "" >> $testdir/results.log
echo "successful compilation" >> $testdir/results.log
