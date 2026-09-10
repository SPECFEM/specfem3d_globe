#!/bin/bash
###################################################
#
# Runs test_gf_source against the HDF5-free build from
# 0.configure.default_make.sh.
#
# No database, no HDF5, no MPI: this is a tier-1 test of the two ways to
# build a source -- from a file, through the solver's own readers, and from
# values, as the C facade does -- which must agree. Runs on every commit in
# CI. See gf3df_integration_plan/testing.md.
#
###################################################

testdir=`pwd`

# executable
var=test_gf_source

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

# needs the kernel objects from 0.configure.default_make.sh
if [ ! -e ./obj/gf_source.gf3d.o ]; then
  echo "skipped: no gf3d kernel objects (see 0.configure.default_make.sh)" >> $testdir/results.log
  echo "skipped: no gf3d kernel objects"
  exit 0
fi

# clean
mkdir -p bin OUTPUT_FILES
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

# The kernels are meant to be free of MPI and of HDF5. An accidental
# `use specfem_par` or `use hdf5` in one of them shows up here long before
# lib/libgf3d.so would notice.
echo "checking that $var links no MPI" >> $testdir/results.log
if nm ./bin/$var | grep -i ' U .*mpi_' >> $testdir/results.log 2>&1; then
  echo "$var references MPI, please check..." >> $testdir/results.log
  exit 1
fi

echo "checking that $var links no HDF5" >> $testdir/results.log
if nm ./bin/$var | grep -i ' U .*h5[a-z]*_' >> $testdir/results.log 2>&1; then
  echo "$var references HDF5, please check..." >> $testdir/results.log
  exit 1
fi

# runs test
echo "run: `date`" >> $testdir/results.log
./bin/$var >> $testdir/results.log 2>$testdir/error.log

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
rm -f ./OUTPUT_FILES/test_gf_source.CMTSOLUTION ./OUTPUT_FILES/test_gf_source.FORCESOLUTION
# done
echo "successfully tested: `date`" >> $testdir/results.log
