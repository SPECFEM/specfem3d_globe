#!/bin/bash
testdir=`pwd`

# executable
var=test_save

# title
echo >> $testdir/results.log
echo "test: $var" >> $testdir/results.log
echo >> $testdir/results.log

echo "directory: `pwd`" >> $testdir/results.log

# clean
mkdir -p bin
rm -f ./bin/$var

mkdir -p OUTPUT_FILES
rm -f OUTPUT_FILES/*

mkdir -p DATABASES_MPI
rm -f DATABASES_MPI/*

# single compilation
echo "compilation: $var" >> $testdir/results.log

make -f $var.makefile $var >> $testdir/results.log 2>&1

echo "" >> $testdir/results.log

# check
if [ ! -e ./bin/$var ]; then
  echo "compilation of $var failed, please check..." >> $testdir/results.log
  exit 1
fi

# runs test
echo "run: `date`" >> $testdir/results.log
OMP_NUM_THREADS=2 mpirun -np 4 ./bin/$var >> $testdir/results.log 2>$testdir/error.log

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
# done
echo "successfully tested: `date`" >> $testdir/results.log
