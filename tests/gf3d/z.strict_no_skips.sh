#!/bin/bash
###################################################
#
# Fails the run when any other runner in this directory skipped.
#
# Every skip here exits 0, so the whole database tier can go dark and the
# job still passes -- which is exactly how runners 5 through 9f went unrun
# in GitHub CI: a skip was a pass. With GF3D_TEST_STRICT=1 a `skipped:`
# line in results.log is a failure instead.
#
# Opt-in, because most jobs have no HDF5 and must keep skipping: unset, this
# only logs the count, so Test 0, the macOS job and the Intel jobs stay
# green. .github/workflows/CI.yml, Test 19, is the job that sets it.
#
# One check site rather than a guard per runner, and in a runner rather than
# a grep in the YAML, so it is reproducible locally and needs no upkeep when
# runners are added.
#
# Named z.* and not with a number: tests/run_tests.sh runs ./*.sh in glob
# order and `z` sorts after every digit, so this stays last as runners are
# added -- it must be, since it reads what they wrote.
#
# GF3D_TEST_STRICT is documented with this tier's other variables in
# gfdb_env.bash.
#
###################################################

testdir=`pwd`

# Before anything below appends to the same file it is reading. The pattern
# allows leading space: test_gf_python.py's skips are indented.
hits=`grep -nE '^[[:space:]]*skipped:' $testdir/results.log`

nskipped=0
if [ -n "$hits" ]; then
  nskipped=`echo "$hits" | grep -c .`
fi

# title
echo >> $testdir/results.log
echo "test: no runner in this directory skipped" >> $testdir/results.log
echo >> $testdir/results.log
echo "directory: `pwd`" >> $testdir/results.log

if [ "${GF3D_TEST_STRICT}" != "1" ]; then
  echo "GF3D_TEST_STRICT is not set: $nskipped skip(s), not a failure" >> $testdir/results.log
  exit 0
fi

# No line below may itself start with `skipped:`, or a second run of this
# script over the same results.log would report its own output. grep -n
# prefixes every hit with a line number, which is why they are safe to log.
if [ -n "$hits" ]; then
  echo "GF3D_TEST_STRICT=1: $nskipped runner(s) skipped, which this job treats as a failure" >> $testdir/results.log
  echo "$hits" >> $testdir/results.log
  echo "GF3D_TEST_STRICT=1: $nskipped runner(s) skipped:"
  echo "$hits"
  exit 1
fi

# done
echo "GF3D_TEST_STRICT=1: nothing was skipped" >> $testdir/results.log
