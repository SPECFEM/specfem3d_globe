#
# Resolves the Green function database the database-backed tests run against.
#
# Sourced, not executed:  . ./gfdb_env.bash   then test the return status.
#
# Named .bash and not .sh deliberately. tests/run_tests.sh:124 executes every
# ./*.sh in a test directory as a test, so a sourced helper called .sh is run
# as one -- with $testdir unset and a `return` at top level, which is an
# error, and the whole gf3d directory then fails. The guard below is the
# second line of defence in case that convention ever changes.
#
# This exists so that nothing in tests/gf3d/ knows about EXAMPLES/. The
# example databases are gitignored products of a solver run -- 786 MB and
# 905 MB -- so a test that reaches into that directory is a test that runs on
# one developer's machine and nowhere else. Instead:
#
#   GF3D_TEST_GFDB        a database to use. Unset: build a synthetic fixture
#                         (make_fixture_db) into a temporary directory, which
#                         is what makes this tier run on a bare checkout.
#   GF3D_TEST_REFERENCE   solver reference values for that database, in the
#                         format of REF_DATA/reference_regional.txt. Unset,
#                         the assertions that compare against a solver run are
#                         reported as not exercised rather than skipped
#                         wholesale -- the rest of each test still runs.
#
#                         The fixture deliberately supplies none. Its geometry
#                         is chosen by this suite, so an "expected" position
#                         computed here would be this suite's own arithmetic
#                         checked against itself. testing.md's coverage table
#                         puts geographic->Cartesian in the "pinned only by
#                         the forward run" column and that stays true: the
#                         solver is the only oracle for it.
#   GF3D_TEST_FORWARD_SAC a forward run's *.sem.sac directory (9c only; no
#                         fixture can stand in for the solver's own output).
#   GF3D_PYTHON           interpreter for the Python-dependent runners. 9f
#                         needs numpy, 9c also obspy, scipy and h5py:
#                         `uv sync --group test --project
#                         utils/green_function` builds one at
#                         utils/green_function/.venv/bin/python.
#   GF3D_TEST_STRICT      =1 makes a skip a failure: z.strict_no_skips.sh
#                         then fails the run if any runner here wrote a
#                         `skipped:` line. Not read by this file; set by the
#                         one CI job that installs HDF5, so that this tier
#                         cannot silently go unrun again.
#
# Exports: GFDB, CMT, FORCE, REFERENCE, FIXTURE_DIR (empty unless we made one).
# Returns non-zero when it cannot produce a database, so the caller keeps its
# own "skipped: ...; exit 0" idiom rather than exiting from inside a source.
#

# executed rather than sourced: say so and do nothing, rather than failing
# halfway through with $testdir unset
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
  echo "gfdb_env.bash is sourced by the database test runners, not run on its own."
  exit 0
fi

# the source files are test data, not example data
CMT="$testdir/REF_DATA/CMTSOLUTION"
FORCE="$testdir/REF_DATA/FORCESOLUTION"

REFERENCE="${GF3D_TEST_REFERENCE:-}"
FIXTURE_DIR=""
GFDB=""

if [ -n "${GF3D_TEST_GFDB}" ]; then

  if [ ! -e "${GF3D_TEST_GFDB}/mesh_info.h5" ]; then
    echo "GF3D_TEST_GFDB is set but has no mesh_info.h5: ${GF3D_TEST_GFDB}" >> $testdir/results.log
    return 1
  fi
  GFDB="${GF3D_TEST_GFDB}"

else

  # no database given: make one. Tiny (~100 KB) and affine, so gf_open,
  # gf_locate_source, the anchor guard and extraction all run; only the
  # solver comparisons sit out.
  if [ ! -e ./bin/make_fixture_db ]; then
    make -f fixture.makefile make_fixture_db >> $testdir/results.log 2>&1
  fi
  if [ ! -e ./bin/make_fixture_db ]; then
    echo "could not build make_fixture_db; see above" >> $testdir/results.log
    return 1
  fi

  FIXTURE_DIR=`mktemp -d "${TMPDIR:-/tmp}/gf3d_fixture.XXXXXX"` || return 1
  if ! ./bin/make_fixture_db "$FIXTURE_DIR" >> $testdir/results.log 2>&1; then
    echo "make_fixture_db failed" >> $testdir/results.log
    rm -rf "$FIXTURE_DIR"
    return 1
  fi

  GFDB="$FIXTURE_DIR/GFDB"
fi

# read_reference <key> -> the value, or empty when there is no reference file.
# '#' comment lines never match because the key is compared against field 1.
read_reference() {
  [ -n "$REFERENCE" ] && [ -e "$REFERENCE" ] || return 0
  awk -v k="$1" '$1 == k { print $2; exit }' "$REFERENCE"
}

gfdb_cleanup() {
  [ -n "$FIXTURE_DIR" ] && rm -rf "$FIXTURE_DIR"
  return 0
}

# Installed here rather than left to each runner: a fixture is a few megabytes
# in $TMPDIR and every runner has half a dozen exit paths, so relying on each
# of them to remember is relying on the one that does not. Sourced, so this
# trap belongs to the runner's own shell.
trap gfdb_cleanup EXIT

export GFDB CMT FORCE REFERENCE FIXTURE_DIR

echo "database:  $GFDB" >> $testdir/results.log
if [ -n "$FIXTURE_DIR" ]; then
  echo "           (synthetic fixture; set GF3D_TEST_GFDB to use a real one)" >> $testdir/results.log
fi
if [ -n "$REFERENCE" ]; then
  echo "reference: $REFERENCE" >> $testdir/results.log
else
  echo "reference: (none; solver comparisons will be reported, not asserted)" >> $testdir/results.log
fi

return 0
