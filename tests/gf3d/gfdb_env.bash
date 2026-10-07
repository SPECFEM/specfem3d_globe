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
#   GF3D_FIXTURE_CHUNKED  =1 makes the fixture's displacement datasets chunked
#                         as the solver's writer chunks them, rather than
#                         contiguous (see make_fixture_db.f90). The contiguous
#                         default exercises the h5dread_f route of the chunk
#                         reader; 9h builds both variants itself.
#   GF3D_TEST_STRICT      =1 makes a skip a failure: z.strict_no_skips.sh
#                         then fails the run if any runner here wrote a
#                         `skipped:` line. Not read by this file; set by the
#                         one CI job that installs HDF5, so that this tier
#                         cannot silently go unrun again.
#
# Exports: GFDB, CMT, FORCE, REFERENCE, FIXTURE_DIR (empty unless we made one).
# Defines gfdb_make_fixture [chunked], for a runner that needs another fixture,
# and gfdb_extra_variants (below).
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
FIXTURE_DIRS=()
GFDB=""

# gfdb_make_fixture [chunked] -> FIXTURE_MADE, the GFDB of a new fixture.
# Every fixture made goes on the FIXTURE_DIRS array, which the EXIT trap
# below removes; an array, so that a $TMPDIR with a space in it cannot split
# an rm -rf argument in two.
# A variable rather than stdout, because $(...) would run this in a subshell
# and the directory would never reach FIXTURE_DIRS.
gfdb_make_fixture() {
  local d
  FIXTURE_MADE=""
  # Small (~7 MB) and affine, so gf_open, gf_locate_source, the anchor guard
  # and extraction all run; only the solver comparisons sit out. Rebuilt when
  # its source is newer, so that a runner started on its own does not test
  # yesterday's fixture.
  if [ ! -e ./bin/make_fixture_db ] || [ make_fixture_db.f90 -nt ./bin/make_fixture_db ]; then
    make -f fixture.makefile make_fixture_db >> $testdir/results.log 2>&1
  fi
  if [ ! -e ./bin/make_fixture_db ]; then
    echo "could not build make_fixture_db; see above" >> $testdir/results.log
    return 1
  fi
  d=`mktemp -d "${TMPDIR:-/tmp}/gf3d_fixture.XXXXXX"` || return 1
  FIXTURE_DIRS+=("$d")
  if ! ./bin/make_fixture_db "$d" $1 >> $testdir/results.log 2>&1; then
    echo "make_fixture_db failed" >> $testdir/results.log
    return 1
  fi
  FIXTURE_MADE="$d/GFDB"
  return 0
}

# gfdb_extra_variants -> EXTRA_GFDBS, the fixture layouts $GFDB is not.
# For a runner whose main run is on $GFDB and which has one section that must
# also run on both layouts of the fixture (contiguous: the h5dread_f route of
# the chunk reader; chunked as the solver writes: the raw H5Dread_chunk
# route). $GFDB a fixture: the other layout, and failing to build it is an
# error. $GFDB a user database: both layouts, and failing to build them only
# leaves the array short, since that database is what was asked for.
gfdb_extra_variants() {
  EXTRA_GFDBS=()
  if [ -n "$FIXTURE_DIR" ]; then
    if [ "${GF3D_FIXTURE_CHUNKED}" = "1" ]; then
      gfdb_make_fixture || return 1
    else
      gfdb_make_fixture chunked || return 1
    fi
    EXTRA_GFDBS+=("$FIXTURE_MADE")
  else
    if gfdb_make_fixture; then EXTRA_GFDBS+=("$FIXTURE_MADE"); fi
    if gfdb_make_fixture chunked; then EXTRA_GFDBS+=("$FIXTURE_MADE"); fi
  fi
  return 0
}

gfdb_cleanup() {
  local d
  for d in "${FIXTURE_DIRS[@]}"; do rm -rf "$d"; done
  return 0
}

# Installed here rather than left to each runner: a fixture is a few megabytes
# in $TMPDIR and every runner has half a dozen exit paths, so relying on each
# of them to remember is relying on the one that does not. Sourced, so this
# trap belongs to the runner's own shell. Installed before the first fixture
# is made, so that a failed build is cleaned up too.
trap gfdb_cleanup EXIT

if [ -n "${GF3D_TEST_GFDB}" ]; then

  if [ ! -e "${GF3D_TEST_GFDB}/mesh_info.h5" ]; then
    echo "GF3D_TEST_GFDB is set but has no mesh_info.h5: ${GF3D_TEST_GFDB}" >> $testdir/results.log
    return 1
  fi
  GFDB="${GF3D_TEST_GFDB}"

else

  # no database given: make one
  if [ "${GF3D_FIXTURE_CHUNKED}" = "1" ]; then
    gfdb_make_fixture chunked || return 1
  else
    gfdb_make_fixture || return 1
  fi
  GFDB="$FIXTURE_MADE"
  FIXTURE_DIR=`dirname "$GFDB"`
fi

# read_reference <key> -> the value, or empty when there is no reference file.
# '#' comment lines never match because the key is compared against field 1.
read_reference() {
  [ -n "$REFERENCE" ] && [ -e "$REFERENCE" ] || return 0
  awk -v k="$1" '$1 == k { print $2; exit }' "$REFERENCE"
}

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
