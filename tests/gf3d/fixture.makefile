# includes default Makefile from previous configuration
include Makefile

# test target
default: make_fixture_db

## compilation directories
O := ./obj
L := ./lib

# The generator writes HDF5 directly and links nothing from src/gf3d/: it is
# the *writer* side of the format, deliberately independent of the reader it
# exists to feed. $(LDFLAGS) carries the HDF5 libraries that
# 5.configure.hdf5_make.sh configured, which is also why this target only
# works after that script has run -- and it is only ever wanted then.
make_fixture_db:
	mkdir -p ./bin
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/make_fixture_db make_fixture_db.f90 $(LDFLAGS)
