# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_chunk_read

## compilation directories
O := ./obj
L := ./lib

OBJECTS = \
	$(EMPTY_MACRO)

# Links against the library that 5.configure.hdf5_make.sh just built, as
# test_gf_cache does. The test also writes an HDF5 file of its own, which
# $(LDFLAGS) already provides for.
test_gf_chunk_read:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_chunk_read \
		gf_manufactured.f90 test_gf_chunk_read.f90 $L/libgf3d.a $(LDFLAGS)
