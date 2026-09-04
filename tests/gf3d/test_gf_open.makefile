# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_open

## compilation directories
O := ./obj
L := ./lib

OBJECTS = \
	$(EMPTY_MACRO)

# Links against the library that 5.configure.hdf5_make.sh just built, rather
# than recompiling its sources, so the test exercises the same archive a
# downstream caller would get.
#
# Serial link: ${FCCOMPILE_CHECK}, not ${MPIFCCOMPILE_CHECK}, and $(LDFLAGS)
# directly rather than $(MPILIBS) -- see src/gf3d/rules.mk for why.
test_gf_open:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_open test_gf_open.f90 $L/libgf3d.a $(LDFLAGS)
