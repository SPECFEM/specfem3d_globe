# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf3d_ext

## compilation directories
L := ./lib
I := ./include

# The downstream caller's view, and deliberately nothing more: -I./include
# for gf3d.mod, ./lib/libgf3d.a for the code. No -I./obj, no gf_*.o, no
# gf_manufactured -- if `use gf3d` does not name everything the program
# needs, this does not compile, which is the assertion.
#
# The static library here, where test_gf_capi links the shared one, so that
# both artefacts of `make gf3d` are exercised. $(LDFLAGS) supplies HDF5, as
# it does for xgf3d itself; serial link, no $(MPILIBS).
test_gf3d_ext:
	${FCCOMPILE_CHECK} -I$I -o ./bin/test_gf3d_ext \
		test_gf3d_ext.f90 $L/libgf3d.a $(LDFLAGS)
