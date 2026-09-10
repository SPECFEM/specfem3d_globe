# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_capi

## compilation directories
O := ./obj
L := ./lib
I := ./include

# The only test in the tree that compiles C, because it is the only one
# testing a C ABI. It sees exactly what a downstream C caller sees:
# include/gf3d.h and lib/libgf3d.so, and nothing from ./obj.
#
# Linked against the *shared* library on purpose -- test_gf3d_ext covers the
# archive -- with an rpath, because the runner executes ./bin/test_gf_capi
# with no LD_LIBRARY_PATH set.
#
# Linked with ${CC}, not the Fortran driver: ifort and ifx add their own
# for_main.o to every link, which defines main() and calls MAIN__, so a C
# main collides with it ("multiple definition of `main'"). Nothing is lost
# by using the C driver, because libgf3d.so records the Fortran runtime and
# HDF5 as DT_NEEDED and the dynamic linker follows that -- which is also
# exactly how ctypes loads it, so this link exercises the same property the
# Python package depends on.
test_gf_capi:
	${CC} -c $(CPPFLAGS) $(CFLAGS) -I$I -o $O/test_gf_capi.o test_gf_capi.c
	${CC} -o ./bin/test_gf_capi $O/test_gf_capi.o \
		$L/libgf3d.so -Wl,-rpath,$(abspath $L) -lm
