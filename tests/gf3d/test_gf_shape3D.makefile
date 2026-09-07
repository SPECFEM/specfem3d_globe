# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_shape3D

## compilation directories
O := ./obj

# The kernel objects, plus recompute_jacobian.shared.o -- which is the point
# of this test: gf_shape3D.F90 is a fork of it, and the original serves as a
# bit-for-bit oracle.
#
# Note this links the *objects*, not $L/libgf3d.a. 0.configure.default_make.sh
# configures without --with-hdf5, so no archive exists; the whole reason
# src/gf3d/rules.mk keeps a gf3d_KERNEL_OBJECTS list separate from the HDF5
# half is that this test can be built and run in CI with nothing installed.
OBJECTS = \
	$O/gf_par.gf3d.o \
	$O/gf_shape3D.gf3d.o \
	$O/hex_nodes.shared.o \
	$O/shared_par.shared_module.o \
	$O/recompute_jacobian.shared.o \
	$(EMPTY_MACRO)

# Serial link: ${FCCOMPILE_CHECK}, not ${MPIFCCOMPILE_CHECK}, and no
# $(MPILIBS) -- see src/gf3d/rules.mk for why. No $(LDFLAGS) either: these
# kernels are HDF5-free, which is what the runner's `nm` check asserts.
test_gf_shape3D:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_shape3D \
		gf_manufactured.f90 test_gf_shape3D.f90 -I./obj $(OBJECTS)
