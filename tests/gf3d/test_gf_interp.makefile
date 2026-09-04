# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_interp

## compilation directories
O := ./obj

# The kernel objects only. This links the *objects*, not $L/libgf3d.a,
# because 0.configure.default_make.sh configures without --with-hdf5 and no
# archive exists there -- which is the whole reason src/gf3d/rules.mk keeps
# gf3d_KERNEL_OBJECTS separate from the HDF5 half.
OBJECTS = \
	$O/gf_par.gf3d.o \
	$O/gf_interp.gf3d.o \
	$O/gll_library.shared.o \
	$O/lagrange_poly.shared.o \
	$O/shared_par.shared_module.o \
	$(EMPTY_MACRO)

# Serial link: ${FCCOMPILE_CHECK}, not ${MPIFCCOMPILE_CHECK}, and no
# $(MPILIBS) -- see src/gf3d/rules.mk for why. No $(LDFLAGS) either: these
# kernels are HDF5-free, which is what the runner's `nm` check asserts.
test_gf_interp:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_interp \
		gf_manufactured.f90 test_gf_interp.f90 -I./obj $(OBJECTS)
