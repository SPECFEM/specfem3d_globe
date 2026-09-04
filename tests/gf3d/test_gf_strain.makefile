# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_strain

## compilation directories
O := ./obj

# Kernel objects only -- no HDF5 archive exists after
# 0.configure.default_make.sh, which is the point of keeping
# gf3d_KERNEL_OBJECTS separate in src/gf3d/rules.mk.
OBJECTS = \
	$O/gf_par.gf3d.o \
	$O/gf_shape3D.gf3d.o \
	$O/gf_interp.gf3d.o \
	$O/gf_strain.gf3d.o \
	$O/gf_moment.gf3d.o \
	$O/gll_library.shared.o \
	$O/lagrange_poly.shared.o \
	$O/hex_nodes.shared.o \
	$O/shared_par.shared_module.o \
	$(EMPTY_MACRO)

# Serial link: ${FCCOMPILE_CHECK}, no $(MPILIBS), no $(LDFLAGS).
test_gf_strain:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_strain \
		gf_manufactured.f90 test_gf_strain.f90 -I./obj $(OBJECTS)
