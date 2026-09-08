# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_sac

## compilation directories
O := ./obj

# Kernel objects only -- no HDF5 archive exists after
# 0.configure.default_make.sh. The SAC writer needs the solver's byte
# writer (binary_c_io.cc.o) and its leap-year rule (calendar.shared.o), both
# part of gf3d_KERNEL_SHARED_OBJECTS in src/gf3d/rules.mk for this reason.
OBJECTS = \
	$O/gf_par.gf3d.o \
	$O/gf_strain.gf3d.o \
	$O/gf_moment.gf3d.o \
	$O/gf_stf.gf3d.o \
	$O/gf_partials.gf3d.o \
	$O/gf_sac.gf3d.o \
	$O/shared_par.shared_module.o \
	$O/binary_c_io.cc.o \
	$O/calendar.shared.o \
	$(EMPTY_MACRO)

# Serial link: ${FCCOMPILE_CHECK}, no $(MPILIBS), no $(LDFLAGS).
test_gf_sac:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_sac \
		gf_manufactured.f90 test_gf_sac.f90 -I./obj $(OBJECTS)
