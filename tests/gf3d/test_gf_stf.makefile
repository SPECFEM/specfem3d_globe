# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_stf

## compilation directories
O := ./obj

# Kernel objects only -- no HDF5 archive exists after
# 0.configure.default_make.sh, which is the point of keeping
# gf3d_KERNEL_OBJECTS separate in src/gf3d/rules.mk.
OBJECTS = \
	$O/gf_par.gf3d.o \
	$O/gf_stf.gf3d.o \
	$O/shared_par.shared_module.o \
	$(EMPTY_MACRO)

# Two solver routines linked as *oracles*: specfem's own erf
# (comp_source_time_function_heavi is built on it) and its analytic source
# spectrum, which states the 1.628 placement independently of the time
# domain. Neither is part of gf3d_kernels; both are built here through the
# pattern rules that `include Makefile` brings in from src/specfem3D/rules.mk,
# and both are free of MPI and HDF5 (netlib_specfun_erf references only
# `exp`). Production gf_stf uses the F2008 erf intrinsic and never links
# these.
ORACLES = \
	$O/netlib_specfun_erf.solver.o \
	$O/comp_source_spectrum.solver.o \
	$(EMPTY_MACRO)

# Serial link: ${FCCOMPILE_CHECK}, no $(MPILIBS), no $(LDFLAGS).
test_gf_stf: $(ORACLES)
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_stf \
		gf_manufactured.f90 test_gf_stf.f90 -I./obj $(OBJECTS) $(ORACLES)
