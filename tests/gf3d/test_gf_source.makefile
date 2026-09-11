# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_source

## compilation directories
O := ./obj

# Kernel objects only -- no HDF5 archive exists after
# 0.configure.default_make.sh.
#
# get_cmt.solver.o and get_force.solver.o are the point of this test: it
# compares the values route in gf_source against the solver's own readers,
# so both have to be here. They build through the $O/%.solver.o pattern rule
# that `include Makefile` brings in (src/specfem3D/rules.mk), and neither
# uses MPI or HDF5; get_cmt calls julian_day(), hence calendar.shared.o.
#
# gf_shared_params.gf3d.o installs the globals both readers depend on:
# NUMBER_OF_SIMULTANEOUS_RUNS and NOISE_TOMOGRAPHY have no initialiser in
# shared_par.f90, and R_PLANET/RHOAV want their Earth defaults left standing.
OBJECTS = \
	$O/gf_par.gf3d.o \
	$O/gf_shared_params.gf3d.o \
	$O/gf_source.gf3d.o \
	$O/shared_par.shared_module.o \
	$O/calendar.shared.o \
	$O/get_cmt.solver.o \
	$O/get_force.solver.o \
	$(EMPTY_MACRO)

# Serial link: ${FCCOMPILE_CHECK}, no $(MPILIBS), no $(LDFLAGS).
test_gf_source: $(OBJECTS)
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_source \
		gf_manufactured.f90 test_gf_source.f90 -I./obj $(OBJECTS)
