# includes default Makefile from previous configuration
include Makefile

# test target
default: test_read

## compilation directories
O := ./obj

OBJECTS = \
	$O/specfem3D_par.solverstatic_module.o \
	$O/initialize_simulation.solverstatic.o \
	$O/read_mesh_parameters.solverstatic.o \
	$(EMPTY_MACRO)

OBJECTS += $(specfem3D_SHARED_OBJECTS)

test_read:
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_read test_read.f90 -I./obj $(OBJECTS) $(MPILIBS)

