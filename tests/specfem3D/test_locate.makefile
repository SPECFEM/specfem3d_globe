# includes default Makefile from previous configuration
include Makefile

# test target
default: test_locate

## compilation directories
O := ./obj

OBJECTS = \
	$(specfem3D_SOLVER_OBJECTS) \
	$(specfem3D_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

test_locate:
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_locate test_locate.f90 -I./obj $(OBJECTS) $(MPILIBS)

