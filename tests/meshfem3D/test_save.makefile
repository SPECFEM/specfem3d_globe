# includes default Makefile from previous configuration
include Makefile

# test target
default: test_save

## compilation directories
O := ./obj

OBJECTS = \
	$(meshfem3D_MESHER_OBJECTS) \
	$(meshfem3D_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

test_save:
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_save test_save.f90 -I./obj $(OBJECTS) $(MPILIBS)

