# includes default Makefile from previous configuration
include Makefile

# test target
default: test_rotate

## compilation directories
O := ./obj

OBJECTS = \
	$O/rotate_tensor.shared.o \
	$(EMPTY_MACRO)

test_rotate:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_rotate test_rotate.f90 -I./obj $(OBJECTS) $(LDFLAGS) $(LIBS)

