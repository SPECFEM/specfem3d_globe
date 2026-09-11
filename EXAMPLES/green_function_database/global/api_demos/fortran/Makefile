#=====================================================================
#
# Two demonstrations of calling libgf3d from Fortran.
#
#     make            build both
#     make FC=ifort   build with another compiler (see the note below)
#     make clean
#
# then
#
#     ./gf3d_api_demo
#     ./gf3d_events_demo
#
# The link line is the whole point, so it is worth reading:
#
#     $(FC) demo.f90 -I<repo>/include -L<repo>/lib -lgf3d -Wl,-rpath,<repo>/lib
#
# and that is all. No HDF5 flags, no list of the shared objects the library
# reuses: lib/libgf3d.so records HDF5 as a dependency of its own, and the
# -rpath means the loader can find it without LD_LIBRARY_PATH being set.
#
# Linking lib/libgf3d.a instead works too, and is what tests/gf3d does, but
# then HDF5 becomes the caller's problem -- an archive carries undefined
# symbols and no idea where to satisfy them.
#
# A note on FC. Fortran module files are compiler-specific and often
# version-specific, so <repo>/include/gf3d.mod can only be read by the
# compiler that wrote it -- the one the tree was configured with. If the
# default below is not that compiler you will get "Cannot open module file
# 'gf3d.mod'", which is unambiguous; pass FC= to fix it.
#
#=====================================================================

# the top of the specfem3d_globe tree, from here
SPECFEM_DIR ?= ../../../../..

# `FC ?= gfortran` would not work: make has a built-in default for FC (f77),
# so FC is always already defined and ?= never fires. Overriding only when
# the value is that built-in leaves `make FC=ifort` and an FC exported in
# the environment both working.
ifeq ($(origin FC),default)
  FC = gfortran
endif

# kept separate from FCFLAGS because the HDF5 environment modules on some
# systems export an include path there, and it should not displace this
OPT ?= -O2

INCDIR = $(SPECFEM_DIR)/include
LIBDIR = $(SPECFEM_DIR)/lib

# absolute, because an -rpath is resolved by the loader at run time and has
# no idea what the working directory was at build time
LIBDIR_ABS = $(abspath $(LIBDIR))

GF3D_FLAGS = -I$(INCDIR) -L$(LIBDIR) -lgf3d -Wl,-rpath,$(LIBDIR_ABS)

PROGRAMS = gf3d_api_demo gf3d_events_demo

all: $(PROGRAMS)

$(PROGRAMS): %: %.f90 $(LIBDIR)/libgf3d.so $(INCDIR)/gf3d.mod
	$(FC) $(OPT) $(FCFLAGS) $< -o $@ $(GF3D_FLAGS)

# A clearer failure than "no rule to make target" when the library has not
# been built yet.
$(LIBDIR)/libgf3d.so $(INCDIR)/gf3d.mod:
	@echo ""
	@echo "$@ is missing: libgf3d has not been built."
	@echo "From the top of the tree ($(SPECFEM_DIR)):"
	@echo "    ./configure --with-hdf5 HDF5_INC=<dir> HDF5_LIBS=-L<dir>"
	@echo "    make gf3d"
	@echo ""
	@false

clean:
	rm -f $(PROGRAMS) *.mod *.o

.PHONY: all clean
