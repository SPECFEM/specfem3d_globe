#=====================================================================
#
#                       S p e c f e m 3 D  G l o b e
#                       ----------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#                and CNRS / University of Marseille, France
#                 (there are currently many more authors!)
# (c) Princeton University and CNRS / University of Marseille, April 2014
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================

#######################################
##
## Green function extraction: lib/libgf3d.a and bin/xgf3d
##
## This is the only *serial* link target in the tree apart from
## xcreate_header_file, whose recipe (src/create_header_file/rules.mk:117)
## is the template followed here: plain ${FCCOMPILE_CHECK}, no $(MPILIBS),
## and an object list that contains neither $O/parallel.sharedmpi.o nor
## $O/exit_mpi.shared.o.
##
## Note $(LDFLAGS) is referenced directly rather than through $(MPILIBS).
## The HDF5 libraries are appended to LDFLAGS (Makefile.in:528) and LDFLAGS
## only reaches ordinary links via `MPILIBS += $(LDFLAGS) @LIBS@`
## (Makefile.in:59), which also drags in @MPILIBS@.
##
#######################################

## the subset that is free of HDF5, MPI and specfem_par, and therefore
## builds and links from a plain ./configure. tests/gf3d/ builds exactly
## this list through the `gf3d_kernels` target below, so that the tests do
## not have to hard-code object file names.
##
## Keeping this list honest is a maintained contract, not a convention:
## every runner in tests/gf3d/ asserts `nm ./bin/<test> | grep ' U .*mpi_'`
## comes back empty.
gf3d_KERNEL_OBJECTS = \
	$O/gf_par.gf3d.o \
	$O/gf_shared_params.gf3d.o \
	$O/gf_mpi_stubs.gf3d.o \
	$O/gf_dirlist.gf3d_cc.o \
	$(EMPTY_MACRO)

## the parts that read the database, and therefore need HDF5
gf3d_HDF5_OBJECTS = \
	$O/gf_hdf5_read.gf3d.o \
	$O/gf_database.gf3d.o \
	$(EMPTY_MACRO)

## library contents (everything except the program itself)
gf3d_OBJECTS = \
	$(gf3d_KERNEL_OBJECTS) \
	$(gf3d_HDF5_OBJECTS) \
	$(EMPTY_MACRO)

## the driver, linked against the library rather than archived into it
gf3d_PROGRAM_OBJECTS = \
	$O/gf3d_main.gf3d.o \
	$(EMPTY_MACRO)

## These files come from the shared directory.
##
## Kept deliberately short: the list grows stage by stage, and a short list
## makes an accidental MPI dependency obvious the moment it appears.
##
## model_topo_bathy.shared.o is here for get_topo_bathy(). It cannot be had
## on its own — src/shared/rules.mk:72 compiles the whole file into one
## object, and a static archive pulls in whole objects — so it arrives
## together with model_topo_bathy_broadcast(), whose references to wtime(),
## bcast_all_i() and exit_MPI() are satisfied by $O/gf_mpi_stubs.gf3d.o.
gf3d_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/binary_c_io.cc.o \
	$O/flush_system.shared.o \
	$O/model_topo_bathy.shared.o \
	$(EMPTY_MACRO)

gf3d_MODULES = \
	$(FC_MODDIR)/gf_par.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_hdf5_read.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_shared_params.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_database.$(FC_MODEXT) \
	$(EMPTY_MACRO)

#######################################

## The database format is HDF5-only, so the products are only built when
## the tree was configured --with-hdf5. Note this gates what is *built by
## default*, not what *can* be built: the pattern rules below are
## unconditional so that tests/gf3d/ can compile the kernels from a plain
## ./configure.
ifeq ($(HDF5), yes)

gf3d_TARGETS = \
	$L/libgf3d.a \
	$E/xgf3d \
	$(EMPTY_MACRO)

else

gf3d_TARGETS = \
	$(EMPTY_MACRO)

## `make gf3d` would otherwise succeed silently with nothing to do
gf3d: gf3d_no_hdf5

gf3d_no_hdf5:
	@echo ""
	@echo "skipping gf3d: the Green function extraction library requires HDF5"
	@echo "  re-run configure with --with-hdf5 to build lib/libgf3d.a and bin/xgf3d"
	@echo ""

.PHONY: gf3d_no_hdf5

endif

## builds only the HDF5-free kernels, for the unit test suite in tests/gf3d/
gf3d_kernels: $(gf3d_KERNEL_OBJECTS)

.PHONY: gf3d_kernels

#######################################

####
#### rules for the library and the executable
####

$L/libgf3d.a: $(gf3d_OBJECTS) $(gf3d_SHARED_OBJECTS)
	@-mkdir -p $L
	-rm -f $@
	$(AR) $(ARFLAGS) $@ $(gf3d_OBJECTS) $(gf3d_SHARED_OBJECTS)
	$(RANLIB) $@

## serial link: no $(MPILIBS), and $(LDFLAGS) referenced directly for HDF5
$E/xgf3d: $(gf3d_PROGRAM_OBJECTS) $L/libgf3d.a
	${FCCOMPILE_CHECK} -o $@ $(gf3d_PROGRAM_OBJECTS) $L/libgf3d.a $(LDFLAGS)

#######################################

## compilation directories
S := ${S_TOP}/src/gf3d
$(gf3d_OBJECTS): S = ${S_TOP}/src/gf3d
$(gf3d_PROGRAM_OBJECTS): S = ${S_TOP}/src/gf3d

####
#### rule for each .o file below
####

## intra-directory dependencies
$O/gf_hdf5_read.gf3d.o: $O/gf_par.gf3d.o
$O/gf_shared_params.gf3d.o: $O/gf_par.gf3d.o
$O/gf_database.gf3d.o: $O/gf_par.gf3d.o $O/gf_hdf5_read.gf3d.o $O/gf_shared_params.gf3d.o
$O/gf3d_main.gf3d.o: $O/gf_par.gf3d.o $O/gf_database.gf3d.o

## unique object suffix: every rules.mk writes into the same $O, so the
## pattern rules of different subdirectories must not collide
$O/%.gf3d.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gf3d.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gf3d_cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $@ $<
