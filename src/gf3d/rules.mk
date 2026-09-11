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
	$O/gf_shape3D.gf3d.o \
	$O/gf_geometry.gf3d.o \
	$O/gf_geo_chain.gf3d.o \
	$O/gf_interp.gf3d.o \
	$O/gf_strain.gf3d.o \
	$O/gf_moment.gf3d.o \
	$O/gf_stf.gf3d.o \
	$O/gf_source.gf3d.o \
	$O/gf_partials.gf3d.o \
	$O/gf_sac.gf3d.o \
	$(EMPTY_MACRO)

## Shared objects the kernels call, which are themselves free of HDF5 and
## MPI. They are listed separately from gf3d_SHARED_OBJECTS so that the
## `gf3d_kernels` target can build them for tests/gf3d/, which links the
## manufactured-solution tests from a plain ./configure with no HDF5.
##
## Every one of these has no `use` statement beyond `constants`, which is
## what makes them safe here; recompute_jacobian.shared.o is present because
## test_gf_shape3D uses it as the oracle for the fork in gf_shape3D.F90.
## binary_c_io.cc.o and calendar.shared.o are the SAC writer's (gf_sac.F90):
## the solver's own byte writer and its leap-year rule, so that
## test_gf_sac can re-read a header from a plain ./configure.
gf3d_KERNEL_SHARED_OBJECTS = \
	$O/gll_library.shared.o \
	$O/lagrange_poly.shared.o \
	$O/hex_nodes.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/binary_c_io.cc.o \
	$O/calendar.shared.o \
	$(EMPTY_MACRO)

## the parts that read the database, and therefore need HDF5
##
## gf3d.gf3d.o is the public Fortran module -- a re-export facade, so that a
## downstream program needs `use gf3d` and nothing else -- and
## gf3d_capi.gf3d.o is the bind(C) facade behind include/gf3d.h and the
## Python package. Both sit here rather than with the kernels because they
## re-export the database side.
gf3d_HDF5_OBJECTS = \
	$O/gf_hdf5_read.gf3d.o \
	$O/gf_database.gf3d.o \
	$O/gf_element_io.gf3d.o \
	$O/gf_locate.gf3d.o \
	$O/gf_seismograms.gf3d.o \
	$O/gf3d.gf3d.o \
	$O/gf3d_capi.gf3d.o \
	$(EMPTY_MACRO)

## From src/specfem3D/: the FORCESOLUTION and CMTSOLUTION readers, reused
## rather than re-ported so that the non-dimensionalisation (scaleF, and
## scaleM from Stage 4) comes from the solver's own source and cannot drift.
## Both are free of MPI and of HDF5.
gf3d_SOLVER_OBJECTS = \
	$O/get_force.solver.o \
	$O/get_cmt.solver.o \
	$(EMPTY_MACRO)

## get_cmt() calls julian_day() to turn the CMTSOLUTION PDE header into a
## Julian day, so calendar.shared.o comes with it. No `use` statements and
## no MPI.

## library contents (everything except the program itself)
gf3d_OBJECTS = \
	$(gf3d_KERNEL_OBJECTS) \
	$(gf3d_HDF5_OBJECTS) \
	$(gf3d_SOLVER_OBJECTS) \
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
##
## make_ellipticity.shared.o arrives the same way model_topo_bathy.shared.o
## does, and it brings the largest tail in this list. All the library calls
## is add_ellipticity_rtheta(), 35 lines with no dependency beyond
## spline_evaluation(); but a static archive pulls in whole objects, so it
## also brings make_ellipticity()/make_ellipticity_r(), which build the
## spline table from a planet's density profile and therefore reference
## prem_density(), sohl_density(), model_vpremoon_density(),
## get_model_sohl_radii() and intgrl() -- hence the four model objects and
## intgrl.shared.o below. None of it is reachable from library code: the
## database ships its own spline table in mesh_info.h5, so gf3d only ever
## *evaluates* the ellipticity, never builds it. The two exit_MPI calls that
## come with it (make_ellipticity.f90:223,331) are satisfied by
## $O/gf_mpi_stubs.gf3d.o, exactly as model_topo_bathy's are.
##
## The alternative would be to inline those 35 lines into gf_geometry.F90 and
## drop five objects. That is rejected on purpose: applying ellipticity to
## the surface radius before subtracting depth, rather than after, is
## precisely the stale convention GF3DF carries and this port exists to avoid
## (see the header of gf_geometry.F90). A convention that subtle should have
## exactly one definition in the tree, and it should be the solver's.
gf3d_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/flush_system.shared.o \
	$O/model_topo_bathy.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/reduce.shared.o \
	$O/make_ellipticity.shared.o \
	$O/spline_routines.shared.o \
	$O/intgrl.shared.o \
	$O/model_prem.shared.o \
	$O/model_Sohl.shared.o \
	$O/model_vpremoon.shared.o \
	$O/heap_sort.shared.o \
	$O/search_kdtree.shared.o \
	$(gf3d_KERNEL_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

#######################################
##
## -fPIC objects, for lib/libgf3d.so
##
## A shared object needs every one of its objects compiled -fPIC, and
## $O/*.shared.o and $O/*.solver.o are not: they belong to the solver's
## build, which must not change (this is PR 1's invariant, and a global
## -fPIC would break it).
##
## Two halves, and they are handled differently:
##
##  - src/gf3d/'s own objects are compiled -fPIC *always*, in place. Nothing
##    but gf3d links them, position-independent code costs nothing
##    measurable on x86-64, and a second copy of every gf_*.o would mean a
##    second copy of every gf_*.mod and a module-file race under `make -j`.
##
##  - the src/shared and src/specfem3D sources the library reuses get a
##    second object suffix, $O/%.gfpic.o, and their module files go into
##    $O/gfpic/ so that they can never race the solver's own compile of the
##    same source.
##
## The twins are only ever built for the .so; lib/libgf3d.a keeps using the
## ordinary $O/*.shared.o, so a Fortran or C caller of the static library
## links exactly what it did before.
gf3d_PICMODDIR = $O/gfpic

## Module output and search flags for the twins. Deliberately *without*
## $(FC_MODINC)$O: if a twin is missing a dependency below, it must fail to
## compile rather than quietly pick up the solver's .mod file from $O.
FCFLAGS_gfpic = -I${SETUP} $(FC_MODOUT)$(gf3d_PICMODDIR) $(FC_MODINC)$(gf3d_PICMODDIR) -I$B

## Overridable, for a compiler that spells it differently
## (e.g. `make FC_PICFLAG=-PIC gf3d` for some non-GNU/Intel front ends).
FC_PICFLAG ?= -fPIC
CC_PICFLAG ?= -fPIC

gf3d_PIC_SHARED_OBJECTS = \
	$O/shared_par.gfpic.o \
	$O/flush_system.gfpic.o \
	$O/model_topo_bathy.gfpic.o \
	$O/rthetaphi_xyz.gfpic.o \
	$O/reduce.gfpic.o \
	$O/make_ellipticity.gfpic.o \
	$O/spline_routines.gfpic.o \
	$O/intgrl.gfpic.o \
	$O/model_prem.gfpic.o \
	$O/model_Sohl.gfpic.o \
	$O/model_vpremoon.gfpic.o \
	$O/heap_sort.gfpic.o \
	$O/search_kdtree.gfpic.o \
	$O/gll_library.gfpic.o \
	$O/lagrange_poly.gfpic.o \
	$O/hex_nodes.gfpic.o \
	$O/recompute_jacobian.gfpic.o \
	$O/calendar.gfpic.o \
	$O/binary_c_io.gfpic_cc.o \
	$O/get_force.gfpic.o \
	$O/get_cmt.gfpic.o \
	$(EMPTY_MACRO)

## everything that goes into the shared object: the gf3d objects (already
## -fPIC) plus the twins, with gf3d_SOLVER_OBJECTS replaced by theirs
gf3d_PIC_OBJECTS = \
	$(gf3d_KERNEL_OBJECTS) \
	$(gf3d_HDF5_OBJECTS) \
	$(gf3d_PIC_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

## the public C header and Fortran module, installed for downstream callers
gf3d_INCDIR = $B/include

gf3d_INCLUDES = \
	$(gf3d_INCDIR)/gf3d.h \
	$(gf3d_INCDIR)/gf3d.$(FC_MODEXT) \
	$(EMPTY_MACRO)

## Appended to gf3d_MODULES so that the central `clean` in Makefile.in
## removes them. It has to be unconditional: gf3d_TARGETS is empty when the
## tree was configured without HDF5, and `realclean` clears $E, $O and $L but
## never include/.
gf3d_CLEAN_EXTRA = \
	$L/libgf3d.so \
	$(gf3d_PIC_SHARED_OBJECTS) \
	$(gf3d_PICMODDIR)/*.$(FC_MODEXT) \
	$(gf3d_INCLUDES) \
	$(EMPTY_MACRO)

#######################################

gf3d_MODULES = \
	$(FC_MODDIR)/gf_par.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_hdf5_read.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_shared_params.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_database.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_shape3d.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_geometry.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_geo_chain.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_element_io.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_locate.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_interp.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_strain.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_moment.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_stf.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_source.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_partials.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_sac.$(FC_MODEXT) \
	$(FC_MODDIR)/gf_seismograms.$(FC_MODEXT) \
	$(FC_MODDIR)/gf3d.$(FC_MODEXT) \
	$(FC_MODDIR)/gf3d_capi.$(FC_MODEXT) \
	$(gf3d_CLEAN_EXTRA) \
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
	$L/libgf3d.so \
	$E/xgf3d \
	$(gf3d_INCLUDES) \
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
gf3d_kernels: $(gf3d_KERNEL_OBJECTS) $(gf3d_KERNEL_SHARED_OBJECTS)

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

####
#### the shared object, for ctypes and for any C caller that dlopen()s us
####

## The HDF5 library directories are baked in as an rpath, so that
## `ctypes.CDLL("libgf3d.so")` resolves them with no LD_LIBRARY_PATH set --
## a Python user has no reason to know where this tree's HDF5 came from.
## $(LDFLAGS) is where configure put them (Makefile.in, COND_HDF5).
comma := ,
gf3d_RPATH = $(patsubst -L%,-Wl$(comma)-rpath$(comma)%,$(filter -L%,$(LDFLAGS)))

## HDF5_USE_SHLIB: when FC is HDF5's own h5fc/h5pfc wrapper, it appends the
## *static* HDF5 archives by default, and those are not compiled -fPIC, so
## the shared link dies on the first archive member the linker actually
## needs (libhdf5_hl, which Makefile.in leaves to the wrapper). Asking the
## wrapper for the shared libraries instead fixes it at the source. The
## variable is meaningless to a plain gfortran/ifort, where HDF5_LIBS names
## the library directory and the .so files carry their own dependencies.
gf3d_SO_ENV = HDF5_USE_SHLIB=yes

## Note the .so is deliberately *not* in DEFAULT: `make` builds xgf3d, and
## `make gf3d` builds this as well. A shared link is the one step here that
## can fail on an unusual toolchain, and it must not take the executable
## down with it.
##
## Note also that with a parallel HDF5 the .so lists libmpi as NEEDED, pulled
## in transitively by libhdf5 itself, exactly as bin/xgf3d already does. That
## is HDF5's dependency and not ours; what is ours is checked on the objects,
## `nm $O/gf_*.o | grep ' U .*mpi_'`, by tests/gf3d/.
$L/libgf3d.so: $(gf3d_PIC_OBJECTS)
	@-mkdir -p $L
	$(gf3d_SO_ENV) ${FCCOMPILE_CHECK} $(FC_PICFLAG) -shared -o $@ $(gf3d_PIC_OBJECTS) $(LDFLAGS) $(gf3d_RPATH)

####
#### the public header and module
####

## include/ is created by configure and is gitignored, so the tracked copy of
## the header lives beside the source it describes and is installed here.
$(gf3d_INCDIR)/gf3d.h: ${S_TOP}/src/gf3d/gf3d.h
	@-mkdir -p $(gf3d_INCDIR)
	cp -f $< $@

$(gf3d_INCDIR)/gf3d.$(FC_MODEXT): $O/gf3d.gf3d.o
	@-mkdir -p $(gf3d_INCDIR)
	cp -f $(FC_MODDIR)/gf3d.$(FC_MODEXT) $@

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
$O/gf_shape3D.gf3d.o: $O/gf_par.gf3d.o
$O/gf_geometry.gf3d.o: $O/gf_par.gf3d.o $O/gf_shape3D.gf3d.o
$O/gf_geo_chain.gf3d.o: $O/gf_par.gf3d.o $O/gf_geometry.gf3d.o
$O/gf_element_io.gf3d.o: $O/gf_par.gf3d.o $O/gf_hdf5_read.gf3d.o
$O/gf_locate.gf3d.o: $O/gf_par.gf3d.o $O/gf_database.gf3d.o $O/gf_element_io.gf3d.o \
                     $O/gf_geometry.gf3d.o $O/gf_shape3D.gf3d.o $O/search_kdtree.shared.o
$O/gf_interp.gf3d.o: $O/gf_par.gf3d.o
$O/gf_strain.gf3d.o: $O/gf_par.gf3d.o
$O/gf_moment.gf3d.o: $O/gf_par.gf3d.o $O/gf_strain.gf3d.o
$O/gf_stf.gf3d.o: $O/gf_par.gf3d.o
$O/gf_source.gf3d.o: $O/gf_par.gf3d.o
$O/gf_partials.gf3d.o: $O/gf_par.gf3d.o $O/gf_strain.gf3d.o $O/gf_moment.gf3d.o $O/gf_stf.gf3d.o
$O/gf_sac.gf3d.o: $O/gf_par.gf3d.o $O/gf_partials.gf3d.o
$O/gf_seismograms.gf3d.o: $O/gf_par.gf3d.o $O/gf_database.gf3d.o $O/gf_element_io.gf3d.o \
                          $O/gf_interp.gf3d.o $O/gf_source.gf3d.o $O/gf_strain.gf3d.o \
                          $O/gf_moment.gf3d.o $O/gf_stf.gf3d.o $O/gf_partials.gf3d.o \
                          $O/gf_geo_chain.gf3d.o
$O/gf3d.gf3d.o: $O/gf_par.gf3d.o $O/gf_database.gf3d.o $O/gf_locate.gf3d.o \
                $O/gf_source.gf3d.o $O/gf_seismograms.gf3d.o $O/gf_partials.gf3d.o \
                $O/gf_sac.gf3d.o $O/gf_shared_params.gf3d.o
$O/gf3d_capi.gf3d.o: $O/gf_par.gf3d.o $O/gf_database.gf3d.o $O/gf_locate.gf3d.o \
                     $O/gf_source.gf3d.o $O/gf_seismograms.gf3d.o $O/gf_partials.gf3d.o \
                     $O/gf_shared_params.gf3d.o
$O/gf3d_main.gf3d.o: $O/gf_par.gf3d.o $O/gf_database.gf3d.o $O/gf_locate.gf3d.o \
                     $O/gf_source.gf3d.o $O/gf_seismograms.gf3d.o $O/gf_partials.gf3d.o \
                     $O/gf_sac.gf3d.o

## unique object suffix: every rules.mk writes into the same $O, so the
## pattern rules of different subdirectories must not collide
##
## $(FC_PICFLAG): these objects go into lib/libgf3d.so as well as into
## lib/libgf3d.a, and nothing outside src/gf3d/ links them -- see the
## -fPIC section above.
$O/%.gf3d.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(FC_PICFLAG) -c -o $@ $<

$O/%.gf3d.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(FC_PICFLAG) -c -o $@ $<

$O/%.gf3d_cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(CC_PICFLAG) -o $@ $<

####
#### the -fPIC twins of the src/shared and src/specfem3D sources
####

$(gf3d_PICMODDIR):
	@-mkdir -p $@

## shared_par.f90 defines `constants` and `shared_parameters`, which every
## other twin uses, so it is named explicitly: the pattern rule below would
## otherwise make it a prerequisite of itself.
$O/shared_par.gfpic.o: ${S_TOP}/src/shared/shared_par.f90 ${SETUP}/constants.h | $(gf3d_PICMODDIR)
	${FCCOMPILE_CHECK} $(FCFLAGS_gfpic) $(FC_PICFLAG) -c -o $@ $<

## the only inter-twin module dependency in the whole set: every other source
## here uses nothing but `constants`/`shared_parameters` or a module defined
## in its own file. (search_kdtree -> heap_sort and get_cmt -> julian_day are
## plain calls, so they matter at link time only.)
## Mirrors src/shared/rules.mk.
$O/make_ellipticity.gfpic.o: $O/model_prem.gfpic.o $O/model_Sohl.gfpic.o $O/model_vpremoon.gfpic.o

$O/%.gfpic.o: ${S_TOP}/src/shared/%.f90 $O/shared_par.gfpic.o | $(gf3d_PICMODDIR)
	${FCCOMPILE_CHECK} $(FCFLAGS_gfpic) $(FC_PICFLAG) -c -o $@ $<

$O/%.gfpic.o: ${S_TOP}/src/specfem3D/%.f90 $O/shared_par.gfpic.o | $(gf3d_PICMODDIR)
	${FCCOMPILE_CHECK} $(FCFLAGS_gfpic) $(FC_PICFLAG) -c -o $@ $<

$O/%.gfpic_cc.o: ${S_TOP}/src/shared/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(CC_PICFLAG) -o $@ $<
