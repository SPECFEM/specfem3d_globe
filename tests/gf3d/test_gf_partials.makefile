# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_partials

## compilation directories
O := ./obj

# Kernel objects only -- no HDF5 archive exists after
# 0.configure.default_make.sh, which is the point of keeping
# gf3d_KERNEL_OBJECTS separate in src/gf3d/rules.mk. Stage 8's sections
# drive the geometry, the geographic chain and the strain kernels directly,
# so those kernels and the shared routines behind them are linked as well.
OBJECTS = \
	$O/gf_par.gf3d.o \
	$O/gf_shape3D.gf3d.o \
	$O/gf_geometry.gf3d.o \
	$O/gf_geo_chain.gf3d.o \
	$O/gf_interp.gf3d.o \
	$O/gf_strain.gf3d.o \
	$O/gf_moment.gf3d.o \
	$O/gf_stf.gf3d.o \
	$O/gf_partials.gf3d.o \
	$O/gf_mpi_stubs.gf3d.o \
	$O/shared_par.shared_module.o \
	$O/gll_library.shared.o \
	$O/lagrange_poly.shared.o \
	$O/hex_nodes.shared.o \
	$(EMPTY_MACRO)

# The geographic chain's shared routines. gf_geometry calls the solver's
# lat_2_geocentric_colat_dble, reduce and add_ellipticity_rtheta; the last
# arrives with make_ellipticity's whole object, whose tail references the
# three planet density models and intgrl (see src/gf3d/rules.mk), and whose
# two exit_MPI calls are satisfied by gf_mpi_stubs. None of these is in
# gf3d_KERNEL_SHARED_OBJECTS, so they are built here through the pattern
# rules `include Makefile` brings in from src/shared/rules.mk -- the models
# first, because make_ellipticity uses their parameter modules.
GEO = \
	$O/model_prem.shared.o \
	$O/model_Sohl.shared.o \
	$O/model_vpremoon.shared.o \
	$O/intgrl.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/reduce.shared.o \
	$O/spline_routines.shared.o \
	$O/make_ellipticity.shared.o \
	$(EMPTY_MACRO)

# Serial link: ${FCCOMPILE_CHECK}, no $(MPILIBS), no $(LDFLAGS).
test_gf_partials: $(GEO)
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_partials \
		gf_manufactured.f90 test_gf_partials.f90 -I./obj $(OBJECTS) $(GEO)
