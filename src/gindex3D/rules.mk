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

gindex3D_TARGETS = \
	$E/xgindex3D \
	$(EMPTY_MACRO)

gindex3D_OBJECTS = \
	$O/create_gindex.gindex.o \
	$O/gindex3D.gindex.o \
	$O/gindex3D_par.gindex_module.o \
	$O/initialize_gindex.gindex.o \
	$(EMPTY_MACRO)

gindex3D_MODULES = \
	$(FC_MODDIR)/gindex_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

# These files come from the specfem3D directory
gindex3D_SHARED_OBJECTS = \
	$O/asdf_data.solverstatic_module.o \
	$O/bcast_mesh_databases.solverstatic.o \
	$O/locate_regular_points.solverstatic.o \
	$O/read_arrays_solver.solverstatic.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/read_mesh_databases.solverstatic.o \
	$O/specfem3D_par.solverstatic_module.o \
	$(EMPTY_MACRO)

# These files come from the shared directory
gindex3D_SHARED_OBJECTS += \
	$O/adios_manager.shared_adios_module.o \
	$O/auto_ner.shared.o \
	$O/broadcast_computed_parameters.shared.o \
	$O/count_elements.shared.o \
	$O/count_points.shared.o \
	$O/create_name_database.shared.o \
	$O/define_all_layers.shared.o \
	$O/euler_angles.shared.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/hex_nodes.shared.o \
	$O/memory_eval.shared.o \
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/reduce.shared.o \
	$O/rotate_tensor.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/save_header_file.shared.o \
	$O/shared_par.shared_module.o \
	$O/ylm.shared.o \
	$(EMPTY_MACRO)

###
### GPU
###

gindex3D_SHARED_OBJECTS += $(gpu_OBJECTS)

###
### ADIOS
###

# conditional adios linking
ifeq ($(ADIOS),yes)
	gindex3D_SHARED_OBJECTS += $(adios_specfem3D_OBJECTS)
	gindex3D_SHARED_OBJECTS += $(adios_specfem3D_SHARED_OBJECTS)
else ifeq ($(ADIOS2),yes)
	gindex3D_SHARED_OBJECTS += $(adios_specfem3D_OBJECTS)
	gindex3D_SHARED_OBJECTS += $(adios_specfem3D_SHARED_OBJECTS)
else
	gindex3D_SHARED_OBJECTS += $(adios_specfem3D_SHARED_STUBS)
endif

###
### ASDF
###

# conditional asdf linking
ifeq ($(ASDF),yes)
	gindex3D_SHARED_OBJECTS += $(asdf_specfem3D_OBJECTS)
	gindex3D_SHARED_OBJECTS += $(asdf_specfem3D_SHARED_OBJECTS)
else
	gindex3D_SHARED_OBJECTS += ${asdf_specfem3D_SHARED_STUBS}
endif

# conditional CEM or EMC model
ifeq ($(CEM),yes)
	gindex3D_SHARED_OBJECTS += $O/read_write_netcdf.checknetcdf.o
else ifeq ($(EMC),yes)
	gindex3D_SHARED_OBJECTS += $O/read_write_netcdf.checknetcdf.o
else ifeq ($(NETCDF),yes)
	gindex3D_SHARED_OBJECTS += $O/read_write_netcdf.checknetcdf.o
endif

###
### VTK
###

# conditional vtk linking
ifeq ($(VTK),yes)
	gindex3D_SHARED_OBJECTS += $(vtk_specfem3D_OBJECTS)
else
	gindex3D_SHARED_OBJECTS += $(vtk_specfem3D_STUBS)
endif

#######################################

####
#### rules for executables
####

${E}/xgindex3D: $(gindex3D_SHARED_OBJECTS) $(gindex3D_OBJECTS)
	${FCLINK} -o $@ $+ $(SPECFEM_LINK_FLAGS)

#######################################

## compilation directories
S := ${S_TOP}/src/gindex3D
$(gindex3D_OBJECTS): S = ${S_TOP}/src/gindex3D

####
#### rule for each .o file below
####

## additional module dependencies
$O/gindex3D_par.gindex_module.o: $O/specfem3D_par.solverstatic_module.o


## general rules

$O/%.gindex_module.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gindex.o: $S/%.f90 $O/shared_par.shared_module.o $O/gindex3D_par.gindex_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.gindex.o: $S/%.F90 $O/shared_par.shared_module.o $O/gindex3D_par.gindex_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
