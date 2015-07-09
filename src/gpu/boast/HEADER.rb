module BOAST

  # gets version number
  def BOAST::get_boast_version
    # note: function latest_spec_for() only works if internet connection available
    #       spec = Gem.latest_spec_for('BOAST')
    # looks for actually installed version
    spec = Gem.loaded_specs['BOAST']
    if spec.nil? then
      v = ""
    else
      v = spec.version.to_s
    end
    return v
  end
  
  # automatic generation notice
  def BOAST::automatic_notice()
    v = get_boast_version()
    var = "//note: please do not modify this file manually!\n"
    var += "//      this file has been generated automatically by BOAST version #{v}\n"
    var += "//      by: make boast_kernels\n"
    var += "\n"
    return var
  end

  # specfem header
  def BOAST::specfem_header_text()
    var = "\
/*\n\
!=====================================================================\n\
!\n\
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0\n\
!          --------------------------------------------------\n\
!\n\
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp\n\
!                        Princeton University, USA\n\
!                and CNRS / University of Marseille, France\n\
!                 (there are currently many more authors!)\n\
! (c) Princeton University and CNRS / University of Marseille, April 2014\n\
!\n\
! This program is free software; you can redistribute it and/or modify\n\
! it under the terms of the GNU General Public License as published by\n\
! the Free Software Foundation; either version 2 of the License, or\n\
! (at your option) any later version.\n\
!\n\
! This program is distributed in the hope that it will be useful,\n\
! but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
! GNU General Public License for more details.\n\
!\n\
! You should have received a copy of the GNU General Public License along\n\
! with this program; if not, write to the Free Software Foundation, Inc.,\n\
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.\n\
!\n\
!=====================================================================\n\
*/\n"
    if BOAST::get_lang == CUDA then
    var += "\n"
    end
    return var
  end

  # file header
  def BOAST::specfem3d_globe_header_info()
    # text with specfem header
    info = ""
    info += automatic_notice()
    info += specfem_header_text()
    return info
  end
  
  def BOAST::protect(name,val)
    BOAST::get_output.puts "#ifndef #{name}"
    BOAST::get_output.puts "#define #{name} #{val}"
    BOAST::get_output.puts "#endif"
  end

  def BOAST::make_specfem3d_header(opts = {})

    # atomicAdd function for single precision (real) values in OpenCL
    if BOAST::get_lang == CL then
      if BOAST::get_default_real_size == 8 or opts[:double] then
        BOAST::get_output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
      end
      if BOAST::get_default_real_size == 8 then
        BOAST::get_output.puts "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable"
      end
      load "./atomicAdd_f.rb"
    end

    # index definitions
    load "./INDEX2.rb"
    load "./INDEX3.rb"
    load "./INDEX4.rb"
    load "./INDEX5.rb"

    BOAST::get_output.puts "\n"      

    # constants defaults
    ndim = opts[:ndim].nil? ? 3 : opts[:ndim]
    ngllx = opts[:ngllx].nil? ? 5 : opts[:ngllx]
    ngll2 = opts[:ngll2].nil? ? 25 : opts[:ngll2]
    ngll3 = opts[:ngll3].nil? ? 125 : opts[:ngll3]
    ngll3_padded = opts[:ngll3_padded].nil? ? 128 : opts[:ngll3_padded]
    n_sls = opts[:n_sls].nil? ? 3 : opts[:n_sls]
    
    iregion_crust_mantle = opts[:iregion_crust_mantle].nil? ? 1 : opts[:iregion_crust_mantle]
    iregion_inner_core = opts[:iregion_inner_core].nil? ? 3 : opts[:iregion_inner_core]
    iflag_in_fictitious_cube = opts[:iflag_in_fictitious_cube].nil? ? 11 : opts[:iflag_in_fictitious_cube]
    
    r_earth_km = opts[:r_earth_km].nil? ? 6371.0 : opts[:r_earth_km]
    coloring_min_nspec_inner_core = opts[:coloring_min_nspec_inner_core].nil? ? 1000 : opts[:coloring_min_nspec_inner_core]
    coloring_min_nspec_outer_core = opts[:coloring_min_nspec_outer_core].nil? ? 1000 : opts[:coloring_min_nspec_outer_core]

    blocksize_transfer = opts[:blocksize_transfer].nil? ? 256 : opts[:blocksize_transfer]

    # constants declarations
    protect("NDIM", ndim)
    protect("NGLLX", ngllx)
    protect("NGLL2", ngll2)
    protect("NGLL3", ngll3)
    protect("NGLL3_PADDED", ngll3_padded)
    protect("N_SLS", n_sls)
    protect("IREGION_CRUST_MANTLE", iregion_crust_mantle)
    protect("IREGION_INNER_CORE", iregion_inner_core)
    protect("IFLAG_IN_FICTITIOUS_CUBE", iflag_in_fictitious_cube)
    protect("R_EARTH_KM", r_earth_km)
    protect("COLORING_MIN_NSPEC_INNER_CORE", coloring_min_nspec_inner_core)
    protect("COLORING_MIN_NSPEC_OUTER_CORE", coloring_min_nspec_outer_core)
    protect("BLOCKSIZE_TRANSFER", blocksize_transfer)
    
    BOAST::get_output.puts "\n"
  end
end
