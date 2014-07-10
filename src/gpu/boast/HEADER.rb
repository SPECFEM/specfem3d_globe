module BOAST
  def BOAST::protect(name,val)
    BOAST::get_output.puts "#ifndef #{name}"
    BOAST::get_output.puts "#define #{name} #{val}"
    BOAST::get_output.puts "#endif"
  end
  def BOAST::make_specfem3d_header(opts = {})
    if BOAST::get_lang == CL then
      if BOAST::get_default_real_size == 8 or opts[:double] then
        BOAST::get_output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
        BOAST::get_output.puts "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable"
      end
      load "./atomicAdd_f.rb"
    end
    load "./INDEX2.rb"
    load "./INDEX3.rb"
    load "./INDEX4.rb"
    load "./INDEX5.rb"
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
  end
end
