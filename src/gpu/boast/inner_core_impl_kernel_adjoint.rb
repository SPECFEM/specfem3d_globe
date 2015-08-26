module BOAST
  def BOAST::inner_core_impl_kernel_adjoint(ref = true, elem_per_thread = 1, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, r_earth_km = 6371.0, coloring_min_nspec_inner_core = 1000, i_flag_in_fictitious_cube = 11)
    return BOAST::impl_kernel(:inner_core, false, ref, elem_per_thread, mesh_coloring, textures_fields, textures_constants, unroll_loops, n_gllx, n_gll2, n_gll3, n_gll3_padded, n_sls, r_earth_km, coloring_min_nspec_inner_core, i_flag_in_fictitious_cube)
  end
end

