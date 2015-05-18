require "./inner_core_impl_kernel_forward.rb"

module BOAST
  def BOAST::crust_mantle_impl_kernel_forward( ref = true, elem_per_thread = 1, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = false, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, r_earth_km = 6371.0, n_sls = 3)
    return BOAST::impl_kernel(:crust_mantle, true, ref, elem_per_thread, mesh_coloring, textures_fields, textures_constants, unroll_loops, n_gllx, n_gll2, n_gll3, n_gll3_padded, n_sls, r_earth_km)
  end
end
