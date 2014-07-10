require './compute_coupling_fluid_CMB_kernel.rb'
module BOAST
  def BOAST::compute_coupling_CMB_fluid_kernel(ref = true, n_dim = 3, n_gllx = 5)
    BOAST::compute_coupling_kernel(ref, :CMB_fluid, n_dim, n_gllx)
  end
end
