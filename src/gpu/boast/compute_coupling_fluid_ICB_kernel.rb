require './compute_coupling_fluid_CMB_kernel.rb'
module BOAST
  def BOAST::compute_coupling_fluid_ICB_kernel(ref=true, n_dim = 3, n_gllx = 5)
    return BOAST::compute_coupling_kernel(ref, :fluid_ICB, n_dim, n_gllx)
  end
end
