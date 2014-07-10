require './compute_coupling_fluid_CMB_kernel.rb'
module BOAST
  def BOAST::compute_coupling_ICB_fluid_kernel(ref=true, n_dim = 3, n_gllx = 5)
    return BOAST::compute_coupling_kernel(ref, :ICB_fluid, n_dim, n_gllx)
  end
end
