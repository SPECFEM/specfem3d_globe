require "./compute_stacey_elastic_kernel"

module BOAST

  def BOAST::compute_stacey_elastic_backward_kernel(ref = true, n_dim = 3, n_gllx = 5, n_gll2 = 25)
    BOAST::compute_stacey_elastic_k(:backward, ref, n_dim, n_gllx, n_gll2)
  end

end
