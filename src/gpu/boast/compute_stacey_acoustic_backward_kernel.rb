require "./compute_stacey_acoustic_kernel.rb"

module BOAST

  def BOAST::compute_stacey_acoustic_backward_kernel(ref = true, n_gllx = 5, n_gll2 = 25)
    BOAST::compute_stacey_acoustic_k(:acoustic_backward, ref, n_gllx, n_gll2)
  end

end
