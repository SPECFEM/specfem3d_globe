module BOAST
  def BOAST::get_maximum_vector_kernel(ref = true, block_size_transfer = 256)
    BOAST::get_maximum_kernel(:vector, ref, block_size_transfer) 
  end
end
