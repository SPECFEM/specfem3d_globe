require "./update_disp_veloc_kernel.rb"

module BOAST

  def BOAST::update_potential_kernel(ref = true)
    BOAST::update_kernel(:potential, ref)
  end

end
