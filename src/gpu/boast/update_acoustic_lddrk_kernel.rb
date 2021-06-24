require "./update_disp_veloc_LDDRK_kernel.rb"

module BOAST

  def BOAST::update_acoustic_lddrk_kernel(ref = true)
    BOAST::update_lddrk_kernel(:acoustic, ref)
  end

end
