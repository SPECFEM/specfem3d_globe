require "./update_disp_veloc_kernel.rb"

module BOAST

  def BOAST::update_accel_acoustic_kernel(ref = true)
    BOAST::update_kernel(:accel_acoustic, ref)
  end

end
