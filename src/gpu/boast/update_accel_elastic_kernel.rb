require "./update_disp_veloc_kernel.rb"

module BOAST

  def BOAST::update_accel_elastic_kernel(ref = true)
    BOAST::update_kernel(:accel_elastic, ref)
  end

end
