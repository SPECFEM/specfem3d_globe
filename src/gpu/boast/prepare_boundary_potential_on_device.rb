require './prepare_boundary_accel_on_device.rb'

module BOAST
  def BOAST::prepare_boundary_potential_on_device(ref = true)
    BOAST::prepare_boundary_on_device(:potential, ref)
  end
end
