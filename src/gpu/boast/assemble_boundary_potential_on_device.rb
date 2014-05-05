require './assemble_boundary_accel_on_device.rb'
module BOAST
  def BOAST::assemble_boundary_potential_on_device(ref = true)
    BOAST::assemble_boundary_on_device(:potential, ref)
  end
end

