module BOAST

  def BOAST::assemble_boundary_accel_on_device(ref = true)
    BOAST::assemble_boundary_on_device(:accel, ref)
  end

  def BOAST::assemble_boundary_on_device(type, ref)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    num_interfaces        = Int("num_interfaces",        :dir => :in)
    max_nibool_interfaces = Int("max_nibool_interfaces", :dir => :in)
    d_nibool_interfaces   = Int("d_nibool_interfaces",   :dir => :in, :dim => [ Dim(num_interfaces) ])
    d_ibool_interfaces    = Int("d_ibool_interfaces",    :dir => :in, :dim => [ Dim(num_interfaces*max_nibool_interfaces) ])

    if type == :accel then
      function_name = "assemble_boundary_accel_on_device"
      d_accel               = Real("d_accel",              :dir => :inout,:dim => [ Dim() ])
      d_send_accel_buffer   = Real("d_send_accel_buffer",  :dir => :in,   :dim => [ Dim(num_interfaces*max_nibool_interfaces*3) ])
      variables = [d_accel,d_send_accel_buffer]
    elsif type == :potential then
      function_name = "assemble_boundary_potential_on_device"
      d_potential_dot_dot_acoustic    = Real("d_potential_dot_dot_acoustic",    :dir => :inout,:dim => [ Dim() ])
      d_send_potential_dot_dot_buffer = Real("d_send_potential_dot_dot_buffer", :dir => :in,   :dim => [ Dim(num_interfaces*max_nibool_interfaces) ])
      variables = [d_potential_dot_dot_acoustic, d_send_potential_dot_dot_buffer]
    else
      raise "Unsupported type : #{type}!"
    end
    variables += [num_interfaces,max_nibool_interfaces,d_nibool_interfaces,d_ibool_interfaces]

    p = Procedure(function_name, variables)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header
      decl p
      id         = Int("id")
      iglob      = Int("iglob")
      iloc       = Int("iloc")
      iinterface = Int("iinterface")
      decl id
      decl iglob
      decl iloc
      decl iinterface
      print id === get_global_id(0)+get_global_size(0)*get_global_id(1)
      print For(iinterface, 0, num_interfaces-1) {
        print If(id<d_nibool_interfaces[iinterface]) {
          print iloc === id + max_nibool_interfaces*iinterface
          print iglob === d_ibool_interfaces[iloc] - 1
          if type == :accel then
            (0..2).each { |i|
               print atomicAdd(d_accel + iglob*3  + i, d_send_accel_buffer[iloc*3+i])
            }
          else
            print atomicAdd(d_potential_dot_dot_acoustic + iglob, d_send_potential_dot_dot_buffer[iloc])
          end
        }
      }
      close p
    else
      raise "Unsupported language!"
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end
end

