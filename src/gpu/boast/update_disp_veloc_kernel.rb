module BOAST

  def BOAST::update_disp_veloc_kernel(ref = true)
    BOAST::update_kernel(:disp_veloc, ref)
  end

  def BOAST::update_kernel(type, ref)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    size = Int("size", :dir => :in)
    deltat         = Real("deltat",       :dir => :in )
    deltatsqover2  = Real("deltatsqover2",:dir => :in )
    deltatover2    = Real("deltatover2",  :dir => :in )
    case type
    when :disp_veloc
      function_name = "update_disp_veloc_kernel"
      quantity      = Real("displ",        :dir => :inout, :dim => [Dim()] )
      first_deriv   = Real("veloc",        :dir => :inout, :dim => [Dim()] )
      second_deriv  = Real("accel",        :dir => :out,   :dim => [Dim()] )
      variables = [ quantity, first_deriv, second_deriv, size, deltat, deltatsqover2, deltatover2 ]
    when :potential
      function_name = "update_potential_kernel"
      quantity      = Real("potential_acoustic",         :dir => :inout, :dim => [Dim()] )
      first_deriv   = Real("potential_dot_acoustic",     :dir => :inout, :dim => [Dim()] )
      second_deriv  = Real("potential_dot_dot_acoustic", :dir => :out,   :dim => [Dim()] )
      variables = [ quantity, first_deriv, second_deriv, size, deltat, deltatsqover2, deltatover2 ]
    when :accel_elastic
      function_name = "update_accel_elastic_kernel"
      accel           = Real("accel",           :dir => :inout, :dim => [Dim()] )
      veloc           = Real("veloc",           :dir => :in,    :dim => [Dim()] )
      two_omega_earth = Real("two_omega_earth", :dir => :in)
      rmassx          = Real("rmassx",          :dir => :in,    :dim => [Dim()] )
      rmassy          = Real("rmassy",          :dir => :in,    :dim => [Dim()] )
      rmassz          = Real("rmassz",          :dir => :in,    :dim => [Dim()] )
      variables = [ accel, veloc, size, two_omega_earth, rmassx, rmassy, rmassz ]
    when :veloc_elastic, :veloc_acoustic
      function_name = "update_#{type}_kernel"
      veloc = Real("veloc", :dir => :inout, :dim => [Dim()] )
      accel = Real("accel", :dir => :in,    :dim => [Dim()] )
      variables = [ veloc, accel, size, deltatover2 ]
    when :accel_acoustic
      function_name = "update_accel_acoustic_kernel"
      accel = Real("accel", :dir => :inout, :dim => [Dim()] )
      rmass = Real("rmass", :dir => :in,    :dim => [Dim()] )
      variables = [ accel, size, rmass ]
    else
      raise "Invalid update type #{type}!"
    end

    p = Procedure(function_name, variables)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header
      decl p
      decl id = Int("id")
      print id === get_global_id(0) + get_group_id(1)*get_global_size(0)
      print If( id < size ) {
        case type
        when :disp_veloc, :potential
          print quantity[id] === quantity[id] + deltat*first_deriv[id] + deltatsqover2*second_deriv[id]
          print first_deriv[id] === first_deriv[id] + deltatover2*second_deriv[id]
          print second_deriv[id] === 0.0
        when :accel_elastic
          print accel[id*3    ] === accel[id*3    ]*rmassx[id] + two_omega_earth*veloc[id*3 + 1]
          print accel[id*3 + 1] === accel[id*3 + 1]*rmassy[id] - two_omega_earth*veloc[id*3    ]
          print accel[id*3 + 2] === accel[id*3 + 2]*rmassz[id]
        when :veloc_elastic
          print veloc[id*3    ] === veloc[id*3    ] + deltatover2*accel[id*3    ]
          print veloc[id*3 + 1] === veloc[id*3 + 1] + deltatover2*accel[id*3 + 1]
          print veloc[id*3 + 2] === veloc[id*3 + 2] + deltatover2*accel[id*3 + 2]
        when :accel_acoustic
          print accel[id] === accel[id]*rmass[id]
        when :veloc_acoustic
          print veloc[id] === veloc[id] + deltatover2*accel[id]
        end
      }

      close p
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end

end
