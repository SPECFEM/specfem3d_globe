module BOAST

  def BOAST::update_lddrk_kernel(type, ref)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    displ  = Real("displ",        :dir => :inout, :dim => [Dim()] )
    veloc  = Real("veloc",        :dir => :inout, :dim => [Dim()] )
    accel  = Real("accel",        :dir => :in,    :dim => [Dim()] )
    displ_lddrk  = Real("displ_lddrk",        :dir => :inout, :dim => [Dim()] )
    veloc_lddrk  = Real("veloc_lddrk",        :dir => :inout, :dim => [Dim()] )
    alpha  = Real("alpha_lddrk", :dir => :in )
    beta   = Real("beta_lddrk",  :dir => :in )
    deltat = Real("deltat",      :dir => :in )
    size = Int("size", :dir => :in)

    case type
    when :acoustic, :elastic
      function_name = "update_#{type}_lddrk_kernel"
      variables = [ displ, veloc, accel, displ_lddrk, veloc_lddrk, alpha, beta, deltat, size ]
    else
      raise "Invalid update type #{type}!"
    end

    p = Procedure(function_name, variables)
    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      make_specfem3d_header
      open p
      decl id = Int("id")
      comment()

      print id === get_global_id(0) + get_group_id(1)*get_global_size(0)
      comment()

      print If(id < size ) {
        case type
        when :acoustic
          print veloc_lddrk[id] === alpha * veloc_lddrk[id] + deltat * accel[id]
          print displ_lddrk[id] === alpha * displ_lddrk[id] + deltat * veloc[id]

          print veloc[id] === veloc[id] + beta * veloc_lddrk[id]
          print displ[id] === displ[id] + beta * displ_lddrk[id]

        when :elastic
          print veloc_lddrk[id*3    ] === alpha * veloc_lddrk[id*3    ] + deltat * accel[id*3    ]
          print veloc_lddrk[id*3 + 1] === alpha * veloc_lddrk[id*3 + 1] + deltat * accel[id*3 + 1]
          print veloc_lddrk[id*3 + 2] === alpha * veloc_lddrk[id*3 + 2] + deltat * accel[id*3 + 2]

          print displ_lddrk[id*3    ] === alpha * displ_lddrk[id*3    ] + deltat * veloc[id*3    ]
          print displ_lddrk[id*3 + 1] === alpha * displ_lddrk[id*3 + 1] + deltat * veloc[id*3 + 1]
          print displ_lddrk[id*3 + 2] === alpha * displ_lddrk[id*3 + 2] + deltat * veloc[id*3 + 2]

          print veloc[id*3    ] === veloc[id*3    ] + beta * veloc_lddrk[id*3    ]
          print veloc[id*3 + 1] === veloc[id*3 + 1] + beta * veloc_lddrk[id*3 + 1]
          print veloc[id*3 + 2] === veloc[id*3 + 2] + beta * veloc_lddrk[id*3 + 2]

          print displ[id*3    ] === displ[id*3    ] + beta * displ_lddrk[id*3    ]
          print displ[id*3 + 1] === displ[id*3 + 1] + beta * displ_lddrk[id*3 + 1]
          print displ[id*3 + 2] === displ[id*3 + 2] + beta * displ_lddrk[id*3 + 2]
        end
      }
      close p
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end

end
