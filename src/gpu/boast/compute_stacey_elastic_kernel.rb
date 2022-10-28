module BOAST

  def BOAST::compute_stacey_elastic_kernel(ref = true, n_dim = 3, n_gllx = 5, n_gll2 = 25)
    BOAST::compute_stacey_elastic_k(:forward, ref, n_dim, n_gllx, n_gll2)
  end

  def BOAST::compute_stacey_elastic_k(type, ref = true, n_dim = 3, n_gllx = 5, n_gll2 = 25)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    if type == :forward then
      function_name = "compute_stacey_elastic_kernel"
      veloc                    = Real("veloc",                   :dir => :in,   :dim => [ Dim() ])
      accel                    = Real("accel",                   :dir => :inout,:dim => [ Dim() ])
      abs_boundary_normal      = Real("abs_boundary_normal",     :dir => :in,   :dim => [ Dim() ])
      abs_boundary_jacobian2Dw = Real("abs_boundary_jacobian2Dw", :dir => :in,   :dim => [ Dim() ])
      rho_vp                   = Real("rho_vp",                  :dir => :in,   :dim => [ Dim() ])
      rho_vs                   = Real("rho_vs",                  :dir => :in,   :dim => [ Dim() ])
      save_stacey              = Int( "SAVE_STACEY",             :dir => :in)
      b_absorb_field           = Real("b_absorb_field",          :dir => :out,  :dim => [ Dim() ])
      variables = [veloc, accel]
    elsif type == :backward then
      function_name = "compute_stacey_elastic_backward_kernel"
      b_accel                 = Real("b_accel",                 :dir => :inout,:dim => [ Dim() ])
      b_absorb_field          = Real("b_absorb_field",          :dir => :in,   :dim => [ Dim() ])
      variables = [b_accel, b_absorb_field]
    else
      raise "Unsupported type : #{type}!"
    end

    num_abs_boundary_faces = Int( "num_abs_boundary_faces", :dir => :in)
    abs_boundary_ispec     = Int( "abs_boundary_ispec",     :dir => :in,   :dim => [ Dim() ])
    abs_boundary_npoin     = Int( "abs_boundary_npoin",     :dir => :in,   :dim => [ Dim() ])
    abs_boundary_ijk       = Int( "abs_boundary_ijk",       :dir => :in,   :dim => [ Dim() ])
    ibool                  = Int( "ibool",                  :dir => :in,   :dim => [ Dim() ])
    variables += [ num_abs_boundary_faces, abs_boundary_ispec, abs_boundary_npoin, abs_boundary_ijk ]
    if type == :forward then
      variables += [ abs_boundary_normal, abs_boundary_jacobian2Dw, ibool, rho_vp, rho_vs, save_stacey, b_absorb_field ]
    elsif type == :backward then
      variables += [ ibool ]
     end

    ndim  = Int("NDIM",  :const => n_dim)
    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)

    p = Procedure(function_name, variables)
    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      make_specfem3d_header( :ndim => n_dim, :ngllx => n_gllx, :ngll2 => n_gll2 )
      open p
      decl npoin = Int("npoin")
      decl igll = Int("igll")
      decl iface = Int("iface")
      decl i = Int("i"), j = Int("j"), k = Int("k")
      decl iglob = Int("iglob")
      decl ispec = Int("ispec")
      if type == :forward then
        decl *v = [Real("vx"), Real("vy"), Real("vz") ]
        decl vn = Real("vn")
        decl *n = [Real("nx"), Real("ny"), Real("nz") ]
        decl rho_vp_temp = Real("rho_vp_temp")
        decl rho_vs_temp = Real("rho_vs_temp")
        decl *t = [Real("tx"), Real("ty"), Real("tz") ]
        decl weight = Real("weight")
      end
      comment()

      print igll === get_local_id(0)
      print iface === get_group_id(0)+get_group_id(1)*get_num_groups(0)
      comment()

      print If(iface < num_abs_boundary_faces ) {

        print npoin === abs_boundary_npoin[iface]
        comment()

        print If(igll < npoin ) {
          print ispec === abs_boundary_ispec[iface]-1
          comment()

          print i === abs_boundary_ijk[INDEX3(3,ngll2,0,igll,iface)] - 1
          print j === abs_boundary_ijk[INDEX3(3,ngll2,1,igll,iface)] - 1
          print k === abs_boundary_ijk[INDEX3(3,ngll2,2,igll,iface)] - 1
          comment()

          print iglob === ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] - 1
          comment()

          if type == :forward then
            (0..2).each { |indx|
              print v[indx] === veloc[iglob*3+indx]
            }
            (0..2).each { |indx|
              print n[indx] === abs_boundary_normal[INDEX3(ndim,ngll2,indx,igll,iface)]
            }
            comment()
            print vn === v[0]*n[0] + v[1]*n[1] + v[2]*n[2]
            print rho_vp_temp === rho_vp[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)]
            print rho_vs_temp === rho_vs[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)]
            comment()
            (0..2).each { |indx|
              print t[indx] === rho_vp_temp*vn*n[indx] + rho_vs_temp*(v[indx] - vn*n[indx])
            }
            comment()
            print weight === abs_boundary_jacobian2Dw[INDEX2(ngll2,igll,iface)]
            comment()
            (0..2).each { |indx|
              print atomicAdd(accel + iglob*3 + indx, -t[indx]*weight)
            }
            comment()
            print If(save_stacey) {
              (0..2).each { |indx|
                print b_absorb_field[INDEX3(ndim,ngll2,indx,igll,iface)] === t[indx]*weight
              }
            }
          else
            (0..2).each { |indx|
              print atomicAdd(b_accel + iglob*3 + indx, -b_absorb_field[INDEX3(ndim,ngll2,indx,igll,iface)])
            }
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

