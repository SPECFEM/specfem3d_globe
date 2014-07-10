module BOAST

  def BOAST::compute_stacey_acoustic_kernel(ref = true, n_gllx = 5, n_gll2 = 25)
    BOAST::compute_stacey_acoustic_k(:acoustic_forward, ref, n_gllx, n_gll2)
  end

  def BOAST::compute_stacey_acoustic_k(type, ref = true, n_gllx = 5, n_gll2 = 25)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    if type == :acoustic_forward then
      function_name = "compute_stacey_acoustic_kernel"
      b_absorb_potential            = Real("b_absorb_potential",         :dir => :out,  :dim => [ Dim() ])
      potential_dot_acoustic        = Real("potential_dot_acoustic",     :dir => :in,   :dim => [ Dim() ])
      potential_dot_dot_acoustic    = Real("potential_dot_dot_acoustic", :dir => :inout,:dim => [ Dim() ])
      abs_boundary_jacobian2D       = Real("abs_boundary_jacobian2D",    :dir => :in,   :dim => [ Dim() ])
      wgllwgll                      = Real("wgllwgll",                   :dir => :in,   :dim => [ Dim() ])
      vpstore                       = Real("vpstore",                    :dir => :in,   :dim => [ Dim() ])
      save_forward                  = Int( "SAVE_FORWARD",               :dir => :in)
      variables = [potential_dot_acoustic, potential_dot_dot_acoustic]
    elsif type == :acoustic_backward then
      function_name = "compute_stacey_acoustic_backward_kernel"
      b_potential_dot_dot_acoustic  = Real("b_potential_dot_dot_acoustic", :dir => :inout,:dim => [ Dim() ])
      b_absorb_potential            = Real("b_absorb_potential",           :dir => :in,  :dim => [ Dim() ])
      variables = [b_potential_dot_dot_acoustic, b_absorb_potential]
    else
      raise "Unsupported type : #{type}!"
    end

    interface_type                = Int( "interface_type",             :dir => :in)
    num_abs_boundary_faces        = Int( "num_abs_boundary_faces",     :dir => :in)
    abs_boundary_ispec            = Int( "abs_boundary_ispec",         :dir => :in,   :dim => [ Dim() ])
    nkmin_xi                      = Int( "nkmin_xi",                   :dir => :in,   :dim => [ Dim() ])
    nkmin_eta                     = Int( "nkmin_eta",                  :dir => :in,   :dim => [ Dim() ])
    njmin                         = Int( "njmin",                      :dir => :in,   :dim => [ Dim() ])
    njmax                         = Int( "njmax",                      :dir => :in,   :dim => [ Dim() ])
    nimin                         = Int( "nimin",                      :dir => :in,   :dim => [ Dim() ])
    nimax                         = Int( "nimax",                      :dir => :in,   :dim => [ Dim() ])
    ibool                         = Int( "ibool",                      :dir => :in,   :dim => [ Dim() ])
    variables += [ interface_type, num_abs_boundary_faces, abs_boundary_ispec, nkmin_xi, nkmin_eta, njmin, njmax, nimin, nimax ]
    if type == :acoustic_forward then
      variables += [ abs_boundary_jacobian2D, wgllwgll, ibool, vpstore, save_forward, b_absorb_potential ]
    elsif type == :acoustic_backward then
      variables += [ ibool ]
    end
    
    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)

    p = Procedure(function_name, variables)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ngllx => n_gllx, :ngll2 => n_gll2 )
      decl p
      decl igll = Int("igll")
      decl iface = Int("iface")
      decl i = Int("i"), j = Int("j"), k = Int("k")
      decl iglob = Int("iglob")
      decl ispec = Int("ispec")
      decl sn = Real("sn") if type == :acoustic_forward
      decl jacobianw = Real("jacobianw") if type == :acoustic_forward
      decl fac1 = Real("fac1") if type == :acoustic_forward
      
      print igll === get_local_id(0)
      print iface === get_group_id(0)+get_group_id(1)*get_num_groups(0)

      print If( iface < num_abs_boundary_faces ) {
        print ispec === abs_boundary_ispec[iface]-1

        print Case( interface_type,
          4, lambda {
            print If( Expression("||", nkmin_xi[INDEX2(2,0,iface)] == 0, njmin[INDEX2(2,0,iface)] == 0) )   { print Return(nil) }
            print i === 0
            print k === igll/ngllx
            print j === igll-k*ngllx
            print If( Expression("||", k < nkmin_xi[INDEX2(2,0,iface)]-1, k > ngllx-1) )                    { print Return(nil) }
            print If( Expression("||", j <    njmin[INDEX2(2,0,iface)]-1, j > njmax[INDEX2(2,0,iface)]-1) ) { print Return(nil) }
            print fac1 === wgllwgll[k*ngllx+j] if type == :acoustic_forward
          },
          5, lambda {
            print If( Expression("||", nkmin_xi[INDEX2(2,1,iface)] == 0, njmin[INDEX2(2,1,iface)] == 0) )   { print Return(nil) }
            print i === ngllx-1
            print k === igll/ngllx
            print j === igll-k*ngllx
            print If( Expression("||", k < nkmin_xi[INDEX2(2,1,iface)]-1, k > ngllx-1) )                    { print Return(nil) }
            print If( Expression("||", j <    njmin[INDEX2(2,1,iface)]-1, j > njmax[INDEX2(2,1,iface)]-1) ) { print Return(nil) }
            print fac1 === wgllwgll[k*ngllx+j] if type == :acoustic_forward
          },
          6, lambda {
            print If( Expression("||", nkmin_eta[INDEX2(2,0,iface)] == 0, nimin[INDEX2(2,0,iface)] == 0) )  { print Return(nil) }
            print j === 0
            print k === igll/ngllx
            print i === igll-k*ngllx
            print If( Expression("||", k < nkmin_eta[INDEX2(2,0,iface)]-1, k > ngllx-1) )                   { print Return(nil) }
            print If( Expression("||", i <     nimin[INDEX2(2,0,iface)]-1, i > nimax[INDEX2(2,0,iface)]-1) ){ print Return(nil) }
            print fac1 === wgllwgll[k*ngllx+i] if type == :acoustic_forward
          },
          7, lambda {
            print If( Expression("||", nkmin_eta[INDEX2(2,1,iface)] == 0, nimin[INDEX2(2,1,iface)] == 0) )  { print Return(nil) }
            print j === ngllx-1
            print k === igll/ngllx
            print i === igll-k*ngllx
            print If( Expression("||", k < nkmin_eta[INDEX2(2,1,iface)]-1, k > ngllx-1) )                   { print Return(nil) }
            print If( Expression("||", i <     nimin[INDEX2(2,1,iface)]-1, i > nimax[INDEX2(2,1,iface)]-1) ){ print Return(nil) }
            print fac1 === wgllwgll[k*ngllx+i] if type == :acoustic_forward
          },
          8, lambda {
            print k === 0
            print j === igll/ngllx
            print i === igll - j*ngllx
            print If( Expression("||", j < 0, j > ngllx-1) )                                                { print Return(nil) }
            print If( Expression("||", i < 0, i > ngllx-1) )                                                { print Return(nil) }
            print fac1 === wgllwgll[j*ngllx+i] if type == :acoustic_forward
          })
      }
      print iglob === ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] - 1
      if type == :acoustic_forward then
        print sn === potential_dot_acoustic[iglob] / vpstore[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)]
        print jacobianw === abs_boundary_jacobian2D[INDEX2(ngll2,igll,iface)]*fac1
        print atomicAdd(potential_dot_dot_acoustic + iglob, -sn*jacobianw)
        print If( save_forward ) {
          print b_absorb_potential[INDEX2(ngll2,igll,iface)] === sn*jacobianw
        }
      elsif type == :acoustic_backward then
        print atomicAdd(b_potential_dot_dot_acoustic + iglob, -b_absorb_potential[INDEX2(ngll2,igll,iface)])
      end
      close p
    else
      raise "Unsupported language!"
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end
end

