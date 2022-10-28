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
      abs_boundary_jacobian2Dw      = Real("abs_boundary_jacobian2Dw",    :dir => :in,   :dim => [ Dim() ])
      vpstore                       = Real("vpstore",                    :dir => :in,   :dim => [ Dim() ])
      save_stacey                   = Int( "SAVE_STACEY",                :dir => :in)
      variables = [potential_dot_acoustic, potential_dot_dot_acoustic]
    elsif type == :acoustic_backward then
      function_name = "compute_stacey_acoustic_backward_kernel"
      b_potential_dot_dot_acoustic  = Real("b_potential_dot_dot_acoustic", :dir => :inout,:dim => [ Dim() ])
      b_absorb_potential            = Real("b_absorb_potential",           :dir => :in,  :dim => [ Dim() ])
      variables = [b_potential_dot_dot_acoustic, b_absorb_potential]
    else
      raise "Unsupported type : #{type}!"
    end

    num_abs_boundary_faces        = Int( "num_abs_boundary_faces",     :dir => :in)
    abs_boundary_ispec            = Int( "abs_boundary_ispec",         :dir => :in,   :dim => [ Dim() ])
    abs_boundary_npoin            = Int( "abs_boundary_npoin",         :dir => :in,   :dim => [ Dim() ])
    abs_boundary_ijk              = Int( "abs_boundary_ijk",           :dir => :in,   :dim => [ Dim() ])
    ibool                         = Int( "ibool",                      :dir => :in,   :dim => [ Dim() ])
    variables += [ num_abs_boundary_faces, abs_boundary_ispec, abs_boundary_npoin, abs_boundary_ijk ]
    if type == :acoustic_forward then
      variables += [ abs_boundary_jacobian2Dw, ibool, vpstore, save_stacey, b_absorb_potential ]
    elsif type == :acoustic_backward then
      variables += [ ibool ]
    end

    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)

    p = Procedure(function_name, variables)
    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      make_specfem3d_header( :ngllx => n_gllx, :ngll2 => n_gll2 )
      open p
      decl npoin = Int("npoin")
      decl igll = Int("igll")
      decl iface = Int("iface")
      decl i = Int("i"), j = Int("j"), k = Int("k")
      decl iglob = Int("iglob")
      decl ispec = Int("ispec")
      decl sn = Real("sn") if type == :acoustic_forward
      decl weight = Real("weight") if type == :acoustic_forward
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

          if type == :acoustic_forward then
            print sn === potential_dot_acoustic[iglob] / vpstore[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)]
            comment()
            print weight === abs_boundary_jacobian2Dw[INDEX2(ngll2,igll,iface)]
            comment()
            print atomicAdd(potential_dot_dot_acoustic + iglob, -sn*weight)
            comment()
            print If(save_stacey) {
              print b_absorb_potential[INDEX2(ngll2,igll,iface)] === sn*weight
            }
          elsif type == :acoustic_backward then
            print atomicAdd(b_potential_dot_dot_acoustic + iglob, -b_absorb_potential[INDEX2(ngll2,igll,iface)])
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

