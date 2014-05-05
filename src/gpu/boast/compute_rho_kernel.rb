
module BOAST

  def BOAST::compute_rho_kernel(ref = true, n_gll3 = 125)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_rho_kernel"
    ibool   = Int( "ibool",   :dir => :in, :dim => [Dim()] )
    accel   = Real("accel",   :dir => :in, :dim => [Dim(3), Dim()] )
    b_displ = Real("b_displ", :dir => :in, :dim => [Dim(3), Dim()] )
    rho_kl  = Real("rho_kl",  :dir => :inout,:dim => [Dim()] )
    nspec   = Int( "NSPEC",   :dir => :in)
    deltat  = Real("deltat",  :dir => :in)

    ngll3 = Int("NGLL3", :const => n_gll3)

    p = Procedure(function_name, [ibool, accel, b_displ, rho_kl, nspec, deltat])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ngll3 => n_gll3 )
      decl p
        decl ispec = Int("ispec")
        decl ijk_ispec = Int("ijk_ispec")
        decl iglob = Int("iglob")
        print ispec === get_group_id(0) + get_group_id(1)*get_num_groups(0)
        print If( ispec < nspec ) {
          print ijk_ispec === get_local_id(0) + ngll3*ispec
          print iglob === ibool[ijk_ispec] - 1
          print rho_kl[ijk_ispec] === rho_kl[ijk_ispec] + deltat * ( accel[0, iglob] * b_displ[0, iglob]\
                                                                   + accel[1, iglob] * b_displ[1, iglob]\
                                                                   + accel[2, iglob] * b_displ[2, iglob])
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
