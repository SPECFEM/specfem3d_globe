
module BOAST

  def BOAST::compute_iso_kernel(ref = true, n_gll3 = 125)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_iso_kernel"
    v = []
    v.push epsilondev_xx          = Real("epsilondev_xx",          :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yy          = Real("epsilondev_yy",          :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xy          = Real("epsilondev_xy",          :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xz          = Real("epsilondev_xz",          :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yz          = Real("epsilondev_yz",          :dir => :in, :dim => [Dim()] )
    v.push epsilon_trace_over_3   = Real("epsilon_trace_over_3",   :dir => :in, :dim => [Dim()] )
    v.push b_epsilondev_xx        = Real("b_epsilondev_xx",        :dir => :in, :dim => [Dim()] )
    v.push b_epsilondev_yy        = Real("b_epsilondev_yy",        :dir => :in, :dim => [Dim()] )
    v.push b_epsilondev_xy        = Real("b_epsilondev_xy",        :dir => :in, :dim => [Dim()] )
    v.push b_epsilondev_xz        = Real("b_epsilondev_xz",        :dir => :in, :dim => [Dim()] )
    v.push b_epsilondev_yz        = Real("b_epsilondev_yz",        :dir => :in, :dim => [Dim()] )
    v.push b_epsilon_trace_over_3 = Real("b_epsilon_trace_over_3", :dir => :in, :dim => [Dim()] )
    v.push mu_kl                  = Real("mu_kl",                  :dir => :inout,:dim => [Dim()] )
    v.push kappa_kl               = Real("kappa_kl",               :dir => :inout,:dim => [Dim()] )
    v.push nspec                  = Int( "NSPEC",                  :dir => :in)
    v.push deltat                 = Real("deltat",                 :dir => :in)

    ngll3 = Int("NGLL3", :const => n_gll3)

    p = Procedure(function_name, v)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ngll3 => n_gll3 )
      decl p
        decl ispec = Int("ispec")
        decl ijk_ispec = Int("ijk_ispec")
        print ispec === get_group_id(0) + get_group_id(1)*get_num_groups(0)
        print If( ispec < nspec ) {
          print ijk_ispec === get_local_id(0) + ngll3*ispec
          print mu_kl[ijk_ispec] === mu_kl[ijk_ispec] + deltat * ( epsilondev_xx[ijk_ispec]*b_epsilondev_xx[ijk_ispec]\
                                                                 + epsilondev_yy[ijk_ispec]*b_epsilondev_yy[ijk_ispec]\
                                                                 + (    epsilondev_xx[ijk_ispec] +   epsilondev_yy[ijk_ispec] )\
                                                                   * (b_epsilondev_xx[ijk_ispec] + b_epsilondev_yy[ijk_ispec] )\
                                                                 + ( epsilondev_xy[ijk_ispec]*b_epsilondev_xy[ijk_ispec]\
                                                                   + epsilondev_xz[ijk_ispec]*b_epsilondev_xz[ijk_ispec]\
                                                                   + epsilondev_yz[ijk_ispec]*b_epsilondev_yz[ijk_ispec])*2 )
          print kappa_kl[ijk_ispec] === kappa_kl[ijk_ispec]  + deltat * ( epsilon_trace_over_3[ijk_ispec] * b_epsilon_trace_over_3[ijk_ispec] * 9 )
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
