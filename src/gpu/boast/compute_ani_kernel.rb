module BOAST

  require './compute_strain_product_helper.rb'

  def BOAST::compute_ani_kernel(ref = true, n_gll3 = 125)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_ani_kernel"
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
    v.push cijkl_kl               = Real("cijkl_kl",               :dir => :inout,:dim => [Dim(21),Dim()] )
    v.push nspec                  = Int( "NSPEC",                  :dir => :in)
    v.push deltat                 = Real("deltat",                 :dir => :in)

    epsilondev = [ epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz ]
    b_epsilondev = [ b_epsilondev_xx, b_epsilondev_yy, b_epsilondev_xy, b_epsilondev_xz, b_epsilondev_yz ]

    ngll3 = Int("NGLL3", :const => n_gll3)

    p = Procedure(function_name, v)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ngll3 => n_gll3 )
      sub_compute_strain_product =  compute_strain_product()
      print sub_compute_strain_product
      decl p
        decl i = Int("i")
        decl ispec = Int("ispec")
        decl ijk_ispec = Int("ijk_ispec")
        decl eps_trace_over_3 = Real("eps_trace_over_3")
        decl b_eps_trace_over_3 = Real("b_eps_trace_over_3")
        decl prod = Real("prod", :dim => [Dim(21)], :allocate => true)
        decl epsdev = Real("epsdev", :dim => [Dim(5)], :allocate => true)
        decl b_epsdev = Real("b_epsdev", :dim => [Dim(5)], :allocate => true)


        print ispec === get_group_id(0) + get_group_id(1)*get_num_groups(0)
        print If( ispec < nspec ) {
          print ijk_ispec === get_local_id(0) + ngll3*ispec
          (0..4).each { |indx|
            print epsdev[indx] === epsilondev[indx][ijk_ispec]
          }
          (0..4).each { |indx|
            print epsdev[indx] === b_epsilondev[indx][ijk_ispec]
          }
          print eps_trace_over_3 === epsilon_trace_over_3[ijk_ispec]
          print b_eps_trace_over_3 === b_epsilon_trace_over_3[ijk_ispec]
          print sub_compute_strain_product.call(prod,eps_trace_over_3,epsdev,b_eps_trace_over_3,b_epsdev)

          print For(i,0,20) {
            print cijkl_kl[i, ijk_ispec] === cijkl_kl[i, ijk_ispec] + deltat * prod[i]
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
