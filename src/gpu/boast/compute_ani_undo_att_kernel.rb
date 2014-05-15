module BOAST

  require './compute_strain_product_helper.rb'
  require './compute_element_strain_undo_att_helper.rb'

  def BOAST::compute_ani_undo_att_kernel(ref = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)
    compute_undo_att_kernel(:ani, ref, n_gllx, n_gll2, n_gll3, n_gll3_padded)
  end

  def BOAST::compute_undo_att_kernel(type, ref = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    if type == :ani then
      function_name = "compute_ani_undo_att_kernel"
    elsif type == :iso then
      function_name = "compute_iso_undo_att_kernel"
    else
      raise "Unsupported type #{type}!"
    end
    v = []
    v.push epsilondev_xx          = Real("epsilondev_xx",        :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yy          = Real("epsilondev_yy",        :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xy          = Real("epsilondev_xy",        :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xz          = Real("epsilondev_xz",        :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yz          = Real("epsilondev_yz",        :dir => :in, :dim => [Dim()] )
    v.push epsilon_trace_over_3   = Real("epsilon_trace_over_3", :dir => :in, :dim => [Dim()] )
    if type == :ani then
      v.push cijkl_kl             = Real("cijkl_kl",             :dir => :inout,:dim => [Dim(21),Dim()] )
    elsif type == :iso then
      v.push mu_kl                = Real("mu_kl",                :dir => :inout,:dim => [Dim()] )
      v.push kappa_kl             = Real("kappa_kl",             :dir => :inout,:dim => [Dim()] )
    end      
    v.push nspec                  = Int( "NSPEC",                :dir => :in)
    v.push deltat                 = Real("deltat",               :dir => :in)
    v.push d_ibool                = Int( "d_ibool",              :dir => :in, :dim => [Dim()] )
    v.push d_b_displ              = Real("d_b_displ",            :dir => :in, :dim => [Dim(3), Dim()] )
    v.push *d_xi =    [d_xix      = Real("d_xix",    :dir => :in, :dim => [Dim()] ), d_xiy    = Real("d_xiy",   :dir => :in, :dim => [Dim()] ), d_xiz    = Real("d_xiz",   :dir => :in, :dim => [Dim()] ) ]
    v.push *d_eta =   [d_etax     = Real("d_etax",                  :dir => :in, :dim => [Dim()] ), d_etay = Real("d_etay",:dir => :in, :dim => [Dim()] ), d_etaz = Real("d_etaz",:dir => :in, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax   = Real("d_gammax",                :dir => :in, :dim => [Dim()] ), d_gammay = Real("d_gammay",:dir => :in, :dim => [Dim()] ), d_gammaz = Real("d_gammaz",:dir => :in, :dim => [Dim()] ) ]
    v.push d_hprime_xx             = Real("d_hprime_xx",             :dir => :in, :dim => [Dim()] )

    epsilondev = [ epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz ]

    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)
    ngll3 = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    p = Procedure(function_name, v)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded )
      sub_compute_element_strain_undo_att = compute_element_strain_undo_att(n_gllx, n_gll2, n_gll3, n_gll3_padded )
      print sub_compute_element_strain_undo_att
      if type == :ani then
        sub_compute_strain_product =  compute_strain_product()
        print sub_compute_strain_product
      end

      decl p

      decl ispec = Int("ispec")
      decl ijk_ispec = Int("ijk_ispec")
      decl tx = Int("tx")
      decl iglob = Int("iglob")

      decl eps_trace_over_3 = Real("eps_trace_over_3")
      decl b_eps_trace_over_3 = Real("b_eps_trace_over_3")
      if type == :ani then 
        decl prod = Real("prod", :dim => [Dim(21)], :allocate => true)
        decl i = Int("i")
      end
      decl epsdev = Real("epsdev", :dim => [Dim(5)], :allocate => true)
      decl b_epsdev = Real("b_epsdev", :dim => [Dim(5)], :allocate => true)

      decl *s_dummy_loc = ["x", "y", "z"].collect { |a|
        Real("s_dummy#{a}_loc", :local => true, :dim => [Dim(ngll3)] )
      }
      decl sh_hprime_xx = Real("sh_hprime_xx",     :local => true, :dim => [Dim(ngll2)] )
      
      print ispec === get_group_id(0) + get_group_id(1)*get_num_groups(0)
      print ijk_ispec === get_local_id(0) + ngll3*ispec
      print tx === get_local_id(0)

      print If(tx < ngll2) {
        print sh_hprime_xx[tx] === d_hprime_xx[tx]
      }

      print If(ispec < nspec) {
        print iglob === d_ibool[ijk_ispec] - 1
        (0..2).each { |indx|
          print s_dummy_loc[indx][tx] === d_b_displ[indx,iglob]
        }
      }
      print barrier(:local)

      print If(ispec < nspec) {
        (0..4).each { |indx|
          print epsdev[indx] === epsilondev[indx][ijk_ispec]
        }

        print eps_trace_over_3 === epsilon_trace_over_3[ijk_ispec]

        print sub_compute_element_strain_undo_att.call(ispec,ijk_ispec,
                                                       d_ibool,
                                                       *s_dummy_loc,
                                                       *d_xi, *d_eta, *d_gamma,
                                                       sh_hprime_xx,
                                                       b_epsdev,b_eps_trace_over_3.address)
        if type == :ani then 
          print sub_compute_strain_product.call(prod, eps_trace_over_3, epsdev, b_eps_trace_over_3, b_epsdev)
          print For(i, 0, 21-1) {
            print cijkl_kl[i, ijk_ispec] === cijkl_kl[i, ijk_ispec] + deltat * prod[i]
          }
        else
          print mu_kl[ijk_ispec] === mu_kl[ijk_ispec] + deltat * \
                                   ( epsdev[0] * b_epsdev[0] + \
                                     epsdev[1] * b_epsdev[1] + \
                                     (epsdev[0] + epsdev[1]) * (b_epsdev[0] + b_epsdev[1]) + \
                                     ( epsdev[2] * b_epsdev[2] + \
                                       epsdev[3] * b_epsdev[3] + \
                                       epsdev[4] * b_epsdev[4]) * 2.0 )

          print kappa_kl[ijk_ispec] === kappa_kl[ijk_ispec] + deltat * ( eps_trace_over_3 * b_eps_trace_over_3 * 9.0 )
        end
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

