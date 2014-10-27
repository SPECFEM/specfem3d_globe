module BOAST

  require './compute_element_strain_undoatt_helper.rb'

  def BOAST::compute_strain_kernel(ref = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_strain_kernel"
    v = []
    v.push d_displ                = Real("d_displ",              :dir => :in, :dim => [Dim(3), Dim()] )
    v.push d_veloc                = Real("d_veloc",              :dir => :in, :dim => [Dim(3), Dim()] )
    v.push epsilondev_xx          = Real("epsilondev_xx",        :dir => :out, :dim => [Dim()] )
    v.push epsilondev_yy          = Real("epsilondev_yy",        :dir => :out, :dim => [Dim()] )
    v.push epsilondev_xy          = Real("epsilondev_xy",        :dir => :out, :dim => [Dim()] )
    v.push epsilondev_xz          = Real("epsilondev_xz",        :dir => :out, :dim => [Dim()] )
    v.push epsilondev_yz          = Real("epsilondev_yz",        :dir => :out, :dim => [Dim()] )
    v.push epsilon_trace_over_3   = Real("epsilon_trace_over_3", :dir => :out, :dim => [Dim()] )
    v.push nspec                  = Int( "NSPEC",                :dir => :in)
    v.push nspec_strain_only      = Int( "NSPEC_STRAIN_ONLY",    :dir => :in)
    v.push deltat                 = Real("deltat",               :dir => :in)
    v.push d_ibool                = Int( "d_ibool",              :dir => :in, :dim => [Dim()] )
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
    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif (get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded )

      sub_compute_element_strain_undoatt = compute_element_strain_undoatt(n_gllx, n_gll2, n_gll3, n_gll3_padded )
      print sub_compute_element_strain_undoatt

      open p

      decl ispec = Int("ispec")
      decl ijk_ispec = Int("ijk_ispec")
      decl tx = Int("tx")
      decl iglob = Int("iglob")

      decl eps_trace_over_3 = Real("eps_trace_over_3")
      decl epsdev = Real("epsdev", :dim => [Dim(5)], :allocate => true)

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
          print s_dummy_loc[indx][tx] === d_displ[indx,iglob] + deltat * d_veloc[indx,iglob]
        }
      }
      print barrier(:local)

      print If(ispec < nspec) {

        print sub_compute_element_strain_undoatt.call(ispec,ijk_ispec,
                                                      d_ibool,
                                                      *s_dummy_loc,
                                                      *d_xi, *d_eta, *d_gamma,
                                                      sh_hprime_xx,
                                                      epsdev,eps_trace_over_3.address)

        print If(nspec_strain_only == 1, lambda {
          print epsilon_trace_over_3[tx] === eps_trace_over_3
        }, lambda {
          print epsilon_trace_over_3[ijk_ispec] === eps_trace_over_3
        })
        
        (0..4).each { |indx|
          print epsilondev[indx][ijk_ispec] === epsdev[indx]
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

