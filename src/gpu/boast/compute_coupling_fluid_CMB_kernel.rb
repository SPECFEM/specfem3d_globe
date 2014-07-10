module BOAST
  def BOAST::compute_coupling_fluid_CMB_kernel(ref = true, n_dim = 3, n_gllx = 5)
    BOAST::compute_coupling_kernel(ref, :fluid_CMB, n_dim, n_gllx)
  end
 
  def BOAST::compute_coupling_kernel(ref = true, type = :fluid_CMB, n_dim = 3, n_gllx = 5)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    if type == :fluid_CMB then
      function_name = "compute_coupling_fluid_CMB_kernel"
      displ                 = Real("displ_crust_mantle",        :dir => :in,    :dim => [ Dim() ])
      accel_out             = Real("accel_outer_core",          :dir => :inout, :dim => [ Dim() ])
      ibool_1               = Int( "ibool_crust_mantle",        :dir => :in,    :dim => [ Dim() ])
      ibelm_1               = Int( "ibelm_bottom_crust_mantle", :dir => :in,    :dim => [ Dim() ])
      normal_outer_core     = Real("normal_top_outer_core",     :dir => :in,    :dim => [ Dim() ])
      jacobian2D_outer_core = Real("jacobian2D_top_outer_core", :dir => :in,    :dim => [ Dim() ])
      ibool_2               = Int( "ibool_outer_core",          :dir => :in,    :dim => [ Dim() ])
      ibelm_2               = Int( "ibelm_top_outer_core",      :dir => :in,    :dim => [ Dim() ])
      nspec2D               = Int( "NSPEC2D_TOP_OC",            :dir => :in)
    elsif type == :fluid_ICB then
      function_name = "compute_coupling_fluid_ICB_kernel"
      displ                 = Real("displ_inner_core",           :dir => :in,    :dim => [ Dim() ])
      accel_out             = Real("accel_outer_core",           :dir => :inout, :dim => [ Dim() ])
      ibool_1               = Int( "ibool_inner_core",           :dir => :in,    :dim => [ Dim() ])
      ibelm_1               = Int( "ibelm_top_inner_core",       :dir => :in,    :dim => [ Dim() ])
      normal_outer_core     = Real("normal_bottom_outer_core",   :dir => :in,    :dim => [ Dim() ])
      jacobian2D_outer_core = Real("jacobian2D_bottom_outer_core", :dir => :in,  :dim => [ Dim() ])
      ibool_2               = Int( "ibool_outer_core",           :dir => :in,    :dim => [ Dim() ])
      ibelm_2               = Int( "ibelm_bottom_outer_core",    :dir => :in,    :dim => [ Dim() ])
      nspec2D               = Int( "NSPEC2D_BOTTOM_OC",          :dir => :in)
    elsif type == :CMB_fluid then
      function_name = "compute_coupling_CMB_fluid_kernel"
      displ                 = Real("displ_crust_mantle",        :dir => :in,    :dim => [ Dim() ])
      accel_out             = Real("accel_crust_mantle",        :dir => :inout, :dim => [ Dim() ])
      accel_in              = Real("accel_outer_core",          :dir => :in,    :dim => [ Dim() ])
      ibool_2               = Int( "ibool_crust_mantle",        :dir => :in,    :dim => [ Dim() ])
      ibelm_2               = Int( "ibelm_bottom_crust_mantle", :dir => :in,    :dim => [ Dim() ])
      normal_outer_core     = Real("normal_top_outer_core",     :dir => :in,    :dim => [ Dim() ])
      jacobian2D_outer_core = Real("jacobian2D_top_outer_core", :dir => :in,    :dim => [ Dim() ])
      ibool_1               = Int( "ibool_outer_core",          :dir => :in,    :dim => [ Dim() ])
      ibelm_1               = Int( "ibelm_top_outer_core",      :dir => :in,    :dim => [ Dim() ])
      rho_oc                = Real("RHO_TOP_OC",                :dir => :in)
      minus_g               = Real("minus_g_cmb",               :dir => :in)
      gravity               = Int( "GRAVITY")
      nspec2D               = Int( "NSPEC2D_BOTTOM_CM",         :dir => :in)
    elsif type == :ICB_fluid then
      function_name = "compute_coupling_ICB_fluid_kernel"
      displ                 = Real("displ_inner_core",          :dir => :in,    :dim => [ Dim() ])
      accel_out             = Real("accel_inner_core",          :dir => :inout, :dim => [ Dim() ])
      accel_in              = Real("accel_outer_core",          :dir => :in,    :dim => [ Dim() ])
      ibool_2               = Int( "ibool_inner_core",          :dir => :in,    :dim => [ Dim() ])
      ibelm_2               = Int( "ibelm_top_inner_core",      :dir => :in,    :dim => [ Dim() ])
      normal_outer_core     = Real("normal_bottom_outer_core",  :dir => :in,    :dim => [ Dim() ])
      jacobian2D_outer_core = Real("jacobian2D_bottom_outer_core",:dir => :in,  :dim => [ Dim() ])
      ibool_1               = Int( "ibool_outer_core",          :dir => :in,    :dim => [ Dim() ])
      ibelm_1               = Int( "ibelm_bottom_outer_core",   :dir => :in,    :dim => [ Dim() ])
      rho_oc                = Real("RHO_BOTTOM_OC",             :dir => :in)
      minus_g               = Real("minus_g_icb",               :dir => :in)
      gravity               = Int( "GRAVITY")
      nspec2D               = Int( "NSPEC2D_TOP_IC",            :dir => :in)
    else
     raise "Unsupported coupling_fluid_type!"
    end
    wgllwgll_xy               = Real("wgllwgll_xy",             :dir => :in,    :dim => [ Dim() ])

    ndim =               Int( "NDIM",              :const => n_dim)
    ngllx =              Int( "NGLLX",             :const => n_gllx)
    if type == :fluid_ICB or type == :fluid_CMB then
      variables = [displ, accel_out, ibool_1, ibelm_1, normal_outer_core, jacobian2D_outer_core, wgllwgll_xy, ibool_2, ibelm_2, nspec2D]
    elsif type == :CMB_fluid or type == :ICB_fluid then
      variables = [displ, accel_out, accel_in, ibool_2, ibelm_2, normal_outer_core, jacobian2D_outer_core, wgllwgll_xy, ibool_1, ibelm_1, rho_oc, minus_g, gravity, nspec2D]
    end
    p = Procedure(function_name, variables)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ndim => n_dim, :ngllx => n_gllx )
      decl p
      decl i = Int("i"), j = Int("j"), k = Int("k")
      decl iface =      Int("iface")
      decl k_corresp =  Int("k_corresp")
      if type == :CMB_fluid then
        decl iglob_1 =  Int("iglob_oc")
        decl iglob_2 =  Int("iglob_cm")
        decl pressure = Real("pressure")
      elsif type == :ICB_fluid then
        decl iglob_1 =  Int("iglob_oc")
        decl iglob_2 =  Int("iglob_ic")
        decl pressure = Real("pressure")
      elsif type == :fluid_CMB then
        decl displ_n = Real("displ_n")
        decl iglob_1 =  Int("iglob_cm")
        decl iglob_2 =  Int("iglob_oc")
      elsif type == :fluid_ICB then
        decl displ_n = Real("displ_n")
        decl iglob_1 =  Int("iglob_ic")
        decl iglob_2 =  Int("iglob_oc")
      end
      decl ispec =      Int("ispec")
      decl ispec_selected = Int("ispec_selected")
      decl *(displ_a = [Real("displ_x"), Real("displ_y"), Real("displ_z")]) unless type == :CMB_fluid or type == :ICB_fluid
      decl *(n = [Real("nx"), Real("ny"), Real("nz")])
      decl weight = Real("weight")

      print i === get_local_id(0)
      print j === get_local_id(1)
      print iface === get_group_id(0) + get_num_groups(0)*get_group_id(1)
      print If( iface < nspec2D ) {
        print ispec === ibelm_2[iface] - 1
        print ispec_selected === ibelm_1[iface] - 1
        if type == :fluid_CMB or type == :ICB_fluid then
          print k === ngllx - 1
          print k_corresp === 0
        elsif type == :fluid_ICB or type == :CMB_fluid then
          print k === 0
          print k_corresp === ngllx - 1
        end
        print iglob_1 === ibool_1[INDEX4(ngllx,ngllx,ngllx,i,j,k_corresp,ispec_selected)] - 1
        displ_a.each_index { |indx| print displ_a[indx] === displ[iglob_1*3+indx] } unless type == :CMB_fluid or type == :ICB_fluid
        n.each_index { |indx| print n[indx] === normal_outer_core[INDEX4(ndim,ngllx,ngllx,indx,i,j,iface)] }
        print displ_n === displ_a[0]*n[0] + displ_a[1]*n[1] + displ_a[2]*n[2] unless type == :CMB_fluid or type == :ICB_fluid
        print weight === jacobian2D_outer_core[INDEX3(ngllx,ngllx,i,j,iface)]*wgllwgll_xy[INDEX2(ngllx,i,j)]
        print iglob_2 === ibool_2[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] - 1
        if type == :fluid_CMB then
          print atomicAdd(accel_out+iglob_2, weight*displ_n)
        elsif type == :fluid_ICB then
          print atomicAdd(accel_out+iglob_2, -weight*displ_n)
        else
          print If( gravity, lambda {
            print pressure === rho_oc * ( minus_g * ( displ[iglob_2*3]*n[0] + displ[iglob_2*3+1]*n[1] + displ[iglob_2*3+2]*n[2] ) - accel_in[iglob_1] )
          } , nil, lambda {
            print pressure === -rho_oc * accel_in[iglob_1]
          } )
          if type == :CMB_fluid then
            n.each_index { |indx| print atomicAdd(accel_out+iglob_2*3+indx, weight*n[indx]*pressure) }
          else
            n.each_index { |indx| print atomicAdd(accel_out+iglob_2*3+indx, -weight*n[indx]*pressure) }
          end
        end
      }
      close p
    else
      raise "Unsupported language!"
    end
    pop_env(:array_start)
    kernel.procedure = p
    return kernel
  end
end
