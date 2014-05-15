module BOAST
  def BOAST::compute_coupling_ocean_kernel(ref = true, n_dim = 3)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_coupling_ocean_kernel"
    accel_crust_mantle    = Real("accel_crust_mantle",      :dir => :out, :dim => [ Dim(3), Dim() ])
    rmassx_crust_mantle   = Real("rmassx_crust_mantle",     :dir => :in,  :dim => [ Dim() ])
    rmassy_crust_mantle   = Real("rmassy_crust_mantle",     :dir => :in,  :dim => [ Dim() ])
    rmassz_crust_mantle   = Real("rmassz_crust_mantle",     :dir => :in,  :dim => [ Dim() ])
    rmass_crust_mantle = [ rmassx_crust_mantle, rmassy_crust_mantle, rmassz_crust_mantle ]
    rmass_ocean_load      = Real("rmass_ocean_load",        :dir => :in,  :dim => [ Dim() ])
    npoin_ocean_load      = Int( "npoin_ocean_load",        :dir => :in)
    ibool_ocean_load      = Int( "ibool_ocean_load",        :dir => :in,  :dim => [ Dim() ])
    normal_ocean_load     = Real("normal_ocean_load",       :dir => :in,  :dim => [ Dim() ])

    ndim =               Int( "NDIM",              :const => n_dim)

    p = Procedure(function_name, [accel_crust_mantle, rmassx_crust_mantle, rmassy_crust_mantle, rmassz_crust_mantle, rmass_ocean_load, npoin_ocean_load, ibool_ocean_load, normal_ocean_load])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ndim => n_dim )
      decl p
      decl ipoin = Int("ipoin")
      decl iglob = Int("iglob")

      decl *(n               = [Real("nx"), Real("ny"), Real("nz")])
      decl rmass             = Real("rmass")
      decl force_normal_comp = Real("force_normal_comp")
      decl *(additional_term = [Real("additional_term_x"), Real("additional_term_y"), Real("additional_term_z")])

      print ipoin === get_global_id(0) + get_global_size(0)*get_global_id(1)
      print If( ipoin < npoin_ocean_load ) {
        print iglob === ibool_ocean_load[ipoin] - 1
        n.each_index { |indx| print n[indx] === normal_ocean_load[INDEX2(ndim,indx,ipoin)] }

        print force_normal_comp === accel_crust_mantle[0,iglob]*n[0] / rmassx_crust_mantle[iglob] \
                                  + accel_crust_mantle[1,iglob]*n[1] / rmassy_crust_mantle[iglob] \
                                  + accel_crust_mantle[2,iglob]*n[2] / rmassz_crust_mantle[iglob]

        print rmass === rmass_ocean_load[ipoin]

        additional_term.each_index { |indx| print additional_term[indx] === (rmass - rmass_crust_mantle[indx][iglob]) * force_normal_comp }

        additional_term.each_index { |indx| print accel_crust_mantle[indx,iglob] === accel_crust_mantle[indx,iglob] + additional_term[indx] * n[indx] }

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
