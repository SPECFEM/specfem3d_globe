
module BOAST

  def BOAST::compute_strength_noise_kernel(ref = true, n_dim = 3, n_gllx = 5, n_gll2 = 25)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_strength_noise_kernel"
    v = []
    v.push displ               = Real("displ",               :dir => :in,    :dim => [Dim(3), Dim()] )
    v.push ibelm_top           = Int( "ibelm_top",           :dir => :in,    :dim => [Dim()] )
    v.push ibool               = Int( "ibool",               :dir => :in,    :dim => [Dim()] )
    v.push noise_surface_movie = Real("noise_surface_movie", :dir => :in,    :dim => [Dim()] )
    v.push *normal_noise       =[Real("normal_x_noise",      :dir => :in,    :dim => [Dim()] ),\
                                 Real("normal_y_noise",      :dir => :in,    :dim => [Dim()] ),\
                                 Real("normal_z_noise",      :dir => :in,    :dim => [Dim()] ) ]
    v.push sigma_kl            = Real("Sigma_kl",            :dir => :inout, :dim => [Dim()] )
    v.push deltat              = Real("deltat",              :dir => :in)
    v.push nspec_top           = Int( "nspec_top",           :dir => :in)

    ndim  = Int("NDIM",  :const => n_dim)
    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)

    p = Procedure(function_name, v)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ndim => n_dim, :ngllx => n_gllx, :ngll2 => n_gll2 )
      decl p
        decl iface = Int("iface")
        decl ispec = Int("ispec")
        decl igll  = Int("igll")
        decl ipoin = Int("ipoin")
        decl i     = Int("i"), j = Int("j"), k = Int("k")
        decl iglob = Int("iglob")
        decl eta   = Real("eta")

        print iface === get_group_id(0) + get_group_id(1)*get_num_groups(0)
        print If( iface < nspec_top ) {
          print ispec === ibelm_top[iface] - 1
          print igll  === get_local_id(0)
          print ipoin === igll + ngll2*iface

          print k === ngllx - 1
          print j === igll/ngllx
          print i === igll - j*ngllx

          print iglob === ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] - 1

          print eta   === noise_surface_movie[INDEX3(ndim,ngll2,0,igll,iface)]*normal_noise[0][ipoin]\
                        + noise_surface_movie[INDEX3(ndim,ngll2,1,igll,iface)]*normal_noise[1][ipoin]\
                        + noise_surface_movie[INDEX3(ndim,ngll2,2,igll,iface)]*normal_noise[2][ipoin]

          print sigma_kl[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] === sigma_kl[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)]\
                                                                  + deltat*eta*( normal_noise[0][ipoin]*displ[0,iglob]\
                                                                               + normal_noise[1][ipoin]*displ[1,iglob]\
                                                                               + normal_noise[2][ipoin]*displ[2,iglob] )
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
