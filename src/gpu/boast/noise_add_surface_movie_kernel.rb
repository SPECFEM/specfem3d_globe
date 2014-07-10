module BOAST

  def BOAST::noise_add_surface_movie_kernel(ref = true, n_dim = 3, n_gllx = 5, n_gll2 = 25)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    function_name = "noise_add_surface_movie_kernel"
 
    accel =               Real("accel",               :dir => :inout,:dim => [Dim()] )
    ibool =               Int( "ibool",               :dir => :in,   :dim => [Dim()] )
    ibelm_top =           Int( "ibelm_top",           :dir => :in,   :dim => [Dim()] )
    nspec_top =           Int( "nspec_top",           :dir => :in)
    noise_surface_movie = Real("noise_surface_movie", :dir => :in,   :dim => [Dim()] )
    normal_x_noise =      Real("normal_x_noise",      :dir => :in,   :dim => [Dim()] )
    normal_y_noise =      Real("normal_y_noise",      :dir => :in,   :dim => [Dim()] )
    normal_z_noise =      Real("normal_z_noise",      :dir => :in,   :dim => [Dim()] )
    mask_noise =          Real("mask_noise",          :dir => :in,   :dim => [Dim()] )
    jacobian2D =          Real("jacobian2D",          :dir => :in,   :dim => [Dim()] )
    wgllwgll   =          Real("wgllwgll",            :dir => :in,   :dim => [Dim()] )

    normal_noise = [normal_x_noise, normal_y_noise, normal_z_noise]
    
    ndim =  Int("NDIM",  :const => n_dim)
    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)

    p = Procedure(function_name, [accel,ibool,ibelm_top,nspec_top,noise_surface_movie,normal_x_noise,normal_y_noise,normal_z_noise,mask_noise,jacobian2D,wgllwgll])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ndim => n_dim, :ngllx => n_gllx, :ngll2 => n_gll2)
      decl p
        decl igll  = Int("igll")
        decl iface  = Int("iface")

        print igll === get_local_id(0)
        print iface === get_group_id(0) + get_group_id(1)*get_num_groups(0)
        print If( iface < nspec_top ) {
          decl i = Int("i")
          decl j = Int("j")
          decl k = Int("k")
          decl ispec = Int("ispec")
          decl iglob = Int("iglob")
          decl ipoin = Int("ipoin")
          decl eta   =Real("eta")
          decl jacobianw = Real("jacobianw")
          decl *(normal = [Real("normal_x"), Real("normal_y"), Real("normal_z")])
          print ispec === ibelm_top[iface] - 1
          print k === ngllx - 1
          print j === igll/ngllx
          print i === igll-j*ngllx
          print iglob === ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)]-1
          print ipoin === ngll2*iface + igll

          normal.each_index { |indx|
            print normal[indx] === normal_noise[indx][ipoin]
          }

          print eta === 0.0
          
          (0..2).each { |indx| 
            print eta === eta + noise_surface_movie[INDEX3(ndim,ngll2,indx,igll,iface)]*normal[indx]
          }

          print jacobianw === wgllwgll[k*ngllx+i]*jacobian2D[igll+ngll2*iface]
          
          (0..2).each { |indx|
            print atomicAdd(accel+ iglob*3 + indx, eta*mask_noise[ipoin]*normal[indx]*jacobianw)
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
