module BOAST

  def BOAST::noise_transfer_surface_to_host_kernel(ref = true, n_dim = 3, n_gllx = 5, n_gll2 = 25)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    function_name = "noise_transfer_surface_to_host_kernel"

    ibelm_top =           Int( "ibelm_top", :dir => :in, :dim => [Dim()] )
    nspec_top =           Int( "nspec_top", :dir => :in)
    ibool =               Int( "ibool",     :dir => :in, :dim => [Dim()] )
    displ =               Real("displ",     :dir => :in, :dim => [Dim()] )
    noise_surface_movie = Real("noise_surface_movie", :dir => :out, :dim => [Dim()] )

    ndim =  Int("NDIM",  :const => n_dim)
    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)

    p = Procedure(function_name, [ibelm_top,nspec_top,ibool,displ,noise_surface_movie])
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
          print ispec === ibelm_top[iface] - 1
          print k === ngllx - 1
          print j === igll/ngllx
          print i === igll-j*ngllx
          print iglob === ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)]-1
          (0..2).each { |indx|
             print noise_surface_movie[INDEX3(ndim,ngll2,indx,igll,iface)] === displ[iglob*3]
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
