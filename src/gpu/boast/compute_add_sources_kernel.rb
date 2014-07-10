module BOAST
  def BOAST::compute_add_sources_kernel( ref = true, n_dim = 3, n_gllx = 5 )
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_add_sources_kernel"
    accel =                  Real("accel",                            :dir => :inout,:dim => [ Dim() ] )
    ibool =                  Int("ibool",                             :dir => :in,   :dim => [ Dim() ] )
    sourcearrays =           Real("sourcearrays",                     :dir => :in,   :dim => [ Dim() ] )
    stf_pre_compute =        Real("stf_pre_compute",      :size => 8, :dir => :in,   :dim => [ Dim() ] )
    myrank =                 Int("myrank",                            :dir => :in)
    islice_selected_source = Int("islice_selected_source",            :dir => :in,   :dim => [ Dim() ] )
    ispec_selected_source =  Int("ispec_selected_source",             :dir => :in,   :dim => [ Dim() ] )
    nsources =               Int("NSOURCES",                          :dir => :in)

    ndim =                   Int("NDIM",                  :const => n_dim)
    ngllx =                  Int("NGLLX",                 :const => n_gllx)
    p = Procedure(function_name, [accel,ibool,sourcearrays,stf_pre_compute,myrank,islice_selected_source,ispec_selected_source,nsources])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CUDA or get_lang == CL) then
      make_specfem3d_header( :ndim => n_dim, :ngllx => n_gllx, :double => true )
      decl p
      ispec =   Int( "ispec")
      iglob =   Int( "iglob")
      stf =     Real("stf")
      isource = Int( "isource")
      i =       Int( "i")
      j =       Int( "j")
      k =       Int( "k")
      decl ispec
      decl iglob
      decl stf
      decl isource
      decl i
      decl j
      decl k
      print i === get_local_id(0)
      print j === get_local_id(1)
      print k === get_local_id(2)
      print isource === get_group_id(0) + get_num_groups(0)*get_group_id(1)

      print If(isource < nsources) {
        print If(myrank == islice_selected_source[isource]) {
          print ispec === ispec_selected_source[isource] - 1
          print stf === stf_pre_compute[isource]
          print iglob === ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] - 1
          (0..2).each { |indx|
            print atomicAdd(accel+iglob*3+indx, sourcearrays[INDEX5(ndim,ngllx,ngllx,ngllx,indx,i,j,k,isource)]*stf)
          }
        }
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
