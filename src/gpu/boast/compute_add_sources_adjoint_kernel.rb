module BOAST
  def BOAST::compute_add_sources_adjoint_kernel(ref = true, n_dim = 3, n_gllx = 5)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_add_sources_adjoint_kernel"
    accel =              Real("accel",             :dir => :inout,:dim => [ Dim() ])
    nrec =               Int( "nrec",              :dir => :in)
    adj_sourcearrays =   Real("adj_sourcearrays",  :dir => :in,   :dim => [ Dim() ])
    ibool =              Int( "ibool",             :dir => :in,   :dim => [ Dim() ])
    ispec_selected_rec = Int( "ispec_selected_rec",:dir => :in,   :dim => [ Dim() ])
    pre_computed_irec =  Int( "pre_computed_irec", :dir => :in,   :dim => [ Dim() ])
    nadj_rec_local =     Int( "nadj_rec_local",    :dir => :in)

    ndim =               Int( "NDIM",              :const => n_dim)
    ngllx =              Int( "NGLLX",             :const => n_gllx)
    p = Procedure(function_name, [accel,nrec,adj_sourcearrays,ibool,ispec_selected_rec,pre_computed_irec,nadj_rec_local])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CUDA or get_lang == CL) then
      make_specfem3d_header( :ndim => n_dim, :ngllx => n_gllx )
      decl p
      ispec =      Int("ispec")
      iglob =      Int("iglob")
      irec_local = Int("irec_local")
      irec =       Int("irec")
      i =          Int("i")
      j =          Int("j")
      k =          Int("k")
      decl ispec
      decl iglob
      decl irec_local
      decl irec
      decl i
      decl j
      decl k

      print irec_local === get_group_id(0) + get_num_groups(0)*get_group_id(1)
      print If(irec_local < nadj_rec_local) {
        print irec === pre_computed_irec[irec_local]
        print ispec === ispec_selected_rec[irec] - 1
        print i === get_local_id(0)
        print j === get_local_id(1)
        print k === get_local_id(2)
        print iglob === ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] - 1
        (0..2).each { |indx|
          print atomicAdd(accel+iglob*3+indx, adj_sourcearrays[INDEX5(ndim,ngllx,ngllx,ngllx,indx,i,j,k,irec_local)])
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
