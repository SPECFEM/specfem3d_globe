module BOAST
  def BOAST::compute_seismograms_kernel(ref = false, n_dim = 3, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    v = []
    function_name = "compute_seismograms_kernel"
    v.push nrec_local             = Int("nrec_local",             :dir => :in)
    v.push displ                  = Real("displ",                 :dir => :in,   :dim => [ Dim() ])
    v.push d_ibool                = Int("d_ibool",                :dir => :in,   :dim => [ Dim() ])
    v.push hxir                   = Real("hxir",                  :dir => :in,   :dim => [ Dim() ])
    v.push hetar                  = Real("hetar",                 :dir => :in,   :dim => [ Dim() ])
    v.push hgammar                = Real("hgammar",               :dir => :in,   :dim => [ Dim() ])
    v.push seismograms            = Real("seismograms",           :dir => :inout,:dim => [ Dim() ])
    v.push nu                     = Real("nu",                    :dir => :in,   :dim => [ Dim() ])
    v.push ispec_selected_rec     = Int("ispec_selected_rec",     :dir => :in,   :dim => [ Dim() ])
    v.push number_receiver_global = Int("number_receiver_global", :dir => :in,   :dim => [ Dim() ])
    v.push scale_displ            = Real("scale_displ",           :dir => :in)
    v.push seismo_current         = Int("seismo_current",         :dir => :in)

    ndim         = Int("NDIM",         :const => n_dim)
    ngllx        = Int("NGLLX",        :const => n_gllx)
    ngll2        = Int("NGLL2",        :const => n_gll2)
    ngll3        = Int("NGLL3",        :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    p = Procedure(function_name, v )
    if(get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      make_specfem3d_header( :ndim => n_dim, :ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded )
      open p
      ispec =      Int("ispec")
      iglob =      Int("iglob")
      irec_local = Int("irec_local")
      irec =       Int("irec")
      tx =         Int("tx")
      lagrange=    Real("lagrange")
      i =          Int("i")
      j =          Int("j")
      k =          Int("k")
      l =          Int("l")
      s =          Int("s")
      idx =        Int("idx")
      decl ispec
      decl iglob
      decl irec_local
      decl irec
      decl tx
      decl lagrange
      decl i
      decl j
      decl k
      decl l
      decl s
      decl idx

      decl sh_dxd = Real("sh_dxd",     :local => true, :dim => [Dim(ngll3_padded)] )
      decl sh_dyd = Real("sh_dyd",     :local => true, :dim => [Dim(ngll3_padded)] )
      decl sh_dzd = Real("sh_dzd",     :local => true, :dim => [Dim(ngll3_padded)] )
      comment()

      print tx === get_local_id(0)
      print irec_local === get_group_id(0) + get_num_groups(0)*get_group_id(1)
      comment()

      print k === tx/ngll2
      print j === (tx - k*ngll2)/ngllx
      print i ===  tx - k*ngll2 - j*ngllx
      comment()

      print If (irec_local < nrec_local) {

        print irec  === number_receiver_global[irec_local] - 1
        print ispec === ispec_selected_rec[irec] - 1
        comment()

        print sh_dxd[tx] === 0
        print sh_dyd[tx] === 0
        print sh_dzd[tx] === 0

        print If (tx < ngll3) {
          print lagrange === hxir[irec_local*ngllx + i] * hetar[irec_local*ngllx + j] * hgammar[irec_local*ngllx + k]
          print iglob ===    d_ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)]-1

          print sh_dxd[tx] === lagrange * displ[iglob*3 + 0]
          print sh_dyd[tx] === lagrange * displ[iglob*3 + 1]
          print sh_dzd[tx] === lagrange * displ[iglob*3 + 2]
        }
        print barrier(:local)
        comment()

        print l === 1
        (1..7).each { |indx1|
          print s === l*2
          print If ( BOAST::Modulo(tx,s) == 0) {
            print sh_dxd[tx] === sh_dxd[tx] + sh_dxd[tx + l]
            print sh_dyd[tx] === sh_dyd[tx] + sh_dyd[tx + l]
            print sh_dzd[tx] === sh_dzd[tx] + sh_dzd[tx + l]
          }
          print barrier(:local)
          print l ===  l*2
        }
        comment()

        print idx === INDEX3(ndim,nrec_local,0,irec_local,seismo_current)
        comment()

        print If (tx == 0) {
          print seismograms[idx + 0] === scale_displ * (
            nu[(irec_local*3)*3 + 0]*sh_dxd[0] + nu[(irec_local*3 + 1)*3 + 0]*sh_dyd[0] + nu[(irec_local*3 + 2)*3 + 0]*sh_dzd[0]
            )
       }
        print If (tx == 1) {
          print seismograms[idx + 1] === scale_displ * (
            nu[(irec_local*3)*3 + 1]*sh_dxd[0] + nu[(irec_local*3 + 1)*3 + 1]*sh_dyd[0] + nu[(irec_local*3 + 2)*3 + 1]*sh_dzd[0]
            )
       }
        print If (tx == 2) {
          print seismograms[idx + 2] === scale_displ * (
            nu[(irec_local*3)*3 + 2]*sh_dxd[0] + nu[(irec_local*3 + 1)*3 + 2]*sh_dyd[0] + nu[(irec_local*3 + 2)*3 + 2]*sh_dzd[0]
            )
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
