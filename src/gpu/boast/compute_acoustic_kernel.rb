module BOAST
  def BOAST::compute_gradient_kernel(n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)
    v = []
    function_name = "compute_gradient_kernel"
    v.push ijk                  = Int("ijk", :dir => :in)
    v.push ispec                = Int( "ispec",                :dir => :in)
    v.push scalar_field         = Real("scalar_field",         :dir => :in, :dim => [Dim()], :local => true)
    v.push vector_field_element = Real("vector_field_element", :dir => :out,:dim => [Dim(3)],:private => true)
    v.push hprime_xx            = Real("hprime_xx",            :dir => :in, :dim => [Dim()] )
    v.push *d_xi =    [d_xix    = Real("d_xix",    :dir => :in, :dim => [Dim()] ), d_xiy    = Real("d_xiy",   :dir => :in, :dim => [Dim()] ), d_xiz    = Real("d_xiz",   :dir => :in, :dim => [Dim()] ) ]
    v.push *d_eta =   [d_etax   = Real("d_etax",   :dir => :in, :dim => [Dim()] ), d_etay   = Real("d_etay",  :dir => :in, :dim => [Dim()] ), d_etaz   = Real("d_etaz",  :dir => :in, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax = Real("d_gammax", :dir => :in, :dim => [Dim()] ), d_gammay = Real("d_gammay",:dir => :in, :dim => [Dim()] ), d_gammaz = Real("d_gammaz",:dir => :in, :dim => [Dim()] ) ]

    ngllx        = Int("NGLLX",        :const => n_gllx)
    ngll2        = Int("NGLL2",        :const => n_gll2)
    ngll3        = Int("NGLL3",        :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    p = Procedure(function_name, v, [],:local => true) {
      decl *templ    = [Real("temp1l"),    Real("temp2l"),    Real("temp3l")]
      decl *hp       = [Real("hp1"),       Real("hp2"),       Real("hp3")]
      decl *xil      = [Real("xixl"),      Real("xiyl"),      Real("xizl")]
      decl *etal     = [Real("etaxl"),     Real("etayl"),     Real("etazl")]
      decl *gammal   = [Real("gammaxl"),   Real("gammayl"),   Real("gammazl")]
      decl l = Int("l")
      decl offset   = Int("offset")
      decl *offsets = [Int("offset1"), Int("offset2"), Int("offset3")]
      decl i = Int("I"), j = Int("J"), k = Int("K")

      print k === ijk/ngll2
      print j === (ijk - k*ngll2)/ngllx
      print i ===  ijk - k*ngll2 - j*ngllx

      (0..2).each { |indx|
        print templ[indx] === 0.0
      }
      print For(l, 0, ngllx-1) {
        print hp[0] === hprime_xx[l*ngllx + i]
        print hp[1] === hprime_xx[l*ngllx + j]
        print hp[2] === hprime_xx[l*ngllx + k]
        print offsets[0] === k*ngll2 + j*ngllx + l
        print offsets[1] === k*ngll2 + l*ngllx + i
        print offsets[2] === l*ngll2 + j*ngllx + i
        (0..2).each { |indx|
          print templ[indx] === templ[indx] + scalar_field[offsets[indx]]*hp[indx]
        }
      }

      print offset === ispec*ngll3_padded + ijk
      (0..2).each { |indx|
        print xil[indx] === d_xi[indx][offset]
      }
      (0..2).each { |indx|
        print etal[indx] === d_eta[indx][offset]
      }
      (0..2).each { |indx|
        print gammal[indx] === d_gamma[indx][offset]
      }

      (0..2).each { |indx|
        print vector_field_element[indx] === templ[0]*xil[indx] + templ[1]*etal[indx] + templ[2]*gammal[indx]
      }
    }
    return p
  end

  def BOAST::compute_acoustic_kernel(ref = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    v = []
    function_name = "compute_acoustic_kernel"
    v.push ibool      = Int( "ibool",     :dir => :in, :dim => [Dim()] )
    v.push rhostore   = Real("rhostore",  :dir => :in, :dim => [Dim()] )
    v.push kappastore = Real("kappastore",:dir => :in, :dim => [Dim()] )
    v.push hprime_xx  = Real("hprime_xx", :dir => :in, :dim => [Dim()] )
    v.push *d_xi =    [d_xix    = Real("d_xix",    :dir => :in, :dim => [Dim()] ), d_xiy    = Real("d_xiy",   :dir => :in, :dim => [Dim()] ), d_xiz    = Real("d_xiz",   :dir => :in, :dim => [Dim()] ) ]
    v.push *d_eta =   [d_etax   = Real("d_etax",   :dir => :in, :dim => [Dim()] ), d_etay   = Real("d_etay",  :dir => :in, :dim => [Dim()] ), d_etaz   = Real("d_etaz",  :dir => :in, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax = Real("d_gammax", :dir => :in, :dim => [Dim()] ), d_gammay = Real("d_gammay",:dir => :in, :dim => [Dim()] ), d_gammaz = Real("d_gammaz",:dir => :in, :dim => [Dim()] ) ]
    v.push potential_dot_dot_acoustic   = Real("potential_dot_dot_acoustic",  :dir => :in, :dim => [ Dim() ] )
    v.push b_potential_acoustic         = Real("b_potential_acoustic",        :dir => :in, :dim => [ Dim() ] )
    v.push b_potential_dot_dot_acoustic = Real("b_potential_dot_dot_acoustic",:dir => :in, :dim => [ Dim() ] )
    v.push rho_ac_kl   = Real("rho_ac_kl",   :dir => :inout,:dim => [Dim()] )
    v.push kappa_ac_kl = Real("kappa_ac_kl", :dir => :inout,:dim => [Dim()] )
    v.push deltat      = Real("deltat",      :dir => :in)
    v.push nspec       = Int( "NSPEC",       :dir => :in)

    ngll3        = Int("NGLL3",        :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    p = Procedure(function_name, v)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded )
      sub_compute_gradient_kernel = compute_gradient_kernel(n_gllx, n_gll2, n_gll3, n_gll3_padded)
      print sub_compute_gradient_kernel
      decl p
        decl ispec            = Int("ispec")
        decl ijk              = Int("ijk")
        decl ijk_ispec        = Int("ijk_ispec")
        decl ijk_ispec_padded = Int("ijk_ispec_padded")
        decl iglob            = Int("iglob")
        decl accel_elm   = Real("accel_elm",   :dim => [Dim(3)], :allocate => true)
        decl b_displ_elm = Real("b_displ_elm", :dim => [Dim(3)], :allocate => true)
        decl rhol        = Real("rhol")
        decl kappal      = Real("kappal")
        decl div_displ   = Real("div_displ")
        decl b_div_displ = Real("b_div_displ")
        decl scalar_field_displ = Real("scalar_field_displ", :local => true, :dim => [Dim(ngll3)])
        decl scalar_field_accel = Real("scalar_field_accel", :local => true, :dim => [Dim(ngll3)])

        print ispec === get_group_id(0) + get_group_id(1)*get_num_groups(0)
        print If( ispec < nspec ) {
          print ijk === get_local_id(0)
          print ijk_ispec        === ijk + ngll3       *ispec
          print ijk_ispec_padded === ijk + ngll3_padded*ispec
          print iglob === ibool[ijk_ispec] - 1

          print scalar_field_displ[ijk] === b_potential_acoustic[iglob]
          print scalar_field_accel[ijk] === potential_dot_dot_acoustic[iglob]
          print barrier(:local)
          print sub_compute_gradient_kernel.call(ijk, ispec,\
                                                 scalar_field_displ, b_displ_elm,\
                                                 hprime_xx,\
                                                 d_xix,    d_xiy,    d_xiz,\
                                                 d_etax,   d_etay,   d_etaz,\
                                                 d_gammax, d_gammay, d_gammaz)
          print sub_compute_gradient_kernel.call(ijk, ispec,\
                                                 scalar_field_accel, accel_elm,\
                                                 hprime_xx,\
                                                 d_xix,    d_xiy,    d_xiz,\
                                                 d_etax,   d_etay,   d_etaz,\
                                                 d_gammax, d_gammay, d_gammaz)
          print rhol   === rhostore[ijk_ispec_padded]
          print rho_ac_kl[ijk_ispec] === rho_ac_kl[ijk_ispec] + deltat * rhol * ( accel_elm[0]*b_displ_elm[0]\
                                                                                + accel_elm[1]*b_displ_elm[1]\
                                                                                + accel_elm[2]*b_displ_elm[2] )
          print kappal === rhol / kappastore[ijk_ispec_padded]
          print div_displ   === kappal * potential_dot_dot_acoustic[iglob]
          print b_div_displ === kappal * b_potential_dot_dot_acoustic[iglob]
          print kappa_ac_kl[ijk_ispec] === kappa_ac_kl[ijk_ispec] + deltat * div_displ * b_div_displ
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
