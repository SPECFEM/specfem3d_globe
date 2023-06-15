module BOAST
  require "./compute_element_att_stress_helper.rb"
  def BOAST::compute_element_ic_att_stress( n_gll3 = 125, n_sls = 3, use_cuda_shared_async = false )
    return BOAST::compute_element_att_stress( :inner_core, n_gll3, n_sls)
  end

  def BOAST::compute_element_cm_att_stress( n_gll3 = 125, n_sls = 3, use_cuda_shared_async = false )
    return BOAST::compute_element_att_stress( :crust_mantle, n_gll3, n_sls, use_cuda_shared_async)
  end

  require "./compute_element_att_memory_helper.rb"
  def BOAST::compute_element_ic_att_memory(n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, use_cuda_shared_async = false )
    return BOAST::compute_element_att_memory(:inner_core, n_gll3, n_gll3_padded, n_sls)
  end

  def BOAST::compute_element_ic_att_memory_lddrk(n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, use_cuda_shared_async = false )
    return BOAST::compute_element_att_memory(:inner_core_lddrk, n_gll3, n_gll3_padded, n_sls)
  end

  def BOAST::compute_element_cm_att_memory( n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, use_cuda_shared_async = false )
    return compute_element_att_memory( :crust_mantle, n_gll3, n_gll3_padded, n_sls, use_cuda_shared_async )
  end

  def BOAST::compute_element_cm_att_memory_lddrk( n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, use_cuda_shared_async = false )
    return compute_element_att_memory( :crust_mantle_lddrk, n_gll3, n_gll3_padded, n_sls, use_cuda_shared_async )
  end

  require "./compute_element_gravity_helper.rb"
  def BOAST::compute_element_ic_gravity( n_gll3 = 125, use_cuda_shared_async = false )
    return BOAST::compute_element_gravity( :inner_core, n_gll3 )
  end

  def BOAST::compute_element_cm_gravity( n_gll3 = 125, use_cuda_shared_async = false )
    return compute_element_gravity( :crust_mantle, n_gll3, use_cuda_shared_async )
  end

#----------------------------------------------------------------------
##
## defines compute element subroutines
##
#----------------------------------------------------------------------

## crust/mantle

  # cuda asynchronuous copies modification
  # to allow the two expressions as function calls below
  #if use_cuda_shared_sync then
  #  register_funccall("compute_offset_sh")
  #  register_funccall("compute_offset")
  #end

  # CUDA_SHARED_ASYNC modifications
  def BOAST::compute_offset(n_gll3 = 125, n_sls = 3)
    function_name = "compute_offset"
    v = []
    v.push tx = Int( "tx", :dir => :in)
    v.push i_sls = Int( "i_sls", :dir => :in)
    v.push working_element = Int( "working_element", :dir => :in)
    ngll3 = Int("NGLL3", :const => n_gll3)
    nsls  = Int("N_SLS", :const => n_sls)

    offset = Int( "offset", :const => tx + ngll3 * (i_sls + nsls * working_element) )
    p = Procedure( function_name, v, :local => true, :qualifiers => "__forceinline__", :return => offset){
      get_output.puts "  // global memory offset"
    }
    return p
  end

  def BOAST::compute_offset_sh(n_gll3 = 125)
    function_name = "compute_offset_sh"
    v = []
    v.push tx = Int( "tx", :dir => :in)
    v.push i_sls = Int( "i_sls", :dir => :in)
    ngll3 = Int("NGLL3", :const => n_gll3)

    offset = Int( "offset", :const => tx + i_sls * ngll3)
    p = Procedure( function_name, v, :local => true, :qualifiers => "__forceinline__", :return => offset){
      get_output.puts "  // shared memory offset"
    }
    return p
  end

  # element routines
  def BOAST::compute_element_cm_aniso(n_gll3 = 125, use_cuda_shared_async = false)
    function_name = "compute_element_cm_aniso"
    v = []
    ngll3 = Int("NGLL3", :const => n_gll3)
    if use_cuda_shared_async then
      v.push tx                     = Int( "tx",                 :dir => :in)
      v.push sh_cstore              = Real("sh_cstore",          :dir => :in, :dim => [Dim(ngll3), Dim()])
    else
      v.push offset                 = Int( "offset", :dir => :in)
      d_cstore = (0..5).collect { |indx1|
                  (0..5).collect { |indx2|
                    if (indx2 < indx1) then
                      nil
                    else
                      Real( "d_c#{indx1+1}#{indx2+1}store", :dir => :in, :dim => [Dim()] )
                    end
                  }
                }
      v.push *(d_cstore.flatten.reject { |e| e.nil? })
    end
    dudl = ["x", "y", "z"].collect { |a1|
      ["x", "y", "z"].collect { |a2|
        Real("du#{a1}d#{a2}l", :dir => :in)
      }
    }
    # elements outside the diagonal are unused
    v.push *(dudl.flatten)
    v.push duxdyl_plus_duydxl = Real("duxdyl_plus_duydxl", :dir => :in)
    v.push duzdxl_plus_duxdzl = Real("duzdxl_plus_duxdzl", :dir => :in)
    v.push duzdyl_plus_duydzl = Real("duzdyl_plus_duydzl", :dir => :in)
    v.push sigma_xx = Real("sigma_xx", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yy = Real("sigma_yy", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_zz = Real("sigma_zz", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xy = Real("sigma_xy", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xz = Real("sigma_xz", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yz = Real("sigma_yz", :dir => :out, :dim => [Dim()], :private => true )

    p = Procedure( function_name, v, :local => true ) {
      c = (0..5).collect { |indx1|
            (0..5).collect { |indx2|
              if (indx2 < indx1) then
                nil
              else
                Real( "c#{indx1+1}#{indx2+1}" )
              end
            }
          }
      decl *(c.flatten.reject { |e| e.nil? })
      comment()

      if use_cuda_shared_async then
        comment("// CUDA asynchronuous memory copies")
        offset = 0
        (0..5).each { |indx1|
          (0..5).each { |indx2|
            print c[indx1][indx2] === sh_cstore[tx,offset] unless indx2 < indx1
            offset += 1 unless indx2 < indx1
          }
        }
      else
        (0..5).each { |indx1|
          (0..5).each { |indx2|
            print c[indx1][indx2] === d_cstore[indx1][indx2][offset] unless indx2 < indx1
          }
        }
      end
      comment()

      # not needed anymore, scaling will be done in prepare_elastic_elements.F90
      #print If(attenuation) {
      #  print minus_sum_beta === one_minus_sum_beta_use - 1.0
      #  print mul === c[3][3] * minus_sum_beta
      #
      #  print c[0][0] === c[0][0] + mul * 1.33333333333333333333
      #  print c[0][1] === c[0][1] - mul * 0.66666666666666666666
      #  print c[0][2] === c[0][2] - mul * 0.66666666666666666666
      #  print c[1][1] === c[1][1] + mul * 1.33333333333333333333
      #  print c[1][2] === c[1][2] - mul * 0.66666666666666666666
      #  print c[2][2] === c[2][2] + mul * 1.33333333333333333333
      #  print c[3][3] === c[3][3] + mul
      #  print c[4][4] === c[4][4] + mul
      #  print c[5][5] === c[5][5] + mul
      #}
      #comment()

      print sigma_xx.dereference === c[0][0]*dudl[0][0] + c[0][5]*duxdyl_plus_duydxl + c[0][1]*dudl[1][1] \
                                   + c[0][4]*duzdxl_plus_duxdzl + c[0][3]*duzdyl_plus_duydzl + c[0][2]*dudl[2][2]
      print sigma_yy.dereference === c[0][1]*dudl[0][0] + c[1][5]*duxdyl_plus_duydxl + c[1][1]*dudl[1][1] \
                                   + c[1][4]*duzdxl_plus_duxdzl + c[1][3]*duzdyl_plus_duydzl + c[1][2]*dudl[2][2]
      print sigma_zz.dereference === c[0][2]*dudl[0][0] + c[2][5]*duxdyl_plus_duydxl + c[1][2]*dudl[1][1] \
                                   + c[2][4]*duzdxl_plus_duxdzl + c[2][3]*duzdyl_plus_duydzl + c[2][2]*dudl[2][2]
      print sigma_xy.dereference === c[0][5]*dudl[0][0] + c[5][5]*duxdyl_plus_duydxl + c[1][5]*dudl[1][1] \
                                   + c[4][5]*duzdxl_plus_duxdzl + c[3][5]*duzdyl_plus_duydzl + c[2][5]*dudl[2][2]
      print sigma_xz.dereference === c[0][4]*dudl[0][0] + c[4][5]*duxdyl_plus_duydxl + c[1][4]*dudl[1][1] \
                                   + c[4][4]*duzdxl_plus_duxdzl + c[3][4]*duzdyl_plus_duydzl + c[2][4]*dudl[2][2]
      print sigma_yz.dereference === c[0][3]*dudl[0][0] + c[3][5]*duxdyl_plus_duydxl + c[1][3]*dudl[1][1] \
                                   + c[3][4]*duzdxl_plus_duxdzl + c[3][3]*duzdyl_plus_duydzl + c[2][3]*dudl[2][2]
    }
    return p
  end

  def BOAST::compute_element_cm_iso(n_gll3 = 125, use_cuda_shared_async = false)
    function_name = "compute_element_cm_iso"
    v = []
    ngll3 = Int("NGLL3", :const => n_gll3)

    if use_cuda_shared_async then
      v.push tx                   = Int( "tx",                 :dir => :in)
      v.push sh_mul               = Real("sh_mul",             :dir => :in, :dim => [Dim(ngll3), Dim()])
    else
      v.push offset               = Int( "offset",             :dir => :in)
      v.push d_kappavstore        = Real("d_kappavstore",      :dir => :in, :dim => [Dim()])
      v.push d_muvstore           = Real("d_muvstore",         :dir => :in, :dim => [Dim()])
    end
    v.push *dudl = ["x", "y", "z"].collect { |a1|
      Real("du#{a1}d#{a1}l", :dir => :in)
    }

    # future check: avoiding these to reduce registers, instead will recompute where needed
    v.push duxdxl_plus_duydyl     = Real("duxdxl_plus_duydyl", :dir => :in)
    v.push duxdxl_plus_duzdzl     = Real("duxdxl_plus_duzdzl", :dir => :in)
    v.push duydyl_plus_duzdzl     = Real("duydyl_plus_duzdzl", :dir => :in)

    v.push duxdyl_plus_duydxl     = Real("duxdyl_plus_duydxl", :dir => :in)
    v.push duzdxl_plus_duxdzl     = Real("duzdxl_plus_duxdzl", :dir => :in)
    v.push duzdyl_plus_duydzl     = Real("duzdyl_plus_duydzl", :dir => :in)
    v.push sigma_xx               = Real("sigma_xx",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yy               = Real("sigma_yy",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_zz               = Real("sigma_zz",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xy               = Real("sigma_xy",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xz               = Real("sigma_xz",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yz               = Real("sigma_yz",           :dir => :out, :dim => [Dim()], :private => true )

    p = Procedure( function_name, v, :local => true ) {
      decl lambdal = Real("lambdal"), mul = Real("mul"), lambdalplus2mul = Real("lambdalplus2mul"), kappal = Real("kappal")
      comment()

      if use_cuda_shared_async then
        print mul === sh_mul[tx,0]
        print kappal === sh_mul[tx,1]
      else
        print mul === d_muvstore[offset]
        print kappal === d_kappavstore[offset]
      end

      # not needed, scaling for attenuation will be done in prepare_elastic_elements.F90
      #print If(attenuation) {
      #  print mul === mul * one_minus_sum_beta_use
      #}
      comment()

      print lambdalplus2mul === kappal + mul * 1.33333333333333333333
      print lambdal === lambdalplus2mul - mul * 2.0

      comment()

      print sigma_xx.dereference === lambdalplus2mul * dudl[0] + lambdal * duydyl_plus_duzdzl
      #print sigma_xx.dereference === lambdalplus2mul * dudl[0] + lambdal * (dudl[1] + dudl[2])  # duydyl_plus_duzdzl

      print sigma_yy.dereference === lambdalplus2mul * dudl[1] + lambdal * duxdxl_plus_duzdzl
      #print sigma_yy.dereference === lambdalplus2mul * dudl[1] + lambdal * (dudl[0] + dudl[2])  # duxdxl_plus_duzdzl

      print sigma_zz.dereference === lambdalplus2mul * dudl[2] + lambdal * duxdxl_plus_duydyl
      #print sigma_zz.dereference === lambdalplus2mul * dudl[2] + lambdal * (dudl[0] + dudl[1])  # duxdxl_plus_duydyl

      print sigma_xy.dereference === mul*duxdyl_plus_duydxl
      print sigma_xz.dereference === mul*duzdxl_plus_duxdzl
      print sigma_yz.dereference === mul*duzdyl_plus_duydzl
    }
    return p
  end

  def BOAST::compute_element_cm_tiso_org(n_gll3 = 125, use_cuda_shared_async = false)
    # not used anymore... left here for reference
    # original routine; computes c11,c12,.. based on muv,kappav,muh,..

    function_name = "compute_element_cm_tiso_org"
    v = []
    ngll3 = Int("NGLL3", :const => n_gll3)

    if use_cuda_shared_async then
      v.push tx                     = Int( "tx",                 :dir => :in)
      v.push sh_mul                 = Real("sh_mul",             :dir => :in, :dim => [Dim(ngll3), Dim()])
    else
      v.push offset                 = Int( "offset", :dir => :in)
      v.push d_kappavstore          = Real("d_kappavstore",          :dir => :in, :dim => [Dim()])
      v.push d_muvstore             = Real("d_muvstore",             :dir => :in, :dim => [Dim()])
      v.push d_kappahstore          = Real("d_kappahstore",          :dir => :in, :dim => [Dim()])
      v.push d_muhstore             = Real("d_muhstore",             :dir => :in, :dim => [Dim()])
      v.push d_eta_anisostore       = Real("d_eta_anisostore",       :dir => :in, :dim => [Dim()])
    end
    dudl = ["x", "y", "z"].collect { |a1|
      ["x", "y", "z"].collect { |a2|
        Real("du#{a1}d#{a2}l", :dir => :in)
      }
    }
    v.push *(dudl.flatten)
    v.push duxdyl_plus_duydxl     = Real("duxdyl_plus_duydxl", :dir => :in)
    v.push duzdxl_plus_duxdzl     = Real("duzdxl_plus_duxdzl", :dir => :in)
    v.push duzdyl_plus_duydzl     = Real("duzdyl_plus_duydzl", :dir => :in)
    v.push iglob                  = Int( "iglob",              :dir => :in)
    if use_cuda_shared_async then
      v.push sh_rstore              = Real("sh_rstore",        :dir => :in,  :dim => [Dim(3), Dim()])
    else
      v.push d_rstore               = Real("d_rstore",         :dir => :in,  :dim => [Dim(3), Dim()])
    end
    v.push sigma_xx               = Real("sigma_xx",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yy               = Real("sigma_yy",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_zz               = Real("sigma_zz",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xy               = Real("sigma_xy",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xz               = Real("sigma_xz",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yz               = Real("sigma_yz",           :dir => :out, :dim => [Dim()], :private => true )

    p = Procedure( function_name, v, :local => true ) {
      decl kappavl = Real("kappavl"), muvl = Real("muvl"), kappahl = Real("kappahl"), muhl = Real("muhl")
      decl rhovpvsq = Real("rhovpvsq"), rhovphsq = Real("rhovphsq"), rhovsvsq = Real("rhovsvsq"), rhovshsq = Real("rhovshsq"), eta_aniso = Real("eta_aniso")
      decl costheta = Real("costheta"), sintheta = Real("sintheta"), cosphi = Real("cosphi"), sinphi = Real("sinphi")
      decl costhetasq = Real("costhetasq"), sinthetasq = Real("sinthetasq"), cosphisq = Real("cosphisq"), sinphisq = Real("sinphisq"), costhetafour = Real("costhetafour"), sinthetafour = Real("sinthetafour"), cosphifour = Real("cosphifour"), sinphifour = Real("sinphifour")
      decl costwotheta = Real("costwotheta"), sintwotheta = Real("sintwotheta"), costwophi = Real("costwophi"), sintwophi = Real("sintwophi"), cosfourtheta = Real("cosfourtheta"), cosfourphi = Real("cosfourphi")
      decl costwothetasq = Real("costwothetasq"), costwophisq = Real("costwophisq"), sintwophisq = Real("sintwophisq")
      decl etaminone = Real("etaminone"), twoetaminone = Real("twoetaminone")
      decl two_eta_aniso = Real("two_eta_aniso"), four_eta_aniso = Real("four_eta_aniso"), six_eta_aniso = Real("six_eta_aniso")
      decl two_rhovsvsq = Real("two_rhovsvsq"), two_rhovshsq = Real("two_rhovshsq")
      decl four_rhovsvsq = Real("four_rhovsvsq"), four_rhovshsq = Real("four_rhovshsq")
      c = (0..5).collect { |indx1|
            (0..5).collect { |indx2|
              if (indx2 < indx1) then
                nil
              else
                Real( "c#{indx1+1}#{indx2+1}" )
              end
            }
          }
      decl *(c.flatten.reject { |e| e.nil? })
      decl theta = Real("theta"), phi = Real("phi")
      comment()

      if use_cuda_shared_async then
        comment("// CUDA asynchronuous memory copies")
        print muvl === sh_mul[tx,0]
        print muhl === sh_mul[tx,1]
        print kappavl === sh_mul[tx,2];
        print kappahl === sh_mul[tx,3];
        print eta_aniso === sh_mul[tx,4];
      else
        print kappavl === d_kappavstore[offset]
        print muvl === d_muvstore[offset]
        print kappahl === d_kappahstore[offset]
        print muhl === d_muhstore[offset]
        print eta_aniso === d_eta_anisostore[offset]
      end
      comment()

      # not needed, scaling for attenuation will be done in prepare_elastic_elements.F90
      #print If(attenuation) {
      #  print muvl === muvl * one_minus_sum_beta_use
      #  print muhl === muhl * one_minus_sum_beta_use
      #}

      print rhovpvsq === kappavl + muvl * 1.33333333333333333333
      print rhovphsq === kappahl + muhl * 1.33333333333333333333

      print rhovsvsq === muvl
      print rhovshsq === muhl
      comment()

      if use_cuda_shared_async then
        print theta === sh_rstore[1,tx]
        print phi === sh_rstore[2,tx]
      else
        print theta === d_rstore[1,iglob]
        print phi === d_rstore[2,iglob]
      end
      comment()

      if (get_lang == CL) then
        print sintheta     === sincos(theta,     costheta.address)
        print sinphi       === sincos(phi,       cosphi.address)
        print sintwotheta  === sincos(theta*2.0, costwotheta.address)
        print sintwophi    === sincos(phi*2.0,   costwophi.address)
        print cosfourtheta ===    cos(theta*4.0)
        print cosfourphi   ===    cos(phi*4.0)
      else
        if (get_default_real_size == 4) then
          print sincosf(theta,     sintheta.address,    costheta.address)
          print sincosf(phi,       sinphi.address,      cosphi.address)
          print sincosf(theta*2.0, sintwotheta.address, costwotheta.address)
          print sincosf(phi*2.0,   sintwophi.address,   costwophi.address)
          print cosfourtheta === cosf(theta*4.0)
          print cosfourphi   === cosf(phi*4.0)
        else
          print costheta     === cos(theta)
          print sintheta     === sin(theta)
          print cosphi       === cos(phi)
          print sinphi       === sin(phi)
          print costwotheta  === cos(theta*2.0)
          print sintwotheta  === sin(theta*2.0)
          print costwophi    === cos(phi  *2.0)
          print sintwophi    === sin(phi  *2.0)
          print cosfourtheta === cos(theta*4.0)
          print cosfourphi   === cos(phi  *4.0)
        end
      end
      comment()

      print costhetasq === costheta * costheta
      print sinthetasq === sintheta * sintheta
      print cosphisq === cosphi * cosphi
      print sinphisq === sinphi * sinphi

      print costhetafour === costhetasq * costhetasq
      print sinthetafour === sinthetasq * sinthetasq
      print cosphifour === cosphisq * cosphisq
      print sinphifour === sinphisq * sinphisq

      print costwothetasq === costwotheta * costwotheta

      print costwophisq === costwophi * costwophi
      print sintwophisq === sintwophi * sintwophi

      print etaminone === eta_aniso - 1.0
      print twoetaminone === eta_aniso * 2.0 - 1.0


      print two_eta_aniso === eta_aniso * 2.0
      print four_eta_aniso === eta_aniso * 4.0
      print six_eta_aniso === eta_aniso * 6.0

      print two_rhovsvsq === rhovsvsq * 2.0
      print two_rhovshsq === rhovshsq * 2.0

      print four_rhovsvsq === rhovsvsq * 4.0
      print four_rhovshsq === rhovshsq * 4.0

      comment()

      print c[0][0] === rhovphsq*sinphifour + cosphisq*sinphisq*
                    (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)*
                    sinthetasq)*2.0 + cosphifour*
                    (rhovphsq*costhetafour + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)*
                    costhetasq*sinthetasq*2.0 + rhovpvsq*sinthetafour)

      print c[0][1] === ((rhovphsq - two_rhovshsq)*(cosfourphi + 3.0)*costhetasq)*0.25 -
                    four_rhovshsq*cosphisq*costhetasq*sinphisq +
                    (rhovphsq*(costwotheta*4.0 + cosfourtheta + 11.0)*sintwophisq)*0.03125 +
                    eta_aniso*(rhovphsq - two_rhovsvsq)*(cosphifour +
                    cosphisq*costhetasq*sinphisq*2.0 + sinphifour)*sinthetasq +
                    rhovpvsq*cosphisq*sinphisq*sinthetafour -
                    rhovsvsq*sintwophisq*sinthetafour

      print c[0][2] === (cosphisq*(rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq -
                    eta_aniso*rhovsvsq*12.0 + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq -
                    four_eta_aniso*rhovsvsq)*cosfourtheta))*0.125 +
                    sinphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq +
                    (rhovphsq - two_rhovshsq)*sinthetasq)

      print c[0][3] === costheta*sinphi*((cosphisq*
                    (-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq +
                    (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq -
                    four_eta_aniso*rhovsvsq)*costwotheta))*0.5 +
                    (etaminone*rhovphsq + (rhovshsq - eta_aniso*rhovsvsq)*2.0)*sinphisq)* sintheta

      print c[0][4] === cosphi*costheta*((cosphisq* (-rhovphsq + rhovpvsq +
                    (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)*
                    costwotheta))*0.5 + etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sintheta

      print c[0][5] === (cosphi*sinphi*(cosphisq* (-rhovphsq + rhovpvsq +
                    (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq -
                    four_eta_aniso*rhovsvsq)*costwotheta) +
                    etaminone*(rhovphsq - two_rhovsvsq)*sinphisq*2.0)*sinthetasq)*0.5

      print c[1][1] === rhovphsq*cosphifour + cosphisq*sinphisq*
                    (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)*
                    sinthetasq)*2.0 + sinphifour*
                    (rhovphsq*costhetafour + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)*
                    costhetasq*sinthetasq*2.0 + rhovpvsq*sinthetafour)

      print c[1][2] === ((rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - eta_aniso*rhovsvsq*12.0 +
                    (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)*
                    cosfourtheta)*sinphisq)*0.125 +
                    cosphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq +
                    (rhovphsq - two_rhovshsq)*sinthetasq)

      print c[1][3] === costheta*sinphi*(etaminone*(rhovphsq - two_rhovsvsq)*cosphisq +
                    ((-rhovphsq + rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq +
                    four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*0.5)*sintheta

      print c[1][4] === cosphi*costheta*((etaminone*rhovphsq + (rhovshsq - eta_aniso*rhovsvsq)*2.0)*
                    cosphisq + ((-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq +
                    (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq -
                    four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*0.5)*sintheta

      print c[1][5] === (cosphi*sinphi*(etaminone*(rhovphsq - two_rhovsvsq)*cosphisq*2.0 +
                    (-rhovphsq + rhovpvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq -
                    four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*sinthetasq)*0.5

      print c[2][2] === rhovpvsq*costhetafour + (eta_aniso*(rhovphsq - two_rhovsvsq) + two_rhovsvsq)*
                    costhetasq*sinthetasq*2.0 + rhovphsq*sinthetafour

      print c[2][3] === -((rhovphsq - rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq -
                    four_eta_aniso*rhovsvsq)*costwotheta)*sinphi*sintwotheta)*0.25

      print c[2][4] === -(cosphi*(rhovphsq - rhovpvsq +
                    (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)*
                    costwotheta)*sintwotheta)*0.25

      print c[2][5] === -((rhovphsq - rhovpvsq - four_rhovshsq + four_rhovsvsq +
                    (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)*
                    costwotheta)*sintwophi*sinthetasq)*0.25

      print c[3][3] === cosphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) +
                    sinphisq*(rhovsvsq*costwothetasq +
                    (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

      print c[3][4] === ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq +
                    four_eta_aniso*rhovsvsq + (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq +
                    etaminone*rhovsvsq*4.0)*costwotheta)*sintwophi*sinthetasq)*0.25

      print c[3][5] === -(cosphi*costheta*((rhovshsq - rhovsvsq)*cosphisq -
                    ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq +
                    four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq +
                    four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*0.5)* sintheta)

      print c[4][4] === sinphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) +
                    cosphisq*(rhovsvsq*costwothetasq +
                    (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

      print c[4][5] === costheta*sinphi*((cosphisq*
                    (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq +
                    four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq +
                    four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta))*0.5 +
                    (-rhovshsq + rhovsvsq)*sinphisq)*sintheta

      print c[5][5] === rhovshsq*costwophisq*costhetasq -
                    (rhovphsq - two_rhovshsq)*cosphisq*costhetasq*sinphisq*2.0 +
                    (rhovphsq*(costwotheta*4.0 + cosfourtheta + 11.0)*sintwophisq)*0.03125 -
                    (rhovsvsq*(cosfourphi*(-2.0) + cos(phi*4.0 - theta*2.0) - costwotheta*2.0 +
                    cos((phi*2.0 + theta)*2.0) - 6.0)*sinthetasq)*0.125 +
                    rhovpvsq*cosphisq*sinphisq*sinthetafour -
                    (eta_aniso*(rhovphsq - two_rhovsvsq)*sintwophisq*sinthetafour)*0.5

      print sigma_xx.dereference === c[0][0]*dudl[0][0] + c[0][5]*duxdyl_plus_duydxl + c[0][1]*dudl[1][1] +
                                     c[0][4]*duzdxl_plus_duxdzl + c[0][3]*duzdyl_plus_duydzl + c[0][2]*dudl[2][2]
      print sigma_yy.dereference === c[0][1]*dudl[0][0] + c[1][5]*duxdyl_plus_duydxl + c[1][1]*dudl[1][1] +
                                     c[1][4]*duzdxl_plus_duxdzl + c[1][3]*duzdyl_plus_duydzl + c[1][2]*dudl[2][2]
      print sigma_zz.dereference === c[0][2]*dudl[0][0] + c[2][5]*duxdyl_plus_duydxl + c[1][2]*dudl[1][1] +
                                     c[2][4]*duzdxl_plus_duxdzl + c[2][3]*duzdyl_plus_duydzl + c[2][2]*dudl[2][2]
      print sigma_xy.dereference === c[0][5]*dudl[0][0] + c[5][5]*duxdyl_plus_duydxl + c[1][5]*dudl[1][1] +
                                     c[4][5]*duzdxl_plus_duxdzl + c[3][5]*duzdyl_plus_duydzl + c[2][5]*dudl[2][2]
      print sigma_xz.dereference === c[0][4]*dudl[0][0] + c[4][5]*duxdyl_plus_duydxl + c[1][4]*dudl[1][1] +
                                     c[4][4]*duzdxl_plus_duxdzl + c[3][4]*duzdyl_plus_duydzl + c[2][4]*dudl[2][2]
      print sigma_yz.dereference === c[0][3]*dudl[0][0] + c[3][5]*duxdyl_plus_duydxl + c[1][3]*dudl[1][1] +
                                     c[3][4]*duzdxl_plus_duxdzl + c[3][3]*duzdyl_plus_duydzl + c[2][3]*dudl[2][2]
    }
    return p
  end


  def BOAST::compute_element_cm_tiso(n_gll3 = 125, use_cuda_shared_async = false)
    # using pre-computed c11,c12,..
    # same computations as fully anisotropic element
    # see: compute_element_cm_aniso() routine of CPU routines

    function_name = "compute_element_cm_tiso"
    v = []
    ngll3 = Int("NGLL3", :const => n_gll3)

    if use_cuda_shared_async then
      v.push tx                     = Int( "tx",                 :dir => :in)
      v.push sh_cstore              = Real("sh_cstore",          :dir => :in, :dim => [Dim(ngll3), Dim()])
    else
      v.push offset                 = Int( "offset", :dir => :in)
      d_cstore = (0..5).collect { |indx1|
                  (0..5).collect { |indx2|
                    if (indx2 < indx1) then
                      nil
                    else
                      Real( "d_c#{indx1+1}#{indx2+1}store", :dir => :in, :dim => [Dim()] )
                    end
                  }
                }
      v.push *(d_cstore.flatten.reject { |e| e.nil? })
    end
    dudl = ["x", "y", "z"].collect { |a1|
      ["x", "y", "z"].collect { |a2|
        Real("du#{a1}d#{a2}l", :dir => :in)
      }
    }
    # elements outside the diagonal are unused
    v.push *(dudl.flatten)
    v.push duxdyl_plus_duydxl = Real("duxdyl_plus_duydxl", :dir => :in)
    v.push duzdxl_plus_duxdzl = Real("duzdxl_plus_duxdzl", :dir => :in)
    v.push duzdyl_plus_duydzl = Real("duzdyl_plus_duydzl", :dir => :in)
    v.push sigma_xx = Real("sigma_xx", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yy = Real("sigma_yy", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_zz = Real("sigma_zz", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xy = Real("sigma_xy", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xz = Real("sigma_xz", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yz = Real("sigma_yz", :dir => :out, :dim => [Dim()], :private => true )

    p = Procedure( function_name, v, :local => true ) {
      c = (0..5).collect { |indx1|
            (0..5).collect { |indx2|
              if (indx2 < indx1) then
                nil
              else
                Real( "c#{indx1+1}#{indx2+1}" )
              end
            }
          }
      decl *(c.flatten.reject { |e| e.nil? })
      comment()

      if use_cuda_shared_async then
        comment("// CUDA asynchronuous memory copies")
        offset = 0
        (0..5).each { |indx1|
          (0..5).each { |indx2|
            print c[indx1][indx2] === sh_cstore[tx,offset] unless indx2 < indx1
            offset += 1 unless indx2 < indx1
          }
        }
      else
        (0..5).each { |indx1|
          (0..5).each { |indx2|
            print c[indx1][indx2] === d_cstore[indx1][indx2][offset] unless indx2 < indx1
          }
        }
      end
      comment()

      # not needed anymore, scaling will be done in prepare_elastic_elements.F90
      #print If(attenuation) {
      #  print minus_sum_beta === one_minus_sum_beta_use - 1.0
      #  print mul === c[3][3] * minus_sum_beta
      #
      #  print c[0][0] === c[0][0] + mul * 1.33333333333333333333
      #  print c[0][1] === c[0][1] - mul * 0.66666666666666666666
      #  print c[0][2] === c[0][2] - mul * 0.66666666666666666666
      #  print c[1][1] === c[1][1] + mul * 1.33333333333333333333
      #  print c[1][2] === c[1][2] - mul * 0.66666666666666666666
      #  print c[2][2] === c[2][2] + mul * 1.33333333333333333333
      #  print c[3][3] === c[3][3] + mul
      #  print c[4][4] === c[4][4] + mul
      #  print c[5][5] === c[5][5] + mul
      #}
      #comment()

      print sigma_xx.dereference === c[0][0]*dudl[0][0] + c[0][5]*duxdyl_plus_duydxl + c[0][1]*dudl[1][1] \
                                   + c[0][4]*duzdxl_plus_duxdzl + c[0][3]*duzdyl_plus_duydzl + c[0][2]*dudl[2][2]
      print sigma_yy.dereference === c[0][1]*dudl[0][0] + c[1][5]*duxdyl_plus_duydxl + c[1][1]*dudl[1][1] \
                                   + c[1][4]*duzdxl_plus_duxdzl + c[1][3]*duzdyl_plus_duydzl + c[1][2]*dudl[2][2]
      print sigma_zz.dereference === c[0][2]*dudl[0][0] + c[2][5]*duxdyl_plus_duydxl + c[1][2]*dudl[1][1] \
                                   + c[2][4]*duzdxl_plus_duxdzl + c[2][3]*duzdyl_plus_duydzl + c[2][2]*dudl[2][2]
      print sigma_xy.dereference === c[0][5]*dudl[0][0] + c[5][5]*duxdyl_plus_duydxl + c[1][5]*dudl[1][1] \
                                   + c[4][5]*duzdxl_plus_duxdzl + c[3][5]*duzdyl_plus_duydzl + c[2][5]*dudl[2][2]
      print sigma_xz.dereference === c[0][4]*dudl[0][0] + c[4][5]*duxdyl_plus_duydxl + c[1][4]*dudl[1][1] \
                                   + c[4][4]*duzdxl_plus_duxdzl + c[3][4]*duzdyl_plus_duydzl + c[2][4]*dudl[2][2]
      print sigma_yz.dereference === c[0][3]*dudl[0][0] + c[3][5]*duxdyl_plus_duydxl + c[1][3]*dudl[1][1] \
                                   + c[3][4]*duzdxl_plus_duxdzl + c[3][3]*duzdyl_plus_duydzl + c[2][3]*dudl[2][2]
    }
    return p
  end


#----------------------------------------------------------------------

## inner core

  def BOAST::compute_element_ic_aniso
    function_name = "compute_element_ic_aniso"
    v = []
    v.push offset = Int( "offset", :dir => :in)
    v.push d_c11store              = Real("d_c11store",              :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_c12store              = Real("d_c12store",              :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_c13store              = Real("d_c13store",              :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_c33store              = Real("d_c33store",              :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_c44store              = Real("d_c44store",              :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push *dudl = ["x", "y", "z"].collect { |a1|
      Real("du#{a1}d#{a1}l", :dir => :in)
    }
    # elements outside the diagonal are unused
    v.push duxdyl_plus_duydxl = Real("duxdyl_plus_duydxl", :dir => :in)
    v.push duzdxl_plus_duxdzl = Real("duzdxl_plus_duxdzl", :dir => :in)
    v.push duzdyl_plus_duydzl = Real("duzdyl_plus_duydzl", :dir => :in)
    v.push sigma_xx = Real("sigma_xx", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yy = Real("sigma_yy", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_zz = Real("sigma_zz", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xy = Real("sigma_xy", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xz = Real("sigma_xz", :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yz = Real("sigma_yz", :dir => :out, :dim => [Dim()], :private => true )

    p = Procedure( function_name, v, :local => true ) {
      decl c11 = Real("c11")
      decl c12 = Real("c12")
      decl c13 = Real("c13")
      decl c33 = Real("c33")
      decl c44 = Real("c44")
      decl c66 = Real("c66")

      print c11 === d_c11store[offset]
      print c12 === d_c12store[offset]
      print c13 === d_c13store[offset]
      print c33 === d_c33store[offset]
      print c44 === d_c44store[offset]
      print c66 === 0.5 * (c11-c12)

      print sigma_xx.dereference === c11*dudl[0] + c12*dudl[1] + c13*dudl[2]
      print sigma_yy.dereference === c12*dudl[0] + c11*dudl[1] + c13*dudl[2]
      print sigma_zz.dereference === c13*dudl[0] + c13*dudl[1] + c33*dudl[2]
      print sigma_xy.dereference === c66*duxdyl_plus_duydxl
      print sigma_xz.dereference === c44*duzdxl_plus_duxdzl
      print sigma_yz.dereference === c44*duzdyl_plus_duydzl
    }
    return p
  end

  def BOAST::compute_element_ic_iso
    function_name = "compute_element_ic_iso"
    v = []
    v.push offset = Int( "offset", :dir => :in)
    v.push d_kappavstore          = Real("d_kappavstore",          :dir => :in, :dim => [Dim()])
    v.push d_muvstore             = Real("d_muvstore",             :dir => :in, :dim => [Dim()])
    v.push *dudl = ["x", "y", "z"].collect { |a1|
      Real("du#{a1}d#{a1}l", :dir => :in)
    }

    # future check: avoiding these to reduce registers, instead will recompute where needed
    v.push duxdxl_plus_duydyl     = Real("duxdxl_plus_duydyl", :dir => :in)
    v.push duxdxl_plus_duzdzl     = Real("duxdxl_plus_duzdzl", :dir => :in)
    v.push duydyl_plus_duzdzl     = Real("duydyl_plus_duzdzl", :dir => :in)

    v.push duxdyl_plus_duydxl     = Real("duxdyl_plus_duydxl", :dir => :in)
    v.push duzdxl_plus_duxdzl     = Real("duzdxl_plus_duxdzl", :dir => :in)
    v.push duzdyl_plus_duydzl     = Real("duzdyl_plus_duydzl", :dir => :in)
    v.push sigma_xx               = Real("sigma_xx",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yy               = Real("sigma_yy",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_zz               = Real("sigma_zz",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xy               = Real("sigma_xy",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xz               = Real("sigma_xz",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yz               = Real("sigma_yz",           :dir => :out, :dim => [Dim()], :private => true )

    p = Procedure( function_name, v, :local => true ) {
      decl lambdal = Real("lambdal"), mul = Real("mul"), lambdalplus2mul = Real("lambdalplus2mul"), kappal = Real("kappal")

      print kappal === d_kappavstore[offset]
      print mul === d_muvstore[offset]

      print lambdalplus2mul === kappal + mul * 1.33333333333333333333
      print lambdal === lambdalplus2mul - mul * 2.0

      print sigma_xx.dereference === lambdalplus2mul * dudl[0] + lambdal * duydyl_plus_duzdzl
      #print sigma_xx.dereference === lambdalplus2mul * dudl[0] + lambdal * (dudl[1] + dudl[2])  # duydyl_plus_duzdzl

      print sigma_yy.dereference === lambdalplus2mul * dudl[1] + lambdal * duxdxl_plus_duzdzl
      #print sigma_yy.dereference === lambdalplus2mul * dudl[1] + lambdal * (dudl[0] + dudl[2])  # duxdxl_plus_duzdzl

      print sigma_zz.dereference === lambdalplus2mul * dudl[2] + lambdal * duxdxl_plus_duydyl
      #print sigma_zz.dereference === lambdalplus2mul * dudl[2] + lambdal * (dudl[0] + dudl[1])  # duxdxl_plus_duydyl

      print sigma_xy.dereference === mul * duxdyl_plus_duydxl
      print sigma_xz.dereference === mul * duzdxl_plus_duxdzl
      print sigma_yz.dereference === mul * duzdyl_plus_duydzl
    }
    return p
  end


#----------------------------------------------------------------------
#
# main kernel
#
#----------------------------------------------------------------------

  def BOAST::inner_core_impl_kernel_forward(ref = true, elem_per_thread = 1, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, coloring_min_nspec_inner_core = 1000, i_flag_in_fictitious_cube = 11)

    # setting default values
    launch_bounds = false
    min_blocks = 7

    # iso/tiso kernel (default)
    kernel = BOAST::impl_kernel(:inner_core, true, ref, elem_per_thread, mesh_coloring, textures_fields, textures_constants, unroll_loops, n_gllx, n_gll2, n_gll3, n_gll3_padded, n_sls, coloring_min_nspec_inner_core, i_flag_in_fictitious_cube, launch_bounds, min_blocks, false)

    # aniso kernel
    kernel_aniso = BOAST::impl_kernel(:inner_core, true, ref, elem_per_thread, mesh_coloring, textures_fields, textures_constants, unroll_loops, n_gllx, n_gll2, n_gll3, n_gll3_padded, n_sls, coloring_min_nspec_inner_core, i_flag_in_fictitious_cube, launch_bounds, min_blocks, true)

    # output both kernels
    if (get_lang == CL) then
      # OpenCL will need for each kernel a full output of subroutines to define a (const char*) variable for each kernel
      kernel_total = [kernel,kernel_aniso]
    else
      # CUDA / HIP
      # to create a single kernel with both procedure's codes
      # (since CUDA async memcpy headers can only appear in a single file)
      str = StringIO::new
      s1 = "" + kernel.to_s + "\n"
      s2 = "" + kernel_aniso.to_s + "\n"
      str << s1 + s2

      kernel_total = CKernel::new(code: str)
      # will need to set a procedure for file outputs
      kernel_total.procedure = [kernel.procedure,kernel_aniso.procedure]
    end

    return kernel_total
  end

  def BOAST::impl_kernel(type, forward, ref = true, elem_per_thread = 1, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = false, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, coloring_min_nspec_inner_core = 1000, i_flag_in_fictitious_cube = 11, launch_bounds = false, min_blocks = 7, only_anisotropy = false)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    # creates function name
    # example: crust_mantle_impl_kernel_forward
    #          crust_mantle_aniso_impl_kernel_forward
    #          crust_mantle_impl_kernel_adjoint
    #          crust_mantle_aniso_impl_kernel_adjoint
    #          inner_core_impl_kernel_forward
    #          ..
    if type == :inner_core then
      function_name = "inner_core_"
    elsif type == :crust_mantle then
      function_name = "crust_mantle_"
    else
      raise "Unsupported_type : #{type}!"
    end
    if only_anisotropy then
      function_name += "aniso_"
    end
    function_name += "impl_kernel_"
    # ending
    if forward then
      function_name += "forward"
    else
      function_name += "adjoint"
    end

    # flag for CUDA_SHARED_ASYNC modifications
    use_cuda_shared_async = false
    if (type == :crust_mantle and get_lang == CUDA and forward) then
      # due to an issue with linking and multiple definitions due to #include <cuda/pipeline>
      # we limit this feature to the crust_mantle_**_forward.cu implementation
      ## un/comment to enable/disable:
      use_cuda_shared_async = true
    end

    v = []
    v.push nb_blocks_to_compute    = Int("nb_blocks_to_compute",     :dir => :in)
    v.push d_ibool                 = Int("d_ibool",                  :dir => :in, :dim => [Dim()] )
    if type == :inner_core then
      v.push d_idoubling           = Int("d_idoubling",              :dir => :in, :dim => [Dim()] )
    elsif type == :crust_mantle then
      v.push d_ispec_is_tiso       = Int("d_ispec_is_tiso",          :dir => :in, :dim => [Dim()] )
    end
    v.push d_phase_ispec_inner     = Int("d_phase_ispec_inner",      :dir => :in, :dim => [Dim()] )
    v.push num_phase_ispec         = Int("num_phase_ispec",          :dir => :in )
    v.push d_iphase                = Int("d_iphase",                 :dir => :in )
    v.push deltat                  = Real("deltat",                  :dir => :in )
    v.push use_mesh_coloring_gpu   = Int("use_mesh_coloring_gpu",    :dir => :in )
    v.push d_displ                 = Real("d_displ",                 :dir => :in, :restrict => true, :dim => [Dim(3), Dim()] )
    v.push d_accel                 = Real("d_accel",                 :dir => :inout, :dim => [Dim(3), Dim()] )
    v.push *d_xi = [d_xix          = Real("d_xix",                   :dir => :in, :restrict => true, :dim => [Dim()] ),
                    d_xiy = Real("d_xiy",:dir => :in, :restrict => true, :dim => [Dim()] ),
                    d_xiz = Real("d_xiz",:dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push *d_eta = [d_etax        = Real("d_etax",                  :dir => :in, :restrict => true, :dim => [Dim()] ),
                     d_etay = Real("d_etay",:dir => :in, :restrict => true, :dim => [Dim()] ),
                     d_etaz = Real("d_etaz",:dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax    = Real("d_gammax",                :dir => :in, :restrict => true, :dim => [Dim()] ),
                       d_gammay = Real("d_gammay",:dir => :in, :restrict => true, :dim => [Dim()] ),
                       d_gammaz = Real("d_gammaz",:dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push d_hprime_xx             = Real("d_hprime_xx",             :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_hprimewgll_xx         = Real("d_hprimewgll_xx",         :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_wgllwgll_xy           = Real("d_wgllwgll_xy",           :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_wgllwgll_xz           = Real("d_wgllwgll_xz",           :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_wgllwgll_yz           = Real("d_wgllwgll_yz",           :dir => :in, :restrict => true, :dim => [Dim()] )

    if only_anisotropy then
      # muv for attenuation case
      v.push d_muvstore              = Real("d_muvstore",              :dir => :in, :restrict => true, :dim => [Dim()] )
    else
      # iso/tiso parameters
      v.push d_kappavstore           = Real("d_kappavstore",           :dir => :in, :restrict => true, :dim => [Dim()] )
      v.push d_muvstore              = Real("d_muvstore",              :dir => :in, :restrict => true, :dim => [Dim()] )
      if type == :crust_mantle then
        v.push d_kappahstore           = Real("d_kappahstore",           :dir => :in, :restrict => true, :dim => [Dim()] )
        v.push d_muhstore              = Real("d_muhstore",              :dir => :in, :restrict => true, :dim => [Dim()] )
        v.push d_eta_anisostore        = Real("d_eta_anisostore",        :dir => :in, :restrict => true, :dim => [Dim()] )
      end
    end

    v.push compute_and_store_strain= Int( "COMPUTE_AND_STORE_STRAIN",:dir => :in)
    v.push epsilondev_xx           = Real("epsilondev_xx",           :dir => :inout, :dim => [Dim()] )
    v.push epsilondev_yy           = Real("epsilondev_yy",           :dir => :inout, :dim => [Dim()] )
    v.push epsilondev_xy           = Real("epsilondev_xy",           :dir => :inout, :dim => [Dim()] )
    v.push epsilondev_xz           = Real("epsilondev_xz",           :dir => :inout, :dim => [Dim()] )
    v.push epsilondev_yz           = Real("epsilondev_yz",           :dir => :inout, :dim => [Dim()] )
    v.push epsilon_trace_over_3    = Real("epsilon_trace_over_3",    :dir => :out, :dim => [Dim()] )

    v.push attenuation             = Int( "ATTENUATION",             :dir => :in)
    v.push partial_phys_dispersion_only = Int( "PARTIAL_PHYS_DISPERSION_ONLY", :dir => :in)
    v.push use_3d_attenuation_arrays    = Int( "USE_3D_ATTENUATION_ARRAYS",    :dir => :in)
    v.push one_minus_sum_beta      = Real("one_minus_sum_beta",      :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push factor_common           = Real("factor_common",           :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push r_xx                    = Real("R_xx",                    :dir => :inout, :dim => [Dim()] )
    v.push r_yy                    = Real("R_yy",                    :dir => :inout, :dim => [Dim()] )
    v.push r_xy                    = Real("R_xy",                    :dir => :inout, :dim => [Dim()] )
    v.push r_xz                    = Real("R_xz",                    :dir => :inout, :dim => [Dim()] )
    v.push r_yz                    = Real("R_yz",                    :dir => :inout, :dim => [Dim()] )
    v.push r_xx_lddrk              = Real("R_xx_lddrk",                    :dir => :inout, :dim => [Dim()] )
    v.push r_yy_lddrk              = Real("R_yy_lddrk",                    :dir => :inout, :dim => [Dim()] )
    v.push r_xy_lddrk              = Real("R_xy_lddrk",                    :dir => :inout, :dim => [Dim()] )
    v.push r_xz_lddrk              = Real("R_xz_lddrk",                    :dir => :inout, :dim => [Dim()] )
    v.push r_yz_lddrk              = Real("R_yz_lddrk",                    :dir => :inout, :dim => [Dim()] )
    v.push alpha_lddrk             = Real("alpha_lddrk",             :dir => :in)
    v.push beta_lddrk              = Real("beta_lddrk",              :dir => :in)
    v.push use_lddrk               = Int( "use_lddrk",               :dir => :in)
    v.push alphaval                = Real("alphaval",                :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push betaval                 = Real("betaval",                 :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push gammaval                = Real("gammaval",                :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push tau_sigmainvval         = Real("tau_sigmainvval",         :dir => :in, :restrict => true, :dim => [Dim()] )

    # anisotropy
    if only_anisotropy then
      # aniso kernel
      #v.push anisotropy              = Int( "ANISOTROPY",              :dir => :in)
      if type == :inner_core then
        v.push d_c11store              = Real("d_c11store",              :dir => :in, :restrict => true, :dim => [Dim()] )
        v.push d_c12store              = Real("d_c12store",              :dir => :in, :restrict => true, :dim => [Dim()] )
        v.push d_c13store              = Real("d_c13store",              :dir => :in, :restrict => true, :dim => [Dim()] )
        v.push d_c33store              = Real("d_c33store",              :dir => :in, :restrict => true, :dim => [Dim()] )
        v.push d_c44store              = Real("d_c44store",              :dir => :in, :restrict => true, :dim => [Dim()] )
      elsif type == :crust_mantle then
        d_cstore = (0..5).collect { |indx1|
          (0..5).collect { |indx2|
            if (indx2 < indx1) then
              nil
            else
              Real( "d_c#{indx1+1}#{indx2+1}store", :dir => :in, :restrict => true, :dim => [Dim()] )
            end
          }
        }
        v.push *(d_cstore.flatten.reject { |e| e.nil? })
      end
    else
      # iso/tiso kernel
      if type == :crust_mantle then
        # new tiso routine needs c11store/c12store/..
        d_cstore = (0..5).collect { |indx1|
          (0..5).collect { |indx2|
            if (indx2 < indx1) then
              nil
            else
              Real( "d_c#{indx1+1}#{indx2+1}store", :dir => :in, :restrict => true, :dim => [Dim()] )
            end
          }
        }
        v.push *(d_cstore.flatten.reject { |e| e.nil? })
        # rstore for iso/tiso kernel
        v.push d_rstore                = Real("d_rstore",                :dir => :in, :restrict => true, :dim => [Dim(3), Dim()] )
      end
    end

    v.push gravity                 = Int( "GRAVITY",                 :dir => :in)
    v.push d_gravity_pre_store     = Real("d_gravity_pre_store",     :dir => :in, :restrict => true, :dim => [Dim(3), Dim()] )
    v.push d_gravity_H             = Real("d_gravity_H",             :dir => :in, :restrict => true, :dim => [Dim(6), Dim()] )
    v.push wgll_cube               = Real("wgll_cube",               :dir => :in, :restrict => true, :dim => [Dim()] )

    if type == :inner_core then
      v.push nspec_strain_only = Int( "NSPEC_INNER_CORE_STRAIN_ONLY",   :dir => :in)
      v.push nspec_inner_core  = Int( "NSPEC_INNER_CORE",               :dir => :in)
    elsif type == :crust_mantle then
      v.push nspec_strain_only = Int( "NSPEC_CRUST_MANTLE_STRAIN_ONLY", :dir => :in)
    end

    # definitions
    ngllx        = Int("NGLLX", :const => n_gllx)
    ngll2        = Int("NGLL2", :const => n_gll2)
    ngll3        = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)
    nsls         = Int("N_SLS", :const => n_sls)
    if type == :inner_core then
      iflag_in_fictitious_cube = Int("IFLAG_IN_FICTITIOUS_CUBE", :const => i_flag_in_fictitious_cube)
    end

    use_mesh_coloring       = Int("USE_MESH_COLORING_GPU",   :const => mesh_coloring)
    use_textures_constants  = Int("USE_TEXTURES_CONSTANTS",  :const => textures_constants)
    use_textures_fields     = Int("USE_TEXTURES_FIELDS",     :const => textures_fields)
    manually_unrolled_loops = Int("MANUALLY_UNROLLED_LOOPS", :const => unroll_loops)
    use_launch_bounds       = Int("USE_LAUNCH_BOUNDS",       :const => launch_bounds)
    launch_min_blocks       = Int("LAUNCH_MIN_BLOCKS",       :const => min_blocks)
    cuda_shared_async       = Int("CUDA_SHARED_ASYNC",       :const => use_cuda_shared_async)

    constants = []
    textures_fields = []
    textures_constants = []

    if type == :inner_core then
      d_displ_tex = Real("d_#{forward ? "":"b_"}displ_ic_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      d_accel_tex = Real("d_#{forward ? "":"b_"}accel_ic_tex", :texture => true, :dir => :in, :dim => [Dim()] )
    elsif type == :crust_mantle then
      d_displ_tex = Real("d_#{forward ? "":"b_"}displ_cm_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      d_accel_tex = Real("d_#{forward ? "":"b_"}accel_cm_tex", :texture => true, :dir => :in, :dim => [Dim()] )
    end
    if get_lang == CL then
      textures_fields.push(d_displ_tex, d_accel_tex)
    end
    if type == :inner_core then
      d_hprime_xx_tex = Real("d_hprime_xx_ic_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      d_hprimewgll_xx_tex = Real("d_hprimewgll_xx_ic_tex", :texture => true, :dir => :in, :dim => [Dim()] )
    elsif type == :crust_mantle then
      d_hprime_xx_tex = Real("d_hprime_xx_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      d_hprimewgll_xx_tex = Real("d_hprimewgll_xx_tex", :texture => true, :dir => :in, :dim => [Dim()] )
    end
    if get_lang == CL then
      textures_constants.push(d_hprime_xx_tex, d_hprimewgll_xx_tex)
    end

    if (get_lang == CUDA or get_lang == HIP) then
      qualifiers = "\n#ifdef #{use_launch_bounds}\n__launch_bounds__(#{ngll3_padded}, #{launch_min_blocks})\n#endif\n"
    elsif(get_lang == CL) then
      qualifiers = "" # "__attribute__((reqd_work_group_size(#{ngll3_padded},1,1))) " # (inefficient)
    end

    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu".gsub("_forward","").gsub("_adjoint",""))
    elsif(get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      # outputs header/define/.. statements into default kernel only,
      # so we can append the aniso kernel afterwards into the same file by calling the same impl_kernel() routine twice
      if (get_lang == CL or ((get_lang == CUDA or get_lang == HIP) and not only_anisotropy)) then
        # header
        make_specfem3d_header(:ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded, :n_sls => n_sls, :coloring_min_nspec_inner_core => coloring_min_nspec_inner_core, :iflag_in_fictitious_cube => i_flag_in_fictitious_cube)

        # texture definitions
        if type == :inner_core then
          #DEACTIVATE USE TEXTURES CONSTANTS
          get_output.puts "#ifdef #{use_textures_constants}"
          get_output.puts "#undef #{use_textures_constants}"
          get_output.puts "#endif"
        end
        #if get_lang == CUDA then
        #  get_output.puts "#ifdef #{use_textures_fields}"
        #    decl d_displ_tex
        #    decl d_accel_tex
        #  get_output.puts "#endif"
        #  get_output.puts "#ifdef #{use_textures_constants}"
        #    decl d_hprime_xx_tex
        #    decl d_hprimewgll_xx_tex
        #  get_output.puts "#endif"
        #end
        comment()

        # compilation pragma info
        if (type == :crust_mantle and (get_lang == CUDA or get_lang == HIP) and forward) then
          get_output.puts "#ifdef #{manually_unrolled_loops}"
          get_output.puts "#pragma message (\"\\n\\nCompiling with: #{manually_unrolled_loops} enabled\\n\")"
          get_output.puts "#endif  // #{manually_unrolled_loops}"
          comment()
        end

        if use_cuda_shared_async then
          # pragma messages
          comment("// CUDA asynchronuous memory copies")
          get_output.puts "#ifdef #{cuda_shared_async}"
          get_output.puts "#pragma message (\"\\n\\nCompiling with: #{cuda_shared_async} enabled\\n\")"
          get_output.puts "#endif  // #{cuda_shared_async}"
          comment()
          # float3 modifications
          #get_output.puts "#ifdef #{cuda_shared_async}"
          #comment("// async memory copies use float3 arrays")
          #get_output.puts "#define USE_FLOAT3"
          #get_output.puts "#endif"
          #get_output.puts "#ifdef USE_FLOAT3"
          #get_output.puts "#pragma message (\"\\n\\nCompiling with: USE_FLOAT3 enabled\\n\")"
          #get_output.puts "#endif"
          #comment()
          # header include files
          comment("// CUDA asynchronuous memory copies")
          get_output.puts "#ifdef #{cuda_shared_async}"
          get_output.puts "#include <cooperative_groups.h>"
          get_output.puts "#include <cuda/barrier>"
          get_output.puts "#include <cuda/pipeline>"
          get_output.puts "#endif  // #{cuda_shared_async}"
          comment()
          # inline routines
          sub_compute_offset = compute_offset(n_gll3, n_sls)
          sub_compute_offset_sh = compute_offset_sh(n_gll3)
          get_output.puts "#ifdef #{cuda_shared_async}"
          print sub_compute_offset
          comment()
          print sub_compute_offset_sh
          get_output.puts "#endif  // #{cuda_shared_async}"
          comment()
        end
      end

      # subroutines
      if type == :inner_core then
        sub_compute_element_att_stress =  compute_element_ic_att_stress(n_gll3, n_sls)
        sub_compute_element_att_memory =  compute_element_ic_att_memory(n_gll3, n_gll3_padded, n_sls)
        sub_compute_element_att_memory_lddrk =  compute_element_ic_att_memory_lddrk(n_gll3, n_gll3_padded, n_sls)
        sub_compute_element_gravity =  compute_element_ic_gravity(n_gll3)
      elsif type == :crust_mantle then
        sub_compute_element_att_stress =  compute_element_cm_att_stress(n_gll3, n_sls, use_cuda_shared_async)
        sub_compute_element_att_memory =  compute_element_cm_att_memory(n_gll3, n_gll3_padded, n_sls)
        sub_compute_element_att_memory_lddrk =  compute_element_cm_att_memory_lddrk(n_gll3, n_gll3_padded, n_sls)
        sub_compute_element_gravity =  compute_element_cm_gravity(n_gll3)
      end
      if use_cuda_shared_async then
        sub_compute_element_att_memory_async = compute_element_cm_att_memory(n_gll3, n_gll3_padded, n_sls, true)
        sub_compute_element_att_memory_lddrk_async = compute_element_cm_att_memory_lddrk(n_gll3, n_gll3_padded, n_sls, true)
        sub_compute_element_gravity_async = compute_element_cm_gravity(n_gll3, true)
      end
      # element aniso/iso/tiso routines
      if type == :inner_core then
        sub_compute_element_ic_aniso = compute_element_ic_aniso
        sub_compute_element_ic_iso = compute_element_ic_iso

      elsif type == :crust_mantle then
        sub_compute_element_cm_aniso = compute_element_cm_aniso
        # old way: computes c11,c12,..
        #sub_compute_element_cm_tiso_org = compute_element_cm_tiso_org
        # new way: uses pre-computed c11,c12,.. - same as in CPU version
        sub_compute_element_cm_tiso = compute_element_cm_tiso
        sub_compute_element_cm_iso = compute_element_cm_iso
        if use_cuda_shared_async then
          # aniso
          sub_compute_element_cm_aniso_async = compute_element_cm_aniso(n_gll3, true)
          # tiso
          sub_compute_element_cm_tiso_async = compute_element_cm_tiso(n_gll3, true)
          # iso
          sub_compute_element_cm_iso_async = compute_element_cm_iso(n_gll3, true)
        end
      end

      # subroutine output
      # - for OpenCL, always output full kernel output, including subroutines.
      # - for CUDA/HIP, outputs subroutines into default kernel only,
      #   so we can append the aniso kernel afterwards into the same file by calling the same impl_kernel() routine twice.
      if (get_lang == CL or ((get_lang == CUDA or get_lang == HIP) and not only_anisotropy)) then
        # function definitions
        # element attenuation routines
        comment()
        print sub_compute_element_att_stress
        comment()
        if use_cuda_shared_async then
          # BOAST doesn't allow to modify the function argument list with #ifdef-statements
          # we thus call the routine twice and put the #ifdef-statement around the whole function
          get_output.puts "#ifdef #{cuda_shared_async}"
          # routine w/ cuda shared async copies
          comment("// CUDA asynchronuous memory copies")
          print sub_compute_element_att_memory_async
          get_output.puts "#else"
          # routine w/out
          comment("//default (w/out asynchronuous memory copies) routine")
          print sub_compute_element_att_memory
          get_output.puts "#endif  // #{cuda_shared_async}"
        else
          print sub_compute_element_att_memory
        end
        comment()
        if use_cuda_shared_async then
          # BOAST doesn't allow to modify the function argument list with #ifdef-statements
          # we thus call the routine twice and put the #ifdef-statement around the whole function
          get_output.puts "#ifdef #{cuda_shared_async}"
          # routine w/ cuda shared async copies
          comment("// CUDA asynchronuous memory copies")
          print sub_compute_element_att_memory_lddrk_async
          get_output.puts "#else"
          # routine w/out
          comment("//default (w/out asynchronuous memory copies) routine")
          print sub_compute_element_att_memory_lddrk
          get_output.puts "#endif  // #{cuda_shared_async}"
        else
          print sub_compute_element_att_memory_lddrk
        end
        comment()
        # element gravity routine
        if use_cuda_shared_async then
          # BOAST doesn't allow to modify the function argument list with #ifdef-statements
          # we thus call the routine twice and put the #ifdef-statement around the whole function
          get_output.puts "#ifdef #{cuda_shared_async}"
          # routine w/ cuda shared async copies
          comment("// CUDA asynchronuous memory copies")
          print sub_compute_element_gravity_async
          get_output.puts "#else"
          # routine w/out
          comment("//default (w/out asynchronuous memory copies) routine")
          print sub_compute_element_gravity
          get_output.puts "#endif  // #{cuda_shared_async}"
        else
          print sub_compute_element_gravity
        end
        comment()
        # element aniso/iso/tiso routines
        if type == :inner_core then
          # aniso
          print sub_compute_element_ic_aniso
          comment()
          # iso
          print sub_compute_element_ic_iso
          comment()
        elsif type == :crust_mantle then
          # aniso
          if use_cuda_shared_async then
            # BOAST doesn't allow to modify the function argument list with #ifdef-statements
            # we thus call the routine twice and put the #ifdef-statement around the whole function
            get_output.puts "#ifdef #{cuda_shared_async}"
            # routine w/ cuda shared async copies
            comment("// CUDA asynchronuous memory copies")
            print sub_compute_element_cm_aniso_async
            get_output.puts "#else"
            # routine w/out
            comment("//default (w/out asynchronuous memory copies) routine")
            print sub_compute_element_cm_aniso
            get_output.puts "#endif  // #{cuda_shared_async}"
          else
            print sub_compute_element_cm_aniso
          end
          comment()
          # iso
          if use_cuda_shared_async then
            # BOAST doesn't allow to modify the function argument list with #ifdef-statements
            # we thus call the routine twice and put the #ifdef-statement around the whole function
            get_output.puts "#ifdef #{cuda_shared_async}"
            # routine w/ cuda shared async copies
            comment("// CUDA asynchronuous memory copies")
            print sub_compute_element_cm_iso_async
            get_output.puts "#else"
            # routine w/out
            comment("//default (w/out asynchronuous memory copies) routine")
            print sub_compute_element_cm_iso
            get_output.puts "#endif  // #{cuda_shared_async}"
          else
            print sub_compute_element_cm_iso
          end
          comment()
          # tiso
          if use_cuda_shared_async then
            # BOAST doesn't allow to modify the function argument list with #ifdef-statements
            # we thus call the routine twice and put the #ifdef-statement around the whole function
            get_output.puts "#ifdef #{cuda_shared_async}"
            # routine w/ cuda shared async copies
            comment("// CUDA asynchronuous memory copies")
            print sub_compute_element_cm_tiso_async
            get_output.puts "#else"
            # routine w/out
            comment("//default (w/out asynchronuous memory copies) routine")
            print sub_compute_element_cm_tiso
            get_output.puts "#endif  // #{cuda_shared_async}"
          else
            # old routine
            #print sub_compute_element_cm_tiso_org
            # new routine
            print sub_compute_element_cm_tiso
          end
          comment()
        end

        comment()
        comment("/*----------------------------------------------*/")
        comment("// main function")
        comment("/*----------------------------------------------*/")
        comment()
      end  # only_anisotropy

      # function name
      if get_lang == CL then
        get_output.puts "#ifdef #{use_textures_fields}"
          get_output.puts "#ifdef #{use_textures_constants}"
            p = Procedure(function_name, v+textures_fields+textures_constants, :constants => constants, :qualifiers => qualifiers)
            open p
            set_indent_level(0)
          get_output.puts "#else"
            p = Procedure(function_name, v+textures_fields, :constants => constants, :qualifiers => qualifiers)
            open p
            set_indent_level(0)
          get_output.puts "#endif"
        get_output.puts "#else"
          get_output.puts "#ifdef #{use_textures_constants}"
            p = Procedure(function_name, v+textures_constants, :constants => constants, :qualifiers => qualifiers)
            open p
            set_indent_level(0)
          get_output.puts "#else"
            p = Procedure(function_name, v, :constants => constants, :qualifiers => qualifiers)
            open p
          get_output.puts "#endif"
        get_output.puts "#endif"

        get_output.puts "#ifdef #{use_textures_fields}"
          decl d_displ_tex.sampler
          decl d_accel_tex.sampler
        get_output.puts "#endif"
        get_output.puts "#ifdef #{use_textures_constants}"
          decl d_hprime_xx_tex.sampler
          decl d_hprimewgll_xx_tex.sampler
        get_output.puts "#endif"
      else
        p = Procedure(function_name, v, :constants => constants, :qualifiers => qualifiers)
        open p
      end

      # declarations
      decl k  = Int("K"), j = Int("J"), i = Int("I")
      # as unsigned short
      #decl k  = Int("K", :size => 2, :signed => false), j = Int("J", :size => 2, :signed => false), i = Int("I", :size => 2, :signed => false)

      # putting declaration of l into for loop: for(int l = 0; ..
      l = Int("l")

      active = (1..elem_per_thread).collect { |e_i| Int("active_#{e_i}", :size => 2, :signed => false) }
      decl *active
      decl offset = Int("offset")
      iglob = (1..elem_per_thread).collect { |e_i| Int("iglob_#{e_i}") }
      decl *iglob
      decl working_element = Int("working_element")
      tempanl = ["x", "y", "z"].collect { |a|
        [ 1, 2, 3 ].collect { |n|
          Real("temp#{a}#{n}l")
        }
      }

      # putting declarations into if (active..){ } scope
      #decl *(tempanl.flatten) # will put these local variable declarations into if-scope {}
      #decl *xil   = [ Real("xixl"),   Real("xiyl"),   Real("xizl")   ]
      #decl *etal  = [ Real("etaxl"),  Real("etayl"),  Real("etazl")  ]
      #decl *gammal= [ Real("gammaxl"),Real("gammayl"),Real("gammazl")]
      #decl jacobianl= Real("jacobianl")

      dudl = ["x", "y", "z"].collect { |a1|
        ["x", "y", "z"].collect { |a2|
          Real("du#{a1}d#{a2}l")
        }
      }

      # putting declarations into if (active..){ } scope
      #decl *(dudl.flatten)
      #decl duxdxl_plus_duydyl = Real("duxdxl_plus_duydyl")
      #decl duxdxl_plus_duzdzl = Real("duxdxl_plus_duzdzl")
      #decl duydyl_plus_duzdzl = Real("duydyl_plus_duzdzl")
      #decl duxdyl_plus_duydxl = Real("duxdyl_plus_duydxl")
      #decl duzdxl_plus_duxdzl = Real("duzdxl_plus_duxdzl")
      #decl duzdyl_plus_duydzl = Real("duzdyl_plus_duydzl")
      #decl templ = Real("templ")

      decl *fac = [1,2,3].collect {|n|
        Real("fac#{n}")
      }

      sigma = ["x", "y", "z"].collect { |a1|
        ["x", "y", "z"].collect { |a2|
          Real("sigma_#{a1}#{a2}")
        }
      }
      decl *(sigma.flatten)
      decl *epsilondev_xx_loc = (1..elem_per_thread).collect{ |e_i| Real("epsilondev_xx_loc_#{e_i}") }
      decl *epsilondev_yy_loc = (1..elem_per_thread).collect{ |e_i| Real("epsilondev_yy_loc_#{e_i}") }
      decl *epsilondev_xy_loc = (1..elem_per_thread).collect{ |e_i| Real("epsilondev_xy_loc_#{e_i}") }
      decl *epsilondev_xz_loc = (1..elem_per_thread).collect{ |e_i| Real("epsilondev_xz_loc_#{e_i}") }
      decl *epsilondev_yz_loc = (1..elem_per_thread).collect{ |e_i| Real("epsilondev_yz_loc_#{e_i}") }

      decl *sum_terms = [1,2,3].collect {|n|
        Real("sum_terms#{n}")
      }
      rho_s_H = (1..elem_per_thread).collect{ |e_i|
        [1,2,3].collect {|n|
          Real("rho_s_H_#{e_i}_#{n}")
        }
      }
      decl *(rho_s_H.flatten)
      comment()

      # shared arrays
      comment("  // shared arrays")
      if use_cuda_shared_async then
        comment("// CUDA asynchronuous memory copies")
        get_output.puts "#ifdef #{cuda_shared_async}"
        decl s_dummy_loc_single = Real("s_dummy_loc", :local => true, :dim => [Dim(3),Dim(ngll3)] )
        get_output.puts "#else"
      end
      decl *s_dummy_loc = ["x", "y", "z"].collect { |a|
        Real("s_dummy#{a}_loc", :local => true, :dim => [Dim(ngll3)] )
      }
      # single array
      #decl s_dummy_loc = Real("s_dummy_loc", :local => true, :dim => [Dim(3),Dim(ngll3)] )
      # or with float3 modification
      #decl s_dummy_loc = Real("s_dummy_loc", :local => true, :dim => [Dim(ngll3)], :size => 4, :vector_length => 3)
      if use_cuda_shared_async then
        get_output.puts "#endif  // #{cuda_shared_async}"
      end
      comment()

      s_temp = ["x", "y", "z"].collect { |a|
        [ 1, 2, 3 ].collect { |n|
          Real("s_temp#{a}#{n}", :local => true, :dim => [Dim(ngll3)] )
        }
      }
      decl *(s_temp.flatten)

      decl sh_hprime_xx     = Real("sh_hprime_xx",     :local => true, :dim => [Dim(ngll2)] )
      decl sh_hprimewgll_xx = Real("sh_hprimewgll_xx", :local => true, :dim => [Dim(ngll2)] )
      comment()

      # possible catch for future use to avoid global loads from d_ispec_is_tiso array
      #decl sh_ispec_is_tiso = Real("sh_ispec_is_tiso", :local => true, :scalar => true )
      #comment()

      # asynchronuous memory copies modification
      if use_cuda_shared_async then
        comment("// CUDA asynchronuous memory copies")
        get_output.puts "#ifdef #{cuda_shared_async}"
        comment("  // attenuation")
        decl sh_factor_common     = Real("sh_factor_common",     :local => true, :dim => [Dim(nsls*ngll3)] )

        decl sh_r_xx     = Real("sh_R_xx",     :local => true, :dim => [Dim(nsls*ngll3)] )
        decl sh_r_yy     = Real("sh_R_yy",     :local => true, :dim => [Dim(nsls*ngll3)] )
        decl sh_r_xy     = Real("sh_R_xy",     :local => true, :dim => [Dim(nsls*ngll3)] )
        decl sh_r_xz     = Real("sh_R_xz",     :local => true, :dim => [Dim(nsls*ngll3)] )
        decl sh_r_yz     = Real("sh_R_yz",     :local => true, :dim => [Dim(nsls*ngll3)] )
        comment()

        decl sh_epsilondev_xx     = Real("sh_epsilondev_xx",     :local => true, :dim => [Dim(ngll3)] )
        decl sh_epsilondev_yy     = Real("sh_epsilondev_yy",     :local => true, :dim => [Dim(ngll3)] )
        decl sh_epsilondev_xy     = Real("sh_epsilondev_xy",     :local => true, :dim => [Dim(ngll3)] )
        decl sh_epsilondev_xz     = Real("sh_epsilondev_xz",     :local => true, :dim => [Dim(ngll3)] )
        decl sh_epsilondev_yz     = Real("sh_epsilondev_yz",     :local => true, :dim => [Dim(ngll3)] )
        comment()

        decl sh_alphaval     = Real("sh_alphaval",     :local => true, :dim => [Dim(nsls)] )
        decl sh_betaval      = Real("sh_betaval",      :local => true, :dim => [Dim(nsls)] )
        decl sh_gammaval     = Real("sh_gammaval",     :local => true, :dim => [Dim(nsls)] )
        comment()

        if only_anisotropy then
          comment("  // shmem muv for attenuation")
          decl sh_mul               = Real("sh_mul",               :local => true, :dim => [Dim(ngll3)] )
          comment("  // shmem order c11store,c12store,.. for aniso")
          decl sh_cstore            = Real("sh_cstore",            :local => true, :dim => [Dim(21*ngll3)] )
        else
          # original tiso routine
          #comment("  // shmem order sh_muv, sh_muh, sh_kappav, sh_kappah, sh_eta_anisostore")
          #decl sh_mul               = Real("sh_mul",               :local => true, :dim => [Dim(5*ngll3)] )
          #comment("  // shmem order rstore_x,rstore_y,rstore_z for original tiso routine")
          #decl sh_rstore            = Real("sh_rstore",            :local => true, :dim => [Dim(3*ngll3)] )
          comment("  // shmem order sh_muv, sh_kappav for iso and attenuation")
          decl sh_mul               = Real("sh_mul",               :local => true, :dim => [Dim(2*ngll3)] )
          comment("  // shmem order c11store,c12store,.. for tiso")
          decl sh_cstore            = Real("sh_cstore",            :local => true, :dim => [Dim(21*ngll3)] )
        end
        comment()

        comment("  // gravity")
        comment("  // d_gravity_pre_store")
        decl sh_gravity_pre_store = Real("sh_gravity_pre_store", :local => true, :dim => [Dim(3*ngll3)] )
        comment("  // d_gravity_H")
        decl sh_gravity_H         = Real("sh_gravity_H",         :local => true, :dim => [Dim(6*ngll3)] )
        decl sh_wgll_cube         = Real("sh_wgll_cube",         :local => true, :dim => [Dim(ngll3)] )
        comment()

        # see: https://nvidia.github.io/libcudacxx/extended_api/synchronization_primitives/barrier.html
        # we'll use
        #   __shared__ block_barrier barrier[2];
        # instead of
        #   __shared__ block_barrier *barrier;
        #   barrier = (block_barrier*) new block_barrier[2];
        # as the former will need fewer registers.
        # nvcc will output a warning about dynamic initialization not supported for __shared__ variables.
        # the pragma here will silence the warning.
        get_output.puts "#pragma diag_suppress static_var_with_dynamic_init"
        get_output.puts "  using block_barrier = cuda::barrier<cuda::thread_scope_block>;"
        get_output.puts "  __shared__ block_barrier barrier[2];"
        get_output.puts "#endif  // #{cuda_shared_async}"
        comment()

        # float3 example
        #decl a = BOAST::CustomType("a", :type_name => "float3", :scalar => true )
        # or
        #decl a = BOAST::Real("a", :size => 4, :vector_length => 3)
      end

      # block id
      decl bx = Int("bx", :const => get_group_id(1) * get_num_groups(0) + get_group_id(0))

      # checks if anything to do
      print If(bx >= nb_blocks_to_compute) { print Return(nil) }
      comment()

      # thread id
      if (elem_per_thread == 1) then
        decl tx = Int("tx", :const => get_local_id(0))
        # as unsigned short
        #decl tx = Int("tx", :size => 2, :signed => false, :const => get_local_id(0))
      else
        # tx needs different thread assignment based on element index
        decl tx = Int("tx")
        # as unsigned short
        #decl tx = Int("tx", :size => 2, :signed => false)
      end
      comment()

      # asynchronuous memory copies modification
      if use_cuda_shared_async then
        # todo: for now requires elem_per_thread == 1
        if elem_per_thread > 1 then
          puts "Error: cuda_shared_async only implemented for elem_per_thread == 1 yet in kernel inner_core_impl_kernel_forward.rb"
          abort "CUDA_SHARED_ASYNC requires elem_per_thread == 1"
        end

        comment("// CUDA asynchronuous memory copies")

        # using plain code for CUDA specific modifications
        # there's no need for BOAST translations to other languages
        endl = "\n"
        cuda_text = "#ifdef #{cuda_shared_async}" + endl
        cuda_text += "  auto block = cooperative_groups::this_thread_block();" + endl
        cuda_text += "  // initializes barriers" + endl
        cuda_text += "  if (tx == 0){" + endl
        cuda_text += "    init(&barrier[0], block.size());" + endl
        cuda_text += "    init(&barrier[1], block.size());" + endl
        cuda_text += "  }" + endl
        cuda_text += "  // creates pipeline" + endl
        cuda_text += "  cuda::pipeline<cuda::thread_scope_thread> pipe = cuda::make_pipeline();" + endl
        cuda_text += "  auto thread = cooperative_groups::this_thread();" + endl
        cuda_text += "  // need this for init the barrier" + endl
        cuda_text += "  block.sync();" + endl
        cuda_text += "  // stage to get displ/hprime/.." + endl
        cuda_text += "  pipe.producer_acquire();" + endl
        cuda_text += "#endif  // #{cuda_shared_async}" + endl + endl
        cuda_text += "" + endl
        get_output.puts cuda_text
      end

      elem_per_thread.times { |elem_index|
        if (elem_per_thread > 1) then
          # for elem_per_thread == 2: a single element needs two passes
          #   kernel blocksize becomes 128 / 2 = 64 as the kernel would be too big otherwise, e.g., for ARM GPUs
          #   first pass -> elem_index == 0 and
          #     tx = threadIdx.x + 128 * 0 / 2 = threadIdx.x
          #   second pass -> elem_index == 1 and
          #     tx = threadIdx.x + 128 * 1 / 2 = threadIdx.x + 64
          print tx === get_local_id(0) + ngll3_padded * elem_index / elem_per_thread
        end
        print active[elem_index] === Ternary( tx < ngll3, 1, 0)
        comment()

        print If(active[elem_index]) {
          if elem_index == 0 then
            get_output.puts "#ifdef #{use_mesh_coloring}"
              print working_element === bx
            get_output.puts "#else"
              print If(use_mesh_coloring_gpu => lambda {
                print working_element === bx
              }, :else => lambda {
                print working_element === d_phase_ispec_inner[bx + num_phase_ispec*(d_iphase-1)]-1
              })
            get_output.puts "#endif"
          end
          __texture_fetch = lambda {
            print iglob[elem_index] === d_ibool[working_element*ngll3 + tx]-1

            get_output.puts "#ifdef #{use_textures_fields}"
              if use_cuda_shared_async then
                get_output.puts "#ifdef #{cuda_shared_async}"
                # single array
                (0..2).each { |indx|
                  print s_dummy_loc_single[indx,tx] === d_displ_tex[iglob[elem_index]*3+indx]
                }
                get_output.puts "#else"
              end
              (0..2).each { |indx|
                print s_dummy_loc[indx][tx] === d_displ_tex[iglob[elem_index]*3+indx]
              }
              if use_cuda_shared_async then
                get_output.puts "#endif  // #{cuda_shared_async}"
              end

            get_output.puts "#else"

              if use_cuda_shared_async then
                get_output.puts "#ifdef #{cuda_shared_async}"
                # single array
                get_output.puts "    cuda::memcpy_async(thread, ((float3 *)s_dummy_loc) + tx, ((float3 *)(d_displ)) + #{iglob[elem_index]}, sizeof(float3), pipe);" + endl
                get_output.puts "#else"
              end
              (0..2).each { |indx|
                print s_dummy_loc[indx][tx] === d_displ[indx, iglob[elem_index]]
                # single array
                #print s_dummy_loc[indx,tx] === d_displ[indx, iglob[elem_index]]
              }
              if use_cuda_shared_async then
                get_output.puts "#endif  // #{cuda_shared_async}"
              end

            get_output.puts "#endif"
          }

          if type == :inner_core then
            print If(d_idoubling[working_element] == iflag_in_fictitious_cube => lambda {
              print active[elem_index] === 0
            }, :else => __texture_fetch )
          elsif type == :crust_mantle then
            __texture_fetch.call
          end
        }
        comment()

        # loading ispec_is_tiso into shared memory, still needs some barriers...
        #print If(tx == 0) {
        #  if use_cuda_shared_async then
        #    get_output.puts "#ifdef #{cuda_shared_async}"
        #    get_output.puts "    cuda::memcpy_async(thread, sh_ispec_is_tiso, d_ispec_is_tiso + working_element, sizeof(int), pipe);"
        #    get_output.puts "#else"
        #  end
        #  print sh_ispec_is_tiso === d_ispec_is_tiso[working_element]
        #  if use_cuda_shared_async then
        #    get_output.puts "#endif  // #{cuda_shared_async}"
        #  end
        #}

        #inner core and crust mantle differ here, but crust mantle implementation though more recent seems odd...
        print If(tx < ngll2) {
          get_output.puts "#ifdef #{use_textures_constants}"
            print sh_hprime_xx[tx] === d_hprime_xx_tex[tx]
            print sh_hprimewgll_xx[tx] === d_hprimewgll_xx_tex[tx]
          get_output.puts "#else"
            if use_cuda_shared_async then
              get_output.puts "#ifdef #{cuda_shared_async}"
              get_output.puts "    cuda::memcpy_async(thread, sh_hprime_xx + tx, d_hprime_xx + tx, sizeof(float), pipe);"
              get_output.puts "    cuda::memcpy_async(thread, sh_hprimewgll_xx + tx, d_hprimewgll_xx + tx, sizeof(float), pipe);"
              get_output.puts "#else"
            end
            print sh_hprime_xx[tx] === d_hprime_xx[tx]
            print sh_hprimewgll_xx[tx] === d_hprimewgll_xx[tx]
            if use_cuda_shared_async then
              get_output.puts "#endif  // #{cuda_shared_async}"
            end
          get_output.puts "#endif"
        }
      }
      comment()

      # barrier
      if use_cuda_shared_async then
        get_output.puts "#ifdef #{cuda_shared_async}"
        get_output.puts "#if defined(USE_TEXTURES_FIELDS) || defined(USE_TEXTURES_CONSTANTS)" + endl
        print barrier(:local)
        get_output.puts "#endif" + endl
        get_output.puts "#else"
      end
      print barrier(:local)
      if use_cuda_shared_async then
        get_output.puts "#endif  // #{cuda_shared_async}"
      end
      comment()

      # asynchronuous memory copies modification
      if use_cuda_shared_async then
        comment("// CUDA asynchronuous memory copies")
        cuda_text = "#ifdef #{cuda_shared_async}" + endl
        cuda_text += "  // commits displ/hprime.. copies to pipeline stage" + endl
        #cuda_text += "  // pipe.producer_commit();" + endl
        cuda_text += "  cuda::pipeline_producer_commit(pipe, barrier[0]);" + endl
        cuda_text += "  auto token0 = barrier[0].arrive();" + endl
        cuda_text += "#endif  // #{cuda_shared_async}" + endl + endl
        get_output.puts cuda_text
      end

      elem_per_thread.times { |elem_index|
        if elem_per_thread > 1 then
          print tx === get_local_id(0) + ngll3_padded * elem_index / elem_per_thread
        end
        print k === tx/ngll2
        print j === (tx-k*ngll2)/ngllx
        print i === tx - k*ngll2 - j*ngllx
        comment()

        # asynchronuous memory copies modification
        if use_cuda_shared_async then
          comment("// CUDA asynchronuous memory copies")
          get_output.puts "#ifdef #{cuda_shared_async}"

          if (elem_index == 0) then
            cuda_text = "  // next stage to read muv/kappav/.." + endl
            cuda_text += "  pipe.producer_acquire();" + endl + endl
            get_output.puts cuda_text
          end

          print If(active[elem_index]) {
            if only_anisotropy then
              # aniso fields c11store/.. not tested yet to move to shared memory
              get_output.puts "    offset = tx + (NGLL3_PADDED * working_element);"
              cuda_text = "    // c11,c12,.. for aniso" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore              + tx, d_c11store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3      + tx, d_c12store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 2  + tx, d_c13store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 3  + tx, d_c14store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 4  + tx, d_c15store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 5  + tx, d_c16store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 6  + tx, d_c22store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 7  + tx, d_c23store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 8  + tx, d_c24store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 9  + tx, d_c25store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 10 + tx, d_c26store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 11 + tx, d_c33store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 12 + tx, d_c34store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 13 + tx, d_c35store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 14 + tx, d_c36store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 15 + tx, d_c44store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 16 + tx, d_c45store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 17 + tx, d_c46store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 18 + tx, d_c55store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 19 + tx, d_c56store + offset, sizeof(float), pipe);" + endl
              cuda_text += "    cuda::memcpy_async(thread, sh_cstore + NGLL3 * 20 + tx, d_c66store + offset, sizeof(float), pipe);" + endl
              cuda_text += "" + endl
              get_output.puts cuda_text

              # muv only needed for attenuation
              print If(Expression("&&", attenuation, !partial_phys_dispersion_only)) {
                cuda_text = "      // muv for attenuation" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_mul + tx, d_muvstore + offset, sizeof(float), pipe);" + endl
                get_output.puts cuda_text
              }
            else
              # iso/tiso fields
            #print If( ! anisotropy) {
              get_output.puts "    offset = tx + (NGLL3_PADDED * working_element);"
              print If( ! d_ispec_is_tiso[working_element] => lambda {
                cuda_text = "      // isotropic" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_mul + tx, d_muvstore + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_mul + NGLL3 + tx, d_kappavstore + offset, sizeof(float), pipe);" + endl
                get_output.puts cuda_text
              }, :else => lambda {
                cuda_text = "      // tiso" + endl

                # original tiso routine only
                #cuda_text += "        cuda::memcpy_async(thread, sh_mul + tx, d_muvstore + offset, sizeof(float), pipe);" + endl
                #cuda_text += "        cuda::memcpy_async(thread, sh_mul + NGLL3 + tx, d_muhstore + offset, sizeof(float), pipe);" + endl
                #cuda_text += "        cuda::memcpy_async(thread, sh_mul + NGLL3 * 2 + tx, d_kappavstore + offset, sizeof(float), pipe);" + endl
                #cuda_text += "        cuda::memcpy_async(thread, sh_mul + NGLL3 * 3 + tx, d_kappahstore + offset, sizeof(float), pipe);" + endl
                #cuda_text += "        cuda::memcpy_async(thread, sh_mul + NGLL3 * 4 + tx, d_eta_anisostore + offset, sizeof(float), pipe);" + endl
                #cuda_text += "" + endl
                #cuda_text += "        cuda::memcpy_async(thread, ((float3 *)sh_rstore) + tx, ((float3 *)d_rstore) + iglob_1, sizeof(float3), pipe);" + endl
                #cuda_text += "" + endl

                cuda_text += "      // c11,c12,.. for tiso" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore              + tx, d_c11store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3      + tx, d_c12store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 2  + tx, d_c13store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 3  + tx, d_c14store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 4  + tx, d_c15store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 5  + tx, d_c16store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 6  + tx, d_c22store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 7  + tx, d_c23store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 8  + tx, d_c24store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 9  + tx, d_c25store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 10 + tx, d_c26store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 11 + tx, d_c33store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 12 + tx, d_c34store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 13 + tx, d_c35store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 14 + tx, d_c36store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 15 + tx, d_c44store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 16 + tx, d_c45store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 17 + tx, d_c46store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 18 + tx, d_c55store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 19 + tx, d_c56store + offset, sizeof(float), pipe);" + endl
                cuda_text += "      cuda::memcpy_async(thread, sh_cstore + NGLL3 * 20 + tx, d_c66store + offset, sizeof(float), pipe);" + endl
                cuda_text += "" + endl
                get_output.puts cuda_text

                # muv only needed for attenuation
                print If(Expression("&&", attenuation, !partial_phys_dispersion_only)) {
                  cuda_text = "        // muv for attenuation" + endl
                  cuda_text += "        cuda::memcpy_async(thread, sh_mul + tx, d_muvstore + offset, sizeof(float), pipe);" + endl
                  get_output.puts cuda_text
                }
              })
            #}
            end # only_anisotropy

            print If(gravity) {
              cuda_text = "      cuda::memcpy_async(thread, sh_wgll_cube + tx, wgll_cube + tx, sizeof(float), pipe);" + endl
              cuda_text += "      cuda::memcpy_async(thread, sh_gravity_H + 6 * tx, d_gravity_H + 6 * iglob_1, 6 * sizeof(float), pipe);" + endl
              cuda_text += "      cuda::memcpy_async(thread, ((float3 *)sh_gravity_pre_store) + tx, ((float3 *)d_gravity_pre_store) + iglob_1, sizeof(float3), pipe);" + endl
              get_output.puts cuda_text
            }

            print If(Expression("&&", attenuation, !partial_phys_dispersion_only)) {
              print If(!use_lddrk) {
                cuda_text = "        cuda::memcpy_async(thread, sh_epsilondev_xx + tx, epsilondev_xx + tx + (NGLL3 * working_element), sizeof(float), pipe);" + endl
                cuda_text += "        cuda::memcpy_async(thread, sh_epsilondev_xy + tx, epsilondev_xy + tx + (NGLL3 * working_element), sizeof(float), pipe);" + endl
                cuda_text += "        cuda::memcpy_async(thread, sh_epsilondev_xz + tx, epsilondev_xz + tx + (NGLL3 * working_element), sizeof(float), pipe);" + endl
                cuda_text += "        cuda::memcpy_async(thread, sh_epsilondev_yy + tx, epsilondev_yy + tx + (NGLL3 * working_element), sizeof(float), pipe);" + endl
                cuda_text += "        cuda::memcpy_async(thread, sh_epsilondev_yz + tx, epsilondev_yz + tx + (NGLL3 * working_element), sizeof(float), pipe);" + endl
                cuda_text += "" + endl
                cuda_text += "        if (tx < N_SLS) {" + endl
                cuda_text += "          cuda::memcpy_async(thread, sh_alphaval + tx, alphaval + tx, sizeof(float), pipe);" + endl
                cuda_text += "          cuda::memcpy_async(thread, sh_betaval + tx, betaval + tx, sizeof(float), pipe);" + endl
                cuda_text += "          cuda::memcpy_async(thread, sh_gammaval + tx, gammaval + tx, sizeof(float), pipe);" + endl
                cuda_text += "        }" + endl
                cuda_text += "      }" + endl
                cuda_text += "" + endl
                cuda_text += "      for (int i_sls = 0; i_sls < N_SLS; i_sls += 1) {" + endl
                cuda_text += "        int offset_sh = compute_offset_sh(tx, i_sls);" + endl
                cuda_text += "        int offset_sls = compute_offset(tx, i_sls, working_element);" + endl
                cuda_text += "        cuda::memcpy_async(thread, sh_R_xx + offset_sh, R_xx + offset_sls, sizeof(float), pipe);" + endl
                cuda_text += "        cuda::memcpy_async(thread, sh_R_xy + offset_sh, R_xy + offset_sls, sizeof(float), pipe);" + endl
                cuda_text += "        cuda::memcpy_async(thread, sh_R_xz + offset_sh, R_xz + offset_sls, sizeof(float), pipe);" + endl
                cuda_text += "        cuda::memcpy_async(thread, sh_R_yz + offset_sh, R_yz + offset_sls, sizeof(float), pipe);" + endl
                cuda_text += "        cuda::memcpy_async(thread, sh_R_yy + offset_sh, R_yy + offset_sls, sizeof(float), pipe);" + endl
                cuda_text += "" + endl
                cuda_text += "        if (USE_3D_ATTENUATION_ARRAYS) {" + endl
                cuda_text += "          cuda::memcpy_async(thread, sh_factor_common + offset_sh, factor_common + offset_sls, sizeof(float), pipe);" + endl
                cuda_text += "        } else {" + endl
                cuda_text += "          cuda::memcpy_async(thread, sh_factor_common + offset_sh, factor_common + i_sls + (N_SLS * working_element), sizeof(float), pipe);" + endl
                cuda_text += "        }" + endl
                get_output.puts cuda_text
              }
            }
          }

          cuda_text = "  // commits muv/.. copies to pipeline stage" + endl
          cuda_text += "  cuda::pipeline_producer_commit(pipe, barrier[1]);" + endl
          cuda_text += "  decltype(barrier[1].arrive()) token1;" + endl + endl
          get_output.puts cuda_text

          print If(active[elem_index] => lambda {
            get_output.puts "    token1 = barrier[1].arrive();"
          }, :else => lambda {
            get_output.puts "    barrier[1].arrive_and_drop();"
          })
          comment()

          cuda_text = "  // makes sure we have displ/hprime/..." + endl
          #cuda_text += "  //  pipe.producer_commit();" + endl
          #cuda_text += "  // cuda::pipeline_consumer_wait_prior<1>(pipe);" + endl
          cuda_text += "  barrier[0].wait(std::move(token0));" + endl
          cuda_text += "  pipe.consumer_release();" + endl + endl
          cuda_text += "  // checks if anything to do" + endl
          get_output.puts cuda_text

          print If( ! active[elem_index]) { print Return(nil) }

          get_output.puts "#endif  // #{cuda_shared_async}"
          get_output.puts ""
        end

        print If(active[elem_index]) {
          decl *(tempanl.flatten)
          (0..2).each { |indx1|
            (0..2).each { |indx2|
              print tempanl[indx1][indx2] === 0.0
            }
          }
          for_loop = For(l, 0, ngllx-1, :declit => true) {
            print fac[0] === sh_hprime_xx[l*ngllx + i]
            #print offset === k*ngll2 + j*ngllx + l

            if use_cuda_shared_async then
              get_output.puts "#ifdef #{cuda_shared_async}"
              (0..2).each { |indx1|
                # single array
                print tempanl[indx1][0] === tempanl[indx1][0] + s_dummy_loc_single[indx1,k*ngll2 + j*ngllx + l]*fac[0]
              }
              get_output.puts "#else"
            end
            (0..2).each { |indx1|
              print tempanl[indx1][0] === tempanl[indx1][0] + s_dummy_loc[indx1][k*ngll2 + j*ngllx + l]*fac[0]
              # single array
              #print tempanl[indx1][0] === tempanl[indx1][0] + s_dummy_loc[indx1,k*ngll2 + j*ngllx + l]*fac[0]
            }
            if use_cuda_shared_async then
              get_output.puts "#endif  // #{cuda_shared_async}"
            end

            print fac[1] === sh_hprime_xx[l*ngllx + j]
            #print offset === k*ngll2 + l*ngllx + i
            if use_cuda_shared_async then
              get_output.puts "#ifdef #{cuda_shared_async}"
              (0..2).each { |indx1|
                # single array
                print tempanl[indx1][1] === tempanl[indx1][1] + s_dummy_loc_single[indx1,k*ngll2 + l*ngllx + i]*fac[1]
              }
              get_output.puts "#else"
            end
            (0..2).each { |indx1|
              print tempanl[indx1][1] === tempanl[indx1][1] + s_dummy_loc[indx1][k*ngll2 + l*ngllx + i]*fac[1]
              # single array
              #print tempanl[indx1][1] === tempanl[indx1][1] + s_dummy_loc[indx1,k*ngll2 + l*ngllx + i]*fac[1]
            }
            if use_cuda_shared_async then
              get_output.puts "#endif  // #{cuda_shared_async}"
            end

            print fac[2] === sh_hprime_xx[l*ngllx + k]
            #print offset === l*ngll2 + j*ngllx + i
            if use_cuda_shared_async then
              get_output.puts "#ifdef #{cuda_shared_async}"
              (0..2).each { |indx1|
                # single array
                print tempanl[indx1][2] === tempanl[indx1][2] + s_dummy_loc_single[indx1,l*ngll2 + j*ngllx + i]*fac[2]
              }
              get_output.puts "#else"
            end
            (0..2).each { |indx1|
              print tempanl[indx1][2] === tempanl[indx1][2] + s_dummy_loc[indx1][l*ngll2 + j*ngllx + i]*fac[2]
              # single array
              #print tempanl[indx1][2] === tempanl[indx1][2] + s_dummy_loc[indx1,l*ngll2 + j*ngllx + i]*fac[2]
            }
            if use_cuda_shared_async then
              get_output.puts "#endif  // #{cuda_shared_async}"
            end
          }
          get_output.puts "#ifdef #{manually_unrolled_loops}"
            print for_loop.unroll
          get_output.puts "#else"
            print for_loop
          get_output.puts "#endif"
          comment()

          decl *xil   = [ Real("xixl"),   Real("xiyl"),   Real("xizl")   ]
          decl *etal  = [ Real("etaxl"),  Real("etayl"),  Real("etazl")  ]
          decl *gammal= [ Real("gammaxl"),Real("gammayl"),Real("gammazl")]
          print offset === working_element*ngll3_padded + tx
          (0..2).each { |indx|
            print xil[indx]    === d_xi[indx][offset]
            print etal[indx]   === d_eta[indx][offset]
            print gammal[indx] === d_gamma[indx][offset]
          }
          comment()

          decl *(dudl.flatten)
          (0..2).each { |indx1|
            (0..2).each { |indx2|
              print dudl[indx1][indx2] === xil[indx2]*tempanl[indx1][0] + etal[indx2]*tempanl[indx1][1] + gammal[indx2]*tempanl[indx1][2]
            }
          }
          comment()

          # future check: avoiding these to reduce registers
          unless only_anisotropy then
            decl duxdxl_plus_duydyl = Real("duxdxl_plus_duydyl")
            decl duxdxl_plus_duzdzl = Real("duxdxl_plus_duzdzl")
            decl duydyl_plus_duzdzl = Real("duydyl_plus_duzdzl")
          end

          decl duxdyl_plus_duydxl = Real("duxdyl_plus_duydxl")
          decl duzdxl_plus_duxdzl = Real("duzdxl_plus_duxdzl")
          decl duzdyl_plus_duydzl = Real("duzdyl_plus_duydzl")

          # future check: avoiding these to reduce registers, instead will recompute where needed
          unless only_anisotropy then
            print duxdxl_plus_duydyl === dudl[0][0] + dudl[1][1]
            print duxdxl_plus_duzdzl === dudl[0][0] + dudl[2][2]
            print duydyl_plus_duzdzl === dudl[1][1] + dudl[2][2]
          end

          print duxdyl_plus_duydxl === dudl[0][1] + dudl[1][0]
          print duzdxl_plus_duxdzl === dudl[2][0] + dudl[0][2]
          print duzdyl_plus_duydzl === dudl[2][1] + dudl[1][2]
          comment()

          print If(compute_and_store_strain) {
            decl templ = Real("templ")
            print templ === (dudl[0][0] + dudl[1][1] + dudl[2][2])*0.33333333333333333333333333
            print epsilondev_xx_loc[elem_index] === dudl[0][0] - templ
            print epsilondev_yy_loc[elem_index] === dudl[1][1] - templ

            print epsilondev_xy_loc[elem_index] === duxdyl_plus_duydxl * 0.5
            #print epsilondev_xy_loc[elem_index] === (dudl[0][1] + dudl[1][0]) * 0.5

            print epsilondev_xz_loc[elem_index] === duzdxl_plus_duxdzl * 0.5
            #print epsilondev_xz_loc[elem_index] === (dudl[2][0] + dudl[0][2]) * 0.5

            print epsilondev_yz_loc[elem_index] === duzdyl_plus_duydzl * 0.5
            #print epsilondev_yz_loc[elem_index] === (dudl[2][1] + dudl[1][2]) * 0.5

            print If(nspec_strain_only > 1) {
              print epsilon_trace_over_3[tx + working_element*ngll3] === templ
            }
          }
          comment()

          # asynchronuous memory copies modification
          if use_cuda_shared_async then
            comment("// CUDA asynchronuous memory copies")
            cuda_text = "#ifdef #{cuda_shared_async}" + endl
            cuda_text += "    // makes sure we have sh_mul/.." + endl
            cuda_text += "    barrier[1].wait(std::move(token1));" + endl
            cuda_text += "    pipe.consumer_release();" + endl
            cuda_text += "#endif  // #{cuda_shared_async}" + endl + endl
            get_output.puts cuda_text
          end

          if type == :inner_core then
            # inner core
            #print If(anisotropy => lambda {
            if only_anisotropy then
              # anisotropic elements
              comment("    // anisotropic elements")
              print sub_compute_element_ic_aniso.call( offset,
                                                       d_c11store,d_c12store,d_c13store,d_c33store,d_c44store,
                                                       dudl[0][0], dudl[1][1], dudl[2][2],
                                                       duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                       sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                       sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
            #}, :else => lambda {
            else
              # isotropic elements
              comment("    // isotropic elements")
              print sub_compute_element_ic_iso.call( offset,
                                                     d_kappavstore,d_muvstore,
                                                     dudl[0][0], dudl[1][1], dudl[2][2],
                                                     duxdxl_plus_duydyl, duxdxl_plus_duzdzl, duydyl_plus_duzdzl,
                                                     duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                     sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                     sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
            #})
            end
          elsif type == :crust_mantle then
            # crust mantle
            #print If(anisotropy => lambda {
            if only_anisotropy then
              # anisotropic elements
              # fully anisotropic element
              comment("    // anisotropic elements")
              if use_cuda_shared_async then
                get_output.puts "#ifdef #{cuda_shared_async}"
                print sub_compute_element_cm_aniso_async.call( tx,
                                                               sh_cstore,
                                                               *(dudl.flatten),
                                                               duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                               sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                               sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
                get_output.puts "#else"
              end
              print sub_compute_element_cm_aniso.call( offset,
                                                       *(d_cstore.flatten.reject { |e| e.nil?}),
                                                       *(dudl.flatten),
                                                       duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                       sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                       sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
              if use_cuda_shared_async then
                get_output.puts "#endif  // #{cuda_shared_async}"
              end
            #}, :else => lambda {
            else
              # iso/tiso elements
              print If( ! d_ispec_is_tiso[working_element] => lambda {
                # isotropic element
                comment("      // isotropic elements")
                if use_cuda_shared_async then
                  get_output.puts "#ifdef #{cuda_shared_async}"
                  print sub_compute_element_cm_iso_async.call( tx,
                                                               sh_mul,
                                                               dudl[0][0], dudl[1][1], dudl[2][2],
                                                               duxdxl_plus_duydyl, duxdxl_plus_duzdzl, duydyl_plus_duzdzl,
                                                               duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                               sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                               sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
                  get_output.puts "#else"
                end
                print sub_compute_element_cm_iso.call( offset,
                                                       d_kappavstore,d_muvstore,
                                                       dudl[0][0], dudl[1][1], dudl[2][2],
                                                       duxdxl_plus_duydyl, duxdxl_plus_duzdzl, duydyl_plus_duzdzl,
                                                       duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                       sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                       sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
                if use_cuda_shared_async then
                  get_output.puts "#endif  // #{cuda_shared_async}"
                end

              }, :else => lambda {
                # tiso element
                comment("      // transversely isotropic")
                # old way: uses kappav,muv,kappah,.. to compute c11,c12,..
                #if use_cuda_shared_async then
                #  get_output.puts "#ifdef #{cuda_shared_async}"
                #  print sub_compute_element_cm_tiso_async.call( tx,
                #                                                sh_mul,
                #                                                *(dudl.flatten),
                #                                                duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                #                                                iglob[elem_index],
                #                                                sh_rstore,
                #                                                sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                #                                                sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
                #  get_output.puts "#else"
                #end
                #print sub_compute_element_cm_tiso.call( offset,
                #                                        d_kappavstore, d_muvstore,
                #                                        d_kappahstore, d_muhstore, d_eta_anisostore,
                #                                        *(dudl.flatten),
                #                                        duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                #                                        iglob[elem_index],
                #                                        d_rstore,
                #                                        sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                #                                        sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
                #if use_cuda_shared_async then
                #  get_output.puts "#endif  // #{cuda_shared_async}"
                #end

                # new way: uses pre-computed c11,c12,..
                # note: CPU routine treats tiso element as fully anisotropic
                #       this requires all cij store arrays and adds quite a bit of memory consumption.
                #
                if use_cuda_shared_async then
                  get_output.puts "#ifdef #{cuda_shared_async}"
                  print sub_compute_element_cm_tiso_async.call( tx,
                                                                sh_cstore,
                                                                *(dudl.flatten),
                                                                duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                                sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                                sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
                  get_output.puts "#else"
                end
                print sub_compute_element_cm_tiso.call( offset,
                                                        *(d_cstore.flatten.reject { |e| e.nil?}),
                                                        *(dudl.flatten),
                                                        duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                        sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                        sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
                if use_cuda_shared_async then
                  get_output.puts "#endif  // #{cuda_shared_async}"
                end
              })
            #})
            end # only_anisotropy
          end
          comment()

          print If(Expression("&&", attenuation, !partial_phys_dispersion_only)) {
            if use_cuda_shared_async then
              get_output.puts "#ifdef #{cuda_shared_async}"
              print sub_compute_element_att_stress.call(tx, working_element,\
                                                        sh_r_xx, sh_r_yy, sh_r_xy, sh_r_xz, sh_r_yz,\
                                                        sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,\
                                                        sigma[0][1].address, sigma[0][2].address, sigma[1][2].address)
              get_output.puts "#else"
            end
            print sub_compute_element_att_stress.call(tx, working_element,\
                                                      r_xx, r_yy, r_xy, r_xz, r_yz,\
                                                      sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,\
                                                      sigma[0][1].address, sigma[0][2].address, sigma[1][2].address)
            if use_cuda_shared_async then
              get_output.puts "#endif  // #{cuda_shared_async}"
            end
          }
          comment()

          print sigma[1][0] === sigma[0][1]
          print sigma[2][0] === sigma[0][2]
          print sigma[2][1] === sigma[1][2]
          comment()

          decl jacobianl = Real("jacobianl")
          print jacobianl === Expression("/", 1.0, xil[0]*(etal[1]*gammal[2] - etal[2]*gammal[1])\
                                                 - xil[1]*(etal[0]*gammal[2] - etal[2]*gammal[0])\
                                                 + xil[2]*(etal[0]*gammal[1] - etal[1]*gammal[0]))
          comment()

          print If(gravity) {
            if use_cuda_shared_async then
              get_output.puts "#ifdef #{cuda_shared_async}"
              print sub_compute_element_gravity_async.call(tx, iglob[elem_index],\
                                     sh_gravity_pre_store,\
                                     sh_gravity_H,\
                                     sh_wgll_cube, \
                                     jacobianl,\
                                     s_dummy_loc_single,\
                                     sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,\
                                     sigma[0][1].address, sigma[1][0].address, sigma[0][2].address,\
                                     sigma[2][0].address, sigma[1][2].address, sigma[2][1].address,\
                                     rho_s_H[elem_index][0].address, rho_s_H[elem_index][1].address, rho_s_H[elem_index][2].address)
              get_output.puts "#else"
            end
            print sub_compute_element_gravity.call(tx, iglob[elem_index],\
                                   d_gravity_pre_store,\
                                   d_gravity_H,\
                                   wgll_cube, \
                                   jacobianl,\
                                   s_dummy_loc[0], s_dummy_loc[1], s_dummy_loc[2], \
                                   sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,\
                                   sigma[0][1].address, sigma[1][0].address, sigma[0][2].address,\
                                   sigma[2][0].address, sigma[1][2].address, sigma[2][1].address,\
                                   rho_s_H[elem_index][0].address, rho_s_H[elem_index][1].address, rho_s_H[elem_index][2].address)
            if use_cuda_shared_async then
              get_output.puts "#endif  // #{cuda_shared_async}"
            end
          }
          comment()

          (0..2).each { |indx|
            print s_temp[indx][0][tx] === jacobianl * (sigma[0][indx]*xil[0] + sigma[1][indx]*xil[1] + sigma[2][indx]*xil[2])
          }
          (0..2).each { |indx|
            print s_temp[indx][1][tx] === jacobianl * (sigma[0][indx]*etal[0] + sigma[1][indx]*etal[1] + sigma[2][indx]*etal[2])
          }
          (0..2).each { |indx|
            print s_temp[indx][2][tx] === jacobianl * (sigma[0][indx]*gammal[0] + sigma[1][indx]*gammal[1] + sigma[2][indx]*gammal[2])
          }
        }
      }
      print barrier(:local)
      comment()

      elem_per_thread.times { |elem_index|
        if elem_per_thread > 1 then
          print tx === get_local_id(0) + ngll3_padded * elem_index / elem_per_thread
          print k === tx/ngll2
          print j === (tx-k*ngll2)/ngllx
          print i === tx - k*ngll2 - j*ngllx
          comment()
        end

        print If(active[elem_index]) {
          decl *(tempanl.flatten)
          (0..2).each { |indx1|
            (0..2).each { |indx2|
              print tempanl[indx1][indx2] === 0.0
            }
          }
          for_loop = For(l, 0, ngllx-1, :declit => true) {
            print fac[0] === sh_hprimewgll_xx[i*ngllx + l]
            print offset === k*ngll2 + j*ngllx + l
            (0..2).each { |indx1|
              print tempanl[indx1][0] === tempanl[indx1][0] + s_temp[indx1][0][offset]*fac[0]
            }
            print fac[1] === sh_hprimewgll_xx[j*ngllx + l]
            print offset === k*ngll2 + l*ngllx + i
            (0..2).each { |indx1|
              print tempanl[indx1][1] === tempanl[indx1][1] + s_temp[indx1][1][offset]*fac[1]
            }
            print fac[2] === sh_hprimewgll_xx[k*ngllx + l]
            print offset === l*ngll2 + j*ngllx + i
            (0..2).each { |indx1|
              print tempanl[indx1][2] === tempanl[indx1][2] + s_temp[indx1][2][offset]*fac[2]
            }
          }
          get_output.puts "#ifdef #{manually_unrolled_loops}"
            print for_loop.unroll
          get_output.puts "#else"
            print for_loop
          get_output.puts "#endif"
          comment()

          print fac[0] === d_wgllwgll_yz[k*ngllx+j]
          print fac[1] === d_wgllwgll_xz[k*ngllx+i]
          print fac[2] === d_wgllwgll_xy[j*ngllx+i]
          (0..2).each { |indx|
            print sum_terms[indx] === -(fac[0]*tempanl[indx][0] + fac[1]*tempanl[indx][1] + fac[2]*tempanl[indx][2])
          }
          comment()

          print If(gravity) {
            (0..2).each { |indx|
              print sum_terms[indx] === sum_terms[indx] + rho_s_H[elem_index][indx]
            }
          }
          comment()

          get_output.puts "#ifdef #{use_mesh_coloring}"
            get_output.puts "#ifdef #{use_textures_fields}"
              (0..2).each { |indx|
                print d_accel[indx,iglob[elem_index]] === d_accel_tex[iglob[elem_index]*3+indx] + sum_terms[indx]
              }
            get_output.puts "#else"
              (0..2).each { |indx|
                print d_accel[indx,iglob[elem_index]] === d_accel[indx,iglob[elem_index]] + sum_terms[indx]
              }
            get_output.puts "#endif"
          get_output.puts "#else"
            if type == :inner_core then
              __accel_update = lambda {
                print If(nspec_inner_core > coloring_min_nspec_inner_core => lambda {
                  get_output.puts "#ifdef #{use_textures_fields}"
                    (0..2).each { |indx|
                      print d_accel[indx,iglob[elem_index]] === d_accel_tex[iglob[elem_index]*3+indx] + sum_terms[indx]
                    }
                  get_output.puts "#else"
                    (0..2).each { |indx|
                      print d_accel[indx,iglob[elem_index]] === d_accel[indx,iglob[elem_index]] + sum_terms[indx]
                    }
                  get_output.puts "#endif"
                }, :else => lambda{
                  (0..2).each { |indx|
                    print atomicAdd(d_accel+ iglob[elem_index]*3 +indx, sum_terms[indx])
                  }
                })
              }
            elsif type == :crust_mantle then
              __accel_update = lambda {
                get_output.puts "#ifdef #{use_textures_fields}"
                  (0..2).each { |indx|
                    print d_accel[indx,iglob[elem_index]] === d_accel_tex[iglob[elem_index]*3+indx] + sum_terms[indx]
                  }
                get_output.puts "#else"
                  (0..2).each { |indx|
                    print d_accel[indx,iglob[elem_index]] === d_accel[indx,iglob[elem_index]] + sum_terms[indx]
                  }
                get_output.puts "#endif"
              }
            end
            print If(use_mesh_coloring_gpu => __accel_update, :else => lambda {
              (0..2).each { |indx|
                print atomicAdd(d_accel + iglob[elem_index]*3 + indx, sum_terms[indx])
              }
            })
          get_output.puts "#endif"
          comment()

          # updates R_memory
          print If(Expression("&&", attenuation, !partial_phys_dispersion_only ) ) {
            print If(!use_lddrk => lambda {
              # non-LDDRK update
              if use_cuda_shared_async then
                get_output.puts "#ifdef #{cuda_shared_async}"
                __params = [tx, working_element,\
                            sh_mul, sh_factor_common,\
                            sh_alphaval, sh_betaval, sh_gammaval,\
                            r_xx, r_yy, r_xy, r_xz, r_yz,\
                            sh_r_xx, sh_r_yy, sh_r_xy, sh_r_xz, sh_r_yz,\
                            sh_epsilondev_xx, sh_epsilondev_yy, sh_epsilondev_xy, sh_epsilondev_xz, sh_epsilondev_yz,\
                            epsilondev_xx_loc[elem_index],\
                            epsilondev_yy_loc[elem_index],\
                            epsilondev_xy_loc[elem_index],\
                            epsilondev_xz_loc[elem_index],\
                            epsilondev_yz_loc[elem_index],\
                            use_3d_attenuation_arrays]
                print sub_compute_element_att_memory_async.call( *__params )
                get_output.puts "#else"
              end
              __params = [tx, working_element,\
                          d_muvstore, factor_common,\
                          alphaval, betaval, gammaval,\
                          r_xx, r_yy, r_xy, r_xz, r_yz,\
                          epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz,\
                          epsilondev_xx_loc[elem_index],\
                          epsilondev_yy_loc[elem_index],\
                          epsilondev_xy_loc[elem_index],\
                          epsilondev_xz_loc[elem_index],\
                          epsilondev_yz_loc[elem_index],\
                          use_3d_attenuation_arrays]
              print sub_compute_element_att_memory.call( *__params )
              if use_cuda_shared_async then
                get_output.puts "#endif  // #{cuda_shared_async}"
              end
            }, :else => lambda {
              # LDDRK update
              if use_cuda_shared_async then
                get_output.puts "#ifdef #{cuda_shared_async}"
                __params = [tx, working_element,\
                             sh_mul, sh_factor_common,\
                             tau_sigmainvval,\
                             r_xx, r_yy, r_xy, r_xz, r_yz,\
                             sh_r_xx, sh_r_yy, sh_r_xy, sh_r_xz, sh_r_yz,\
                             r_xx_lddrk, r_yy_lddrk, r_xy_lddrk, r_xz_lddrk, r_yz_lddrk,\
                             alpha_lddrk,beta_lddrk,deltat,\
                             epsilondev_xx_loc[elem_index],\
                             epsilondev_yy_loc[elem_index],\
                             epsilondev_xy_loc[elem_index],\
                             epsilondev_xz_loc[elem_index],\
                             epsilondev_yz_loc[elem_index],\
                             use_3d_attenuation_arrays]
                print sub_compute_element_att_memory_lddrk_async.call( *__params )
                get_output.puts "#else"
              end
              __params = [tx, working_element,\
                           d_muvstore, factor_common,\
                           tau_sigmainvval,\
                           r_xx, r_yy, r_xy, r_xz, r_yz,\
                           r_xx_lddrk, r_yy_lddrk, r_xy_lddrk, r_xz_lddrk, r_yz_lddrk,\
                           alpha_lddrk,beta_lddrk,deltat,\
                           epsilondev_xx_loc[elem_index],\
                           epsilondev_yy_loc[elem_index],\
                           epsilondev_xy_loc[elem_index],\
                           epsilondev_xz_loc[elem_index],\
                           epsilondev_yz_loc[elem_index],\
                           use_3d_attenuation_arrays]
              print sub_compute_element_att_memory_lddrk.call( *__params )
              if use_cuda_shared_async then
                get_output.puts "#endif  // #{cuda_shared_async}"
              end
            })
          }
          comment()

          print If(compute_and_store_strain ) {
            print epsilondev_xx[tx + working_element*ngll3] === epsilondev_xx_loc[elem_index]
            print epsilondev_yy[tx + working_element*ngll3] === epsilondev_yy_loc[elem_index]
            print epsilondev_xy[tx + working_element*ngll3] === epsilondev_xy_loc[elem_index]
            print epsilondev_xz[tx + working_element*ngll3] === epsilondev_xz_loc[elem_index]
            print epsilondev_yz[tx + working_element*ngll3] === epsilondev_yz_loc[elem_index]
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
