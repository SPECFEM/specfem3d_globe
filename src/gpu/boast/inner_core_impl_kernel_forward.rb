module BOAST
  require "./compute_element_att_stress_helper.rb"
  def BOAST::compute_element_ic_att_stress( n_gll3 = 125, n_sls = 3 )
    return BOAST::compute_element_att_stress( :inner_core, n_gll3, n_sls)
  end

  def BOAST::compute_element_cm_att_stress( n_gll3 = 125, n_sls = 3 )
    return BOAST::compute_element_att_stress( :crust_mantle, n_gll3, n_sls)
  end

  require "./compute_element_att_memory_helper.rb"
  def BOAST::compute_element_ic_att_memory(n_gll3 = 125, n_gll3_padded = 128, n_sls = 3 )
    return BOAST::compute_element_att_memory(:inner_core, n_gll3, n_gll3_padded, n_sls)
  end

  def BOAST::compute_element_cm_att_memory( n_gll3 = 125, n_gll3_padded = 128, n_sls = 3 )
    return compute_element_att_memory( :crust_mantle, n_gll3, n_gll3_padded, n_sls )
  end

  require "./compute_element_gravity_helper.rb"
  def BOAST::compute_element_ic_gravity( n_gll3 = 125, r_earth_km = 6371.0 )
    return BOAST::compute_element_gravity( :inner_core, n_gll3, r_earth_km )
  end

  def BOAST::compute_element_cm_gravity( n_gll3 = 125, r_earth_km = 6371.0 )
    return compute_element_gravity( :crust_mantle, n_gll3, r_earth_km )
  end

  def BOAST::compute_element_cm_aniso
    function_name = "compute_element_cm_aniso"
    v = []
    v.push offset = Int( "offset", :dir => :in)
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
    v.push attenuation = Int( "ATTENUATION", :dir => :in)
    v.push one_minus_sum_beta_use = Real("one_minus_sum_beta_use", :dir => :in )
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

    p = Procedure( function_name, v, [], :local => true ) {
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
      decl mul = Real("mul")
      decl minus_sum_beta = Real("minus_sum_beta")

      (0..5).each { |indx1|
        (0..5).each { |indx2|
          print c[indx1][indx2] === d_cstore[indx1][indx2][offset] unless indx2 < indx1
        }
      }
      print If(attenuation) {
        print minus_sum_beta === one_minus_sum_beta_use - 1.0
        print mul === c[3][3] * minus_sum_beta
        
        print c[0][0] === c[0][0] + mul * 1.33333333333333333333
        print c[0][1] === c[0][1] - mul * 0.66666666666666666666
        print c[0][2] === c[0][2] - mul * 0.66666666666666666666
        print c[1][1] === c[1][1] + mul * 1.33333333333333333333
        print c[1][2] === c[1][2] - mul * 0.66666666666666666666
        print c[2][2] === c[2][2] + mul * 1.33333333333333333333
        print c[3][3] === c[3][3] + mul
        print c[4][4] === c[4][4] + mul
        print c[5][5] === c[5][5] + mul
      }
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

  def BOAST::compute_element_cm_iso
    function_name = "compute_element_cm_iso"
    v = []
    v.push offset = Int( "offset", :dir => :in)
    v.push d_kappavstore          = Real("d_kappavstore",          :dir => :in, :dim => [Dim()])
    v.push d_muvstore             = Real("d_muvstore",             :dir => :in, :dim => [Dim()])
    v.push attenuation            = Int( "ATTENUATION",            :dir => :in)
    v.push one_minus_sum_beta_use = Real("one_minus_sum_beta_use", :dir => :in )
    v.push *dudl = ["x", "y", "z"].collect { |a1|
      Real("du#{a1}d#{a1}l", :dir => :in)
    }
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

    p = Procedure( function_name, v, [], :local => true ) {
      decl lambdal = Real("lambdal"), mul = Real("mul"), lambdalplus2mul = Real("lambdalplus2mul"), kappal = Real("kappal")

      print kappal === d_kappavstore[offset]
      print mul === d_muvstore[offset]
      print If(attenuation) {
        print mul === mul * one_minus_sum_beta_use
      }

      print lambdalplus2mul === kappal + mul * 1.33333333333333333333
      print lambdal === lambdalplus2mul - mul * 2.0

      print sigma_xx.dereference === lambdalplus2mul*dudl[0] + lambdal*duydyl_plus_duzdzl
      print sigma_yy.dereference === lambdalplus2mul*dudl[1] + lambdal*duxdxl_plus_duzdzl
      print sigma_zz.dereference === lambdalplus2mul*dudl[2] + lambdal*duxdxl_plus_duydyl

      print sigma_xy.dereference === mul*duxdyl_plus_duydxl
      print sigma_xz.dereference === mul*duzdxl_plus_duxdzl
      print sigma_yz.dereference === mul*duzdyl_plus_duydzl
    }
    return p
  end

  def BOAST::compute_element_cm_tiso
    function_name = "compute_element_cm_tiso"
    v = []
    v.push offset = Int( "offset", :dir => :in)
    v.push d_kappavstore          = Real("d_kappavstore",          :dir => :in, :dim => [Dim()])
    v.push d_muvstore             = Real("d_muvstore",             :dir => :in, :dim => [Dim()])
    v.push d_kappahstore          = Real("d_kappahstore",          :dir => :in, :dim => [Dim()])
    v.push d_muhstore             = Real("d_muhstore",             :dir => :in, :dim => [Dim()])
    v.push d_eta_anisostore       = Real("d_eta_anisostore",       :dir => :in, :dim => [Dim()])
    v.push attenuation            = Int( "ATTENUATION",            :dir => :in)
    v.push one_minus_sum_beta_use = Real("one_minus_sum_beta_use", :dir => :in )
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
#    v.push nglob                  = Int( "NGLOB",              :dir => :in)
    v.push d_ystore               = Real("d_ystore",           :dir => :in,  :dim => [Dim()])
    v.push d_zstore               = Real("d_zstore",           :dir => :in,  :dim => [Dim()])
    v.push sigma_xx               = Real("sigma_xx",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yy               = Real("sigma_yy",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_zz               = Real("sigma_zz",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xy               = Real("sigma_xy",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_xz               = Real("sigma_xz",           :dir => :out, :dim => [Dim()], :private => true )
    v.push sigma_yz               = Real("sigma_yz",           :dir => :out, :dim => [Dim()], :private => true )

    p = Procedure( function_name, v, [], :local => true ) {
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

      print kappavl === d_kappavstore[offset]
      print muvl === d_muvstore[offset]

      print kappahl === d_kappahstore[offset]
      print muhl === d_muhstore[offset]

      print If(attenuation) {
        print muvl === muvl * one_minus_sum_beta_use
        print muhl === muhl * one_minus_sum_beta_use
      }
      print rhovpvsq === kappavl + muvl * 1.33333333333333333333
      print rhovphsq === kappahl + muhl * 1.33333333333333333333


      print rhovsvsq === muvl
      print rhovshsq === muhl

      print eta_aniso === d_eta_anisostore[offset]

      print theta === d_ystore[iglob]
      print phi === d_zstore[iglob]

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
                    (rhovsvsq*(cosfourphi*-2.0 + cos(phi*4.0 - theta*2.0) - costwotheta*2.0 +
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

  def BOAST::inner_core_impl_kernel_forward(ref = true, elem_per_thread = 1, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, r_earth_km = 6371.0, coloring_min_nspec_inner_core = 1000, i_flag_in_fictitious_cube = 11)
    return BOAST::impl_kernel(:inner_core, true, ref, elem_per_thread, mesh_coloring, textures_fields, textures_constants, unroll_loops, n_gllx, n_gll2, n_gll3, n_gll3_padded, n_sls, r_earth_km, coloring_min_nspec_inner_core, i_flag_in_fictitious_cube)
  end

  def BOAST::impl_kernel(type, forward, ref = true, elem_per_thread = 1, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = false, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, r_earth_km = 6371.0, coloring_min_nspec_inner_core = 1000, i_flag_in_fictitious_cube = 11, launch_bounds = false, min_blocks = 7)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    v = []
    if type == :inner_core then
      function_name = "inner_core_impl_kernel"
    elsif type == :crust_mantle then
      function_name = "crust_mantle_impl_kernel"
    else
      raise "Unsupported_type : #{type}!"
    end
    if forward then
      function_name += "_forward"
    else
      function_name += "_adjoint"
    end
    v.push nb_blocks_to_compute    = Int("nb_blocks_to_compute",     :dir => :in)
    #v.push nglob                   = Int("NGLOB",                    :dir => :in)
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
    v.push *d_xi = [d_xix          = Real("d_xix",                   :dir => :in, :restrict => true, :dim => [Dim()] ), d_xiy = Real("d_xiy",:dir => :in, :restrict => true, :dim => [Dim()] ), d_xiz = Real("d_xiz",:dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push *d_eta = [d_etax        = Real("d_etax",                  :dir => :in, :restrict => true, :dim => [Dim()] ), d_etay = Real("d_etay",:dir => :in, :restrict => true, :dim => [Dim()] ), d_etaz = Real("d_etaz",:dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax    = Real("d_gammax",                :dir => :in, :restrict => true, :dim => [Dim()] ), d_gammay = Real("d_gammay",:dir => :in, :restrict => true, :dim => [Dim()] ), d_gammaz = Real("d_gammaz",:dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push d_hprime_xx             = Real("d_hprime_xx",             :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_hprimewgll_xx         = Real("d_hprimewgll_xx",         :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_wgllwgll_xy           = Real("d_wgllwgll_xy",           :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_wgllwgll_xz           = Real("d_wgllwgll_xz",           :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_wgllwgll_yz           = Real("d_wgllwgll_yz",           :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_kappavstore           = Real("d_kappavstore",           :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_muvstore              = Real("d_muvstore",              :dir => :in, :restrict => true, :dim => [Dim()] )
    if type == :crust_mantle then
      v.push d_kappahstore           = Real("d_kappahstore",           :dir => :in, :restrict => true, :dim => [Dim()] )
      v.push d_muhstore              = Real("d_muhstore",              :dir => :in, :restrict => true, :dim => [Dim()] )
      v.push d_eta_anisostore        = Real("d_eta_anisostore",        :dir => :in, :restrict => true, :dim => [Dim()] )
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
    v.push alphaval                = Real("alphaval",                :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push betaval                 = Real("betaval",                 :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push gammaval                = Real("gammaval",                :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push anisotropy              = Int( "ANISOTROPY",              :dir => :in)
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
    v.push gravity                 = Int( "GRAVITY",                 :dir => :in)
    v.push *d_store = [d_xstore    = Real("d_xstore",                :dir => :in, :restrict => true, :dim => [Dim()] ), d_ystore = Real("d_ystore",:dir => :in, :restrict => true, :dim => [Dim()] ), d_zstore = Real("d_zstore",:dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push d_minus_gravity_table   = Real("d_minus_gravity_table",   :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_minus_deriv_gravity_table = Real("d_minus_deriv_gravity_table", :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_density_table         = Real("d_density_table",         :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push wgll_cube               = Real("wgll_cube",               :dir => :in, :restrict => true, :dim => [Dim()] )
    if type == :inner_core then
      v.push nspec_strain_only = Int( "NSPEC_INNER_CORE_STRAIN_ONLY", :dir => :in)
      v.push nspec_inner_core        = Int( "NSPEC_INNER_CORE",        :dir => :in)
    elsif type == :crust_mantle then
      v.push nspec_strain_only = Int( "NSPEC_CRUST_MANTLE_STRAIN_ONLY", :dir => :in)
    end

    ngllx        = Int("NGLLX", :const => n_gllx)
    ngll2        = Int("NGLL2", :const => n_gll2)
    ngll3        = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)
    if type == :inner_core then
      iflag_in_fictitious_cube = Int("IFLAG_IN_FICTITIOUS_CUBE", :const => i_flag_in_fictitious_cube)
    end

    use_mesh_coloring       = Int("USE_MESH_COLORING_GPU",   :const => mesh_coloring)
    use_textures_constants  = Int("USE_TEXTURES_CONSTANTS",  :const => textures_constants)
    use_textures_fields     = Int("USE_TEXTURES_FIELDS",     :const => textures_fields)
    manually_unrolled_loops = Int("MANUALLY_UNROLLED_LOOPS", :const => unroll_loops)
    use_launch_bounds       = Int("USE_LAUNCH_BOUNDS",       :const => launch_bounds)
    launch_min_blocks       = Int("LAUNCH_MIN_BLOCKS",       :const => min_blocks)

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

    if (get_lang == CUDA) then
      qualifiers = "\n#ifdef #{use_launch_bounds}\n__launch_bounds__(#{ngll3_padded}, #{launch_min_blocks})\n#endif\n"
    elsif(get_lang == CL) then
      qualifiers = "" # "__attribute__((reqd_work_group_size(#{ngll3_padded},1,1))) " # (inefficient)
    end

    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu".gsub("_forward","").gsub("_adjoint",""))
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header(:ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded, :n_sls => n_sls, :r_earth_km => r_earth_km, :coloring_min_nspec_inner_core => coloring_min_nspec_inner_core, :iflag_in_fictitious_cube => i_flag_in_fictitious_cube)
      if type == :inner_core then
        #DEACTIVATE USE TEXTURES CONSTANTS
        get_output.puts "#ifdef #{use_textures_constants}"
        get_output.puts "#undef #{use_textures_constants}"
        get_output.puts "#endif"
      end
#      if get_lang == CUDA then
#        get_output.puts "#ifdef #{use_textures_fields}"
#          decl d_displ_tex
#          decl d_accel_tex
#        get_output.puts "#endif"
#        get_output.puts "#ifdef #{use_textures_constants}"
#          decl d_hprime_xx_tex
#          decl d_hprimewgll_xx_tex
#        get_output.puts "#endif"
#      end
      if type == :inner_core then
        sub_compute_element_att_stress =  compute_element_ic_att_stress(n_gll3, n_sls)
        sub_compute_element_att_memory =  compute_element_ic_att_memory(n_gll3, n_gll3_padded, n_sls)
        sub_compute_element_gravity =  compute_element_ic_gravity(n_gll3, r_earth_km)
      elsif type == :crust_mantle then
        sub_compute_element_att_stress =  compute_element_cm_att_stress(n_gll3, n_sls)
        sub_compute_element_att_memory =  compute_element_cm_att_memory(n_gll3, n_gll3_padded, n_sls)
        sub_compute_element_gravity =  compute_element_cm_gravity(n_gll3, r_earth_km)
      end
      print sub_compute_element_att_stress
      print sub_compute_element_att_memory
      print sub_compute_element_gravity
      if type == :crust_mantle then
        sub_compute_element_cm_aniso = compute_element_cm_aniso
        print sub_compute_element_cm_aniso
        sub_compute_element_cm_iso = compute_element_cm_iso
        print sub_compute_element_cm_iso
        sub_compute_element_cm_tiso = compute_element_cm_tiso
        print sub_compute_element_cm_tiso
      end

      if get_lang == CL then
        get_output.puts "#ifdef #{use_textures_fields}"
          get_output.puts "#ifdef #{use_textures_constants}"
            p = Procedure(function_name, v+textures_fields+textures_constants, constants, :qualifiers => qualifiers)
            open p
            set_indent_level(0)
          get_output.puts "#else"
            p = Procedure(function_name, v+textures_fields, constants, :qualifiers => qualifiers)
            open p
            set_indent_level(0)
          get_output.puts "#endif"
        get_output.puts "#else"
          get_output.puts "#ifdef #{use_textures_constants}"
            p = Procedure(function_name, v+textures_constants, constants, :qualifiers => qualifiers)
            open p
            set_indent_level(0)
          get_output.puts "#else"
            p = Procedure(function_name, v, constants, :qualifiers => qualifiers)
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
        p = Procedure(function_name, v, constants, :qualifiers => qualifiers)
        open p
      end

        decl bx = Int("bx")
        decl tx = Int("tx")
        decl k  = Int("K"), j = Int("J"), i = Int("I")
        get_output.puts "#ifndef #{manually_unrolled_loops}"
          decl l = Int("l")
        get_output.puts "#endif"
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
        decl *(tempanl.flatten)
        decl *xil   = [ Real("xixl"),   Real("xiyl"),   Real("xizl")   ]
        decl *etal  = [ Real("etaxl"),  Real("etayl"),  Real("etazl")  ]
        decl *gammal= [ Real("gammaxl"),Real("gammayl"),Real("gammazl")]
        decl jacobianl= Real("jacobianl")
        dudl = ["x", "y", "z"].collect { |a1|
          ["x", "y", "z"].collect { |a2|
            Real("du#{a1}d#{a2}l")
          }
        }
        decl *(dudl.flatten)
        decl duxdxl_plus_duydyl = Real("duxdxl_plus_duydyl")
        decl duxdxl_plus_duzdzl = Real("duxdxl_plus_duzdzl")
        decl duydyl_plus_duzdzl = Real("duydyl_plus_duzdzl")
        decl duxdyl_plus_duydxl = Real("duxdyl_plus_duydxl")
        decl duzdxl_plus_duxdzl = Real("duzdxl_plus_duxdzl")
        decl duzdyl_plus_duydzl = Real("duzdyl_plus_duydzl")
        decl templ = Real("templ")
        decl *fac = [1,2,3].collect {|n|
          Real("fac#{n}")
        }
        if type == :inner_core then
          decl lambdal = Real("lambdal")
          decl mul = Real("mul")
          decl lambdalplus2mul = Real("lambdalplus2mul")
          decl kappal = Real("kappal")
          decl mul_iso = Real("mul_iso")
          decl mul_aniso = Real("mul_aniso")
        elsif type == :crust_mantle then
          decl one_minus_sum_beta_use = Real("one_minus_sum_beta_use")
        end

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
        if type == :inner_core then
          decl c11 = Real("c11")
          decl c12 = Real("c12")
          decl c13 = Real("c13")
          decl c33 = Real("c33")
          decl c44 = Real("c44")
        end
        decl *sum_terms = [1,2,3].collect {|n|
          Real("sum_terms#{n}")
        }
        rho_s_H = (1..elem_per_thread).collect{ |e_i|
          [1,2,3].collect {|n|
            Real("rho_s_H_#{e_i}_#{n}")
          }
        }
        decl *(rho_s_H.flatten)
  
        decl *s_dummy_loc = ["x", "y", "z"].collect { |a|
          Real("s_dummy#{a}_loc", :local => true, :dim => [Dim(ngll3)] )
        }
  
        s_temp = ["x", "y", "z"].collect { |a|
          [ 1, 2, 3 ].collect { |n|
            Real("s_temp#{a}#{n}", :local => true, :dim => [Dim(ngll3)] )
          }
        }
        decl *(s_temp.flatten)
  
        decl sh_hprime_xx     = Real("sh_hprime_xx",     :local => true, :dim => [Dim(ngll2)] )
        decl sh_hprimewgll_xx = Real("sh_hprimewgll_xx", :local => true, :dim => [Dim(ngll2)] )
  
        print bx === get_group_id(1)*get_num_groups(0)+get_group_id(0)
elem_per_thread.times { |elem_index| 
        print tx === get_local_id(0) + ngll3_padded * elem_index / elem_per_thread
        print active[elem_index] === Ternary( Expression("&&", tx < ngll3, bx < nb_blocks_to_compute), 1, 0)

        print If(active[elem_index]) {
  if elem_index == 0 then
          get_output.puts "#ifdef #{use_mesh_coloring}"
            print working_element === bx
          get_output.puts "#else"
            print If(use_mesh_coloring_gpu, lambda {
              print working_element === bx
            }, lambda {
              print working_element === d_phase_ispec_inner[bx + num_phase_ispec*(d_iphase-1)]-1
            })
          get_output.puts "#endif"
  end
          __texture_fetch = lambda {
            print iglob[elem_index] === d_ibool[working_element*ngll3 + tx]-1

            get_output.puts "#ifdef #{use_textures_fields}"
              (0..2).each { |indx|
                print s_dummy_loc[indx][tx] === d_displ_tex[iglob[elem_index]*3+indx]
              }
            get_output.puts "#else"
              (0..2).each { |indx|
                print s_dummy_loc[indx][tx] === d_displ[indx, iglob[elem_index]]
              }
            get_output.puts "#endif"
          }
          if type == :inner_core then
            print If(d_idoubling[working_element] == iflag_in_fictitious_cube, lambda {
              print active[elem_index] === 0
            }, __texture_fetch )
          elsif type == :crust_mantle then
            __texture_fetch.call
          end
        }
        #inner core and crust mantle differ here, but crust mantle implementation though more reccent seems odd...
        print If(tx < ngll2) {
          get_output.puts "#ifdef #{use_textures_constants}"
            print sh_hprime_xx[tx] === d_hprime_xx_tex[tx]
            print sh_hprimewgll_xx[tx] === d_hprimewgll_xx_tex[tx]
          get_output.puts "#else"
            print sh_hprime_xx[tx] === d_hprime_xx[tx]
            print sh_hprimewgll_xx[tx] === d_hprimewgll_xx[tx]
          get_output.puts "#endif"
        }
}
        print barrier(:local)
  
elem_per_thread.times { |elem_index|
  if elem_per_thread > 1 then
        print tx === get_local_id(0) + ngll3_padded * elem_index / elem_per_thread
  end
        print k === tx/ngll2
        print j === (tx-k*ngll2)/ngllx
        print i === tx - k*ngll2 - j*ngllx

        print If(active[elem_index]) {
          (0..2).each { |indx1|
            (0..2).each { |indx2|
              print tempanl[indx1][indx2] === 0.0
            }
          }
          for_loop = For(l, 0, ngllx-1) {
            print fac[0] === sh_hprime_xx[l*ngllx + i]
            #print offset === k*ngll2 + j*ngllx + l
            (0..2).each { |indx1|
              print tempanl[indx1][0] === tempanl[indx1][0] + s_dummy_loc[indx1][k*ngll2 + j*ngllx + l]*fac[0]
            }
            print fac[1] === sh_hprime_xx[l*ngllx + j]
            #print offset === k*ngll2 + l*ngllx + i
            (0..2).each { |indx1|
              print tempanl[indx1][1] === tempanl[indx1][1] + s_dummy_loc[indx1][k*ngll2 + l*ngllx + i]*fac[1]
            }
            print fac[2] === sh_hprime_xx[l*ngllx + k]
            #print offset === l*ngll2 + j*ngllx + i
            (0..2).each { |indx1|
              print tempanl[indx1][2] === tempanl[indx1][2] + s_dummy_loc[indx1][l*ngll2 + j*ngllx + i]*fac[2]
            }
          }
          get_output.puts "#ifdef #{manually_unrolled_loops}"
            for_loop.unroll
          get_output.puts "#else"
            print for_loop
          get_output.puts "#endif"

          print offset === working_element*ngll3_padded + tx
          (0..2).each { |indx|
            print xil[indx]    === d_xi[indx][offset]
            print etal[indx]   === d_eta[indx][offset]
            print gammal[indx] === d_gamma[indx][offset]
          }
  
          (0..2).each { |indx1|
            (0..2).each { |indx2|
              print dudl[indx1][indx2] === xil[indx2]*tempanl[indx1][0] + etal[indx2]*tempanl[indx1][1] + gammal[indx2]*tempanl[indx1][2]
            }
          }
          print duxdxl_plus_duydyl === dudl[0][0] + dudl[1][1]
          print duxdxl_plus_duzdzl === dudl[0][0] + dudl[2][2]
          print duydyl_plus_duzdzl === dudl[1][1] + dudl[2][2]
          print duxdyl_plus_duydxl === dudl[0][1] + dudl[1][0]
          print duzdxl_plus_duxdzl === dudl[2][0] + dudl[0][2]
          print duzdyl_plus_duydzl === dudl[2][1] + dudl[1][2]
  
          print If(compute_and_store_strain) {
            print templ === (dudl[0][0] + dudl[1][1] + dudl[2][2])*0.33333333333333333333333333
            print epsilondev_xx_loc[elem_index] === dudl[0][0] - templ
            print epsilondev_yy_loc[elem_index] === dudl[1][1] - templ
            print epsilondev_xy_loc[elem_index] === duxdyl_plus_duydxl * 0.5
            print epsilondev_xz_loc[elem_index] === duzdxl_plus_duxdzl * 0.5
            print epsilondev_yz_loc[elem_index] === duzdyl_plus_duydzl * 0.5
            print If(nspec_strain_only == 1, lambda {
              print epsilon_trace_over_3[tx] === templ
            }, lambda {
              print epsilon_trace_over_3[tx + working_element*ngll3] === templ
            })
          }
  
          if type == :inner_core then
            print kappal === d_kappavstore[offset]
            print mul === d_muvstore[offset]
            print If(attenuation, lambda {
              print If(use_3d_attenuation_arrays, lambda {
                print mul_iso  === mul * one_minus_sum_beta[tx+working_element*ngll3]
                print mul_aniso === mul * ( one_minus_sum_beta[tx+working_element*ngll3] - 1.0 )
              }, lambda {
                print mul_iso  === mul * one_minus_sum_beta[working_element]
                print mul_aniso === mul * ( one_minus_sum_beta[working_element] - 1.0 )
              })
            }, lambda {
              print mul_iso === mul
            })
            print If(anisotropy, lambda {
              print c11 === d_c11store[offset]
              print c12 === d_c12store[offset]
              print c13 === d_c13store[offset]
              print c33 === d_c33store[offset]
              print c44 === d_c44store[offset]
              print If(attenuation) {
                print c11 === c11 + mul_aniso * 1.33333333333333333333
                print c12 === c12 - mul_aniso * 0.66666666666666666666
                print c13 === c13 - mul_aniso * 0.66666666666666666666
                print c33 === c33 + mul_aniso * 1.33333333333333333333
                print c44 === c44 + mul_aniso
              }
              print sigma[0][0] === c11*dudl[0][0] + c12*dudl[1][1] + c13*dudl[2][2]
              print sigma[1][1] === c12*dudl[0][0] + c11*dudl[1][1] + c13*dudl[2][2]
              print sigma[2][2] === c13*dudl[0][0] + c13*dudl[1][1] + c33*dudl[2][2]
              print sigma[0][1] === (c11-c12)*duxdyl_plus_duydxl*0.5
              print sigma[0][2] === c44*duzdxl_plus_duxdzl
              print sigma[1][2] === c44*duzdyl_plus_duydzl
            }, lambda {
              print lambdalplus2mul === kappal + mul_iso * 1.33333333333333333333
              print lambdal === lambdalplus2mul - mul_iso * 2.0
    
              print sigma[0][0] === lambdalplus2mul*dudl[0][0] + lambdal*duydyl_plus_duzdzl
              print sigma[1][1] === lambdalplus2mul*dudl[1][1] + lambdal*duxdxl_plus_duzdzl
              print sigma[2][2] === lambdalplus2mul*dudl[2][2] + lambdal*duxdxl_plus_duydyl
    
              print sigma[0][1] === mul*duxdyl_plus_duydxl
              print sigma[0][2] === mul*duzdxl_plus_duxdzl
              print sigma[1][2] === mul*duzdyl_plus_duydzl
            })
          elsif type == :crust_mantle then
            print If(attenuation) {
              print If(use_3d_attenuation_arrays, lambda {
                print one_minus_sum_beta_use === one_minus_sum_beta[tx+working_element*ngll3]
              }, lambda {
                print one_minus_sum_beta_use === one_minus_sum_beta[working_element]
              })
            }
            print If(anisotropy, lambda {
              print sub_compute_element_cm_aniso.call( offset,
                                                      *(d_cstore.flatten.reject { |e| e.nil?}),
                                                      attenuation,
                                                      one_minus_sum_beta_use,
                                                      *(dudl.flatten),
                                                      duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                      sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                      sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
            }, lambda {
              print If( ! d_ispec_is_tiso[working_element], lambda {
                print sub_compute_element_cm_iso.call( offset,
                                                       d_kappavstore,d_muvstore,
                                                       attenuation,
                                                       one_minus_sum_beta_use,
                                                       dudl[0][0], dudl[1][1], dudl[2][2],
                                                       duxdxl_plus_duydyl, duxdxl_plus_duzdzl, duydyl_plus_duzdzl,
                                                       duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                       sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                       sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
              }, lambda {
                print sub_compute_element_cm_tiso.call( offset,
                                                   d_kappavstore, d_muvstore,
                                                   d_kappahstore, d_muhstore, d_eta_anisostore,
                                                   attenuation,
                                                   one_minus_sum_beta_use,
                                                   *(dudl.flatten),
                                                   duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,
                                                   iglob[elem_index], #nglob,
                                                   d_store[1], d_store[2],
                                                   sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,
                                                   sigma[0][1].address, sigma[0][2].address, sigma[1][2].address )
              })
            })
          end


  
          print If(Expression("&&", attenuation, !partial_phys_dispersion_only)) {
            print sub_compute_element_att_stress.call(tx, working_element,\
                                   r_xx, r_yy, r_xy, r_xz, r_yz,\
                                   sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,\
                                   sigma[0][1].address, sigma[0][2].address, sigma[1][2].address)
          }
          print sigma[1][0] === sigma[0][1]
          print sigma[2][0] === sigma[0][2]
          print sigma[2][1] === sigma[1][2]
  
          print jacobianl === Expression("/", 1.0, xil[0]*(etal[1]*gammal[2] - etal[2]*gammal[1])\
                                                 - xil[1]*(etal[0]*gammal[2] - etal[2]*gammal[0])\
                                                 + xil[2]*(etal[0]*gammal[1] - etal[1]*gammal[0]))
          print If(gravity) {
            print sub_compute_element_gravity.call(tx, iglob[elem_index],\
                                   d_store[0], d_store[1], d_store[2],\
                                   d_minus_gravity_table, d_minus_deriv_gravity_table, d_density_table,\
                                   wgll_cube, jacobianl,\
                                   s_dummy_loc[0], s_dummy_loc[1], s_dummy_loc[2],\
                                   sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,\
                                   sigma[0][1].address, sigma[1][0].address, sigma[0][2].address,\
                                   sigma[2][0].address, sigma[1][2].address, sigma[2][1].address,\
                                   rho_s_H[elem_index][0].address, rho_s_H[elem_index][1].address, rho_s_H[elem_index][2].address)
          }
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
elem_per_thread.times { |elem_index|
  if elem_per_thread > 1 then
        print tx === get_local_id(0) + ngll3_padded * elem_index / elem_per_thread
        print k === tx/ngll2
        print j === (tx-k*ngll2)/ngllx
        print i === tx - k*ngll2 - j*ngllx
  end

        print If(active[elem_index]) {
          (0..2).each { |indx1|
            (0..2).each { |indx2|
              print tempanl[indx1][indx2] === 0.0
            }
          }
          for_loop = For(l, 0, ngllx-1) {
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
            for_loop.unroll
          get_output.puts "#else"
            print for_loop
          get_output.puts "#endif"
          print fac[0] === d_wgllwgll_yz[k*ngllx+j]
          print fac[1] === d_wgllwgll_xz[k*ngllx+i]
          print fac[2] === d_wgllwgll_xy[j*ngllx+i]
          (0..2).each { |indx|
            print sum_terms[indx] === -(fac[0]*tempanl[indx][0] + fac[1]*tempanl[indx][1] + fac[2]*tempanl[indx][2])
          }
  
          print If(gravity) {
            (0..2).each { |indx|
              print sum_terms[indx] === sum_terms[indx] + rho_s_H[elem_index][indx]
            }
          }
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
                print If(nspec_inner_core > coloring_min_nspec_inner_core, lambda {
                  get_output.puts "#ifdef #{use_textures_fields}"
                    (0..2).each { |indx|
                      print d_accel[indx,iglob[elem_index]] === d_accel_tex[iglob[elem_index]*3+indx] + sum_terms[indx]
                    }
                  get_output.puts "#else"
                    (0..2).each { |indx|
                      print d_accel[indx,iglob[elem_index]] === d_accel[indx,iglob[elem_index]] + sum_terms[indx]
                    }
                  get_output.puts "#endif"
                }, lambda{
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
            print If(use_mesh_coloring_gpu, __accel_update, lambda {
              (0..2).each { |indx|
                print atomicAdd(d_accel + iglob[elem_index]*3 + indx, sum_terms[indx])
              }
            })
          get_output.puts "#endif"
          print If(Expression("&&", attenuation, !partial_phys_dispersion_only ) ) {
            __params = [tx, working_element,\
                        d_muvstore, factor_common,\
                        alphaval, betaval, gammaval,\
                        r_xx, r_yy, r_xy, r_xz, r_yz,\
                        epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz,\
                        epsilondev_xx_loc[elem_index],\
                        epsilondev_yy_loc[elem_index],\
                        epsilondev_xy_loc[elem_index],\
                        epsilondev_xz_loc[elem_index],\
                        epsilondev_yz_loc[elem_index]]
            if type == :crust_mantle then
              __params += [d_cstore[3][3], anisotropy]
            end
            __params.push use_3d_attenuation_arrays
            print sub_compute_element_att_memory.call( *__params )
          }
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
