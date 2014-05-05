module BOAST
  def BOAST::compute_element_strain_undo_att(n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)

    function_name = "compute_element_strain_undo_att"

    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)
    ngll3 = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    v = []
    v.push ispec = Int("ispec", :dir => :in)
    v.push ijk_ispec = Int("ijk_ispec", :dir => :in)
    v.push d_ibool   = Int( "d_ibool",  :dir => :in, :dim => [Dim()] ) #unused
    v.push *s_dummy_loc = ["x", "y", "z"].collect { |a|
      Real("s_dummy#{a}_loc", :dir => :in, :dim => [Dim(ngll3)], :local => true )
    }
    v.push *d_xi =    [d_xix      = Real("d_xix",    :dir => :in, :dim => [Dim()] ), d_xiy    = Real("d_xiy",   :dir => :in, :dim => [Dim()] ), d_xiz    = Real("d_xiz",   :dir => :in, :dim => [Dim()] ) ]
    v.push *d_eta =   [d_etax     = Real("d_etax",                  :dir => :in, :dim => [Dim()] ), d_etay = Real("d_etay",:dir => :in, :dim => [Dim()] ), d_etaz = Real("d_etaz",:dir => :in, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax   = Real("d_gammax",                :dir => :in, :dim => [Dim()] ), d_gammay = Real("d_gammay",:dir => :in, :dim => [Dim()] ), d_gammaz = Real("d_gammaz",:dir => :in, :dim => [Dim()] ) ]
    v.push sh_hprime_xx = Real("sh_hprime_xx", :dir => :in, :dim => [Dim(ngll2)], :local => true )
    v.push epsilondev_loc     = Real("epsilondev_loc",          :dir => :out, :dim => [Dim(5)], :register => true )
    v.push epsilon_trace_over_3  = Real("epsilon_trace_over_3",          :dir => :out, :dim => [Dim(1)], :register => true )


    sub = Procedure(function_name, v, [], :local => true) {
      decl tx = Int("tx")
      decl k = Int("K") 
      decl j = Int("J") 
      decl i = Int("I")
      decl l = Int("l")
      decl offset = Int("offset")
      tempanl = ["x", "y", "z"].collect { |a|
        [ 1, 2, 3 ].collect { |n|
          Real("temp#{a}#{n}l")
        }
      }
      decl *(tempanl.flatten)
      decl *xil   = [ Real("xixl"),   Real("xiyl"),   Real("xizl")   ]
      decl *etal  = [ Real("etaxl"),  Real("etayl"),  Real("etazl")  ]
      decl *gammal= [ Real("gammaxl"),Real("gammayl"),Real("gammazl")]
      dudl = ["x", "y", "z"].collect { |a1|
        ["x", "y", "z"].collect { |a2|
          Real("du#{a1}d#{a2}l")
        }
      }
      decl *(dudl.flatten)
      decl templ= Real("templ")
      decl *fac = (1..3).collect { |n| Real("fac#{n}") }

      print tx === get_local_id(0)
      print k === tx/ngll2
      print j === (tx-k*ngll2)/ngllx
      print i === tx - k*ngll2 - j*ngllx

      tempanl.flatten.each { |t|
        print t === 0.0
      }
      print For(l, 0, ngllx - 1) {
        print fac[0] === sh_hprime_xx[l*ngllx + i]
        (0..2).each { |indx|
          print tempanl[indx][0] === tempanl[indx][0] + s_dummy_loc[indx][k*ngll2 + j*ngllx + l]*fac[0]
        }
        print fac[1] === sh_hprime_xx[l*ngllx + j]
        (0..2).each { |indx|
          print tempanl[indx][1] === tempanl[indx][1] + s_dummy_loc[indx][k*ngll2 + l*ngllx + i]*fac[1]
        }
        print fac[2] === sh_hprime_xx[l*ngllx + k]
        (0..2).each { |indx|
          print tempanl[indx][2] === tempanl[indx][2] + s_dummy_loc[indx][l*ngll2 + j*ngllx + i]*fac[2]
        }
      }
      print offset === ispec*ngll3_padded + tx
      (0..2).each { |indx|
        print xil[indx] === d_xi[indx][offset]
        print etal[indx] === d_eta[indx][offset]
        print gammal[indx] === d_gamma[indx][offset]
      }
      (0..2).each { |indx1|
        (0..2).each { |indx2|
          print dudl[indx1][indx2] === xil[indx2]*tempanl[indx1][0] + etal[indx2]*tempanl[indx1][1] + gammal[indx2]*tempanl[indx1][2]
        }
      }
      print templ === (dudl[0][0] + dudl[1][1] + dudl[2][2]) * 0.33333333333333333333
      print epsilondev_loc[0] === dudl[0][0] - templ;
      print epsilondev_loc[1] === dudl[1][1] - templ;
      print epsilondev_loc[2] === ( dudl[0][1] + dudl[1][0] ) * 0.5;
      print epsilondev_loc[3] === ( dudl[2][0] + dudl[0][2] ) * 0.5;
      print epsilondev_loc[4] === ( dudl[2][1] + dudl[1][2] ) * 0.5;
      print epsilon_trace_over_3.dereference === templ
    }
    return sub
  end
end
