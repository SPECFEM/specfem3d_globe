module BOAST
  def BOAST::compute_vector_gradient_kernel(n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)

    function_name = "compute_vector_gradient_kernel"

    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)
    ngll3 = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    v = []
    v.push ispec = Int("ispec", :dir => :in)
    v.push *s_f = ["x", "y", "z"].collect { |a|
      Real("f#{a}", :dir => :in, :dim => [Dim(ngll3)], :local => true )
    }

    v.push *d_xi = [d_xix = Real("d_xix", :dir => :in, :dim => [Dim()]),
                    d_xiy = Real("d_xiy", :dir => :in, :dim => [Dim()]),
                    d_xiz = Real("d_xiz", :dir => :in, :dim => [Dim()])]
    v.push *d_eta = [d_etax = Real("d_etax", :dir => :in, :dim => [Dim()]),
                     d_etay = Real("d_etay", :dir => :in, :dim => [Dim()]),
                     d_etaz = Real("d_etaz", :dir => :in, :dim => [Dim()])]
    v.push *d_gamma = [d_gammax = Real("d_gammax", :dir => :in, :dim => [Dim()]),
                       d_gammay = Real("d_gammay", :dir => :in, :dim => [Dim()]),
                       d_gammaz = Real("d_gammaz", :dir => :in, :dim => [Dim()])]
    v.push sh_hprime_xx = Real("sh_hprime_xx", :dir => :in, :dim => [Dim(ngll2)], :local => true)
    v.push fgrad_loc = Real("fgrad_loc", :dir => :out, :dim => [Dim(9)], :register => true)

    sub = Procedure(function_name, v, :local => true) {
      decl tx = Int("tx")
      decl k = Int("K")
      decl j = Int("J")
      decl i = Int("I")
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
      dfdl = ["x", "y", "z"].collect { |a1|
        ["x", "y", "z"].collect { |a2|
          Real("df#{a1}d#{a2}l")
        }
      }
      decl *(dfdl.flatten)
      decl *fac = (1..3).collect { |n| Real("fac#{n}") }
      comment()

      print tx === get_local_id(0)
      print k === tx/ngll2
      print j === (tx-k*ngll2)/ngllx
      print i === tx - k*ngll2 - j*ngllx
      comment()

      tempanl.flatten.each { |t|
        print t === 0.0
      }
      comment()

      l = Int("l")
      print For(l, 0, ngllx, :operator => "<", :declit => true) {
        print fac[0] === sh_hprime_xx[l*ngllx + i]
        (0..2).each { |indx|
          print tempanl[indx][0] === tempanl[indx][0] + s_f[indx][k*ngll2 + j*ngllx + l]*fac[0]
        }
        print fac[1] === sh_hprime_xx[l*ngllx + j]
        (0..2).each { |indx|
          print tempanl[indx][1] === tempanl[indx][1] + s_f[indx][k*ngll2 + l*ngllx + i]*fac[1]
        }
        print fac[2] === sh_hprime_xx[l*ngllx + k]
        (0..2).each { |indx|
          print tempanl[indx][2] === tempanl[indx][2] + s_f[indx][l*ngll2 + j*ngllx + i]*fac[2]
        }
      }
      comment()

      print offset === ispec*ngll3_padded + tx
      (0..2).each { |indx|
        print xil[indx] === d_xi[indx][offset]
        print etal[indx] === d_eta[indx][offset]
        print gammal[indx] === d_gamma[indx][offset]
      }
      comment()

      (0..2).each { |indx1|
        (0..2).each { |indx2|
          print dfdl[indx1][indx2] === xil[indx2]*tempanl[indx1][0] + etal[indx2]*tempanl[indx1][1] + gammal[indx2]*tempanl[indx1][2]
        }
      }
      comment()

      print fgrad_loc[0] === dfdl[0][0];
      print fgrad_loc[1] === dfdl[0][1];
      print fgrad_loc[2] === dfdl[0][2];
      print fgrad_loc[3] === dfdl[1][0];
      print fgrad_loc[4] === dfdl[1][1];
      print fgrad_loc[5] === dfdl[1][2];
      print fgrad_loc[6] === dfdl[2][0];
      print fgrad_loc[7] === dfdl[2][1];
      print fgrad_loc[8] === dfdl[2][2];
    }
    return sub
  end

  def BOAST::compute_kappa_mu_hess_kernel(ref = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    function_name = "compute_kappa_mu_hess_kernel"

    v = []
    v.push d_ibool = Int( "d_ibool", :dir => :in, :dim => [Dim()])
    v.push d_veloc = Real("d_veloc", :dir => :in, :dim => [Dim(3), Dim()])
    v.push d_b_veloc = Real("d_b_veloc", :dir => :in, :dim => [Dim(3), Dim()])

    v.push *d_xi = [d_xix   = Real("d_xix", :dir => :in, :dim => [Dim()] ),
                    d_xiy   = Real("d_xiy", :dir => :in, :dim => [Dim()] ),
                    d_xiz   = Real("d_xiz", :dir => :in, :dim => [Dim()] ) ]
    v.push *d_eta = [d_etax = Real("d_etax", :dir => :in, :dim => [Dim()] ),
                     d_etay = Real("d_etay", :dir => :in, :dim => [Dim()] ),
                     d_etaz = Real("d_etaz", :dir => :in, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax = Real("d_gammax", :dir => :in, :dim => [Dim()] ),
                       d_gammay = Real("d_gammay", :dir => :in, :dim => [Dim()] ),
                       d_gammaz = Real("d_gammaz", :dir => :in, :dim => [Dim()] ) ]
    v.push d_hprime_xx = Real("d_hprime_xx", :dir => :in, :dim => [Dim()] )
    v.push deltat = Real("deltat", :dir => :in)
    v.push hess_rho_kl   = Real("hess_rho_kl", :dir => :inout, :dim => [Dim()] )
    v.push hess_kappa_kl = Real("hess_kappa_kl", :dir => :inout, :dim => [Dim()] )
    v.push hess_mu_kl    = Real("hess_mu_kl", :dir => :inout, :dim => [Dim()] )

    v.push nspec = Int("NSPEC", :dir => :in)
    v.push use_source_receiver_hess = Int("USE_SOURCE_RECEIVER_HESSIAN", :dir => :in)

    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)
    ngll3 = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    p = Procedure(function_name, v)

    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif (get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      make_specfem3d_header( :ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded )

      sub_compute_vector_gradient_kernel = compute_vector_gradient_kernel(n_gllx, n_gll2, n_gll3, n_gll3_padded )
      print sub_compute_vector_gradient_kernel

      open p

      decl ispec = Int("ispec")
      decl ijk_ispec = Int("ijk_ispec")
      decl tx = Int("tx")
      decl iglob = Int("iglob")

      decl vgrad = Real("vgrad", :dim => [Dim(9)], :allocate => true)
      decl b_vgrad = Real("b_vgrad", :dim => [Dim(9)], :allocate => true)

      decl *sh_veloc = ["x", "y", "z"].collect { |a|
        Real("sh_veloc#{a}", :local => true, :dim => [Dim(ngll3)] )
      }
      decl *sh_b_veloc = ["x", "y", "z"].collect { |a|
        Real("sh_b_veloc#{a}", :local => true, :dim => [Dim(ngll3)] )
      }
      decl sh_hprime_xx = Real("sh_hprime_xx", :local => true, :dim => [Dim(ngll2)] )

      # local variable
      decl hess_rhol = Real("hess_rhol")
      decl hess_kappal = Real("hess_kappal")
      decl hess_mul = Real("hess_mul")

      comment()

      print ispec === get_group_id(0) + get_group_id(1)*get_num_groups(0)
      print ijk_ispec === get_local_id(0) + ngll3*ispec
      print tx === get_local_id(0)
      comment()

      print If(tx < ngll2) {
        print sh_hprime_xx[tx] === d_hprime_xx[tx]
      }
      comment()

      print If(ispec < nspec) {
        print iglob === d_ibool[ijk_ispec] - 1
        (0..2).each { |indx|
          print sh_veloc[indx][tx] === d_veloc[indx, iglob]
          print sh_b_veloc[indx][tx] === d_b_veloc[indx, iglob]
        }
      }
      print barrier(:local)
      comment()

      print If(ispec < nspec) {

        print sub_compute_vector_gradient_kernel.call(ispec,
                                                      *sh_veloc,
                                                      *d_xi, *d_eta, *d_gamma,
                                                      sh_hprime_xx,
                                                      vgrad)
        comment()

        print If(use_source_receiver_hess => lambda {
          print sub_compute_vector_gradient_kernel.call(ispec, *sh_b_veloc, *d_xi, *d_eta, *d_gamma, sh_hprime_xx,  b_vgrad)
        }, :else => lambda {
          l = Int("l")
          print For(l, 0, 8, :declit => true) {
            print b_vgrad[l] === vgrad[l]
          }
        })
        comment()

        print hess_rhol === (vgrad[0]*b_vgrad[0] + vgrad[4]*b_vgrad[4] + vgrad[8]*b_vgrad[8])
        print hess_kappal === 3.0 * (vgrad[0] + vgrad[4] + vgrad[8]) * (b_vgrad[0] + b_vgrad[4] + b_vgrad[8])
        print hess_mul === (vgrad[1] + vgrad[3]) * (b_vgrad[1] + b_vgrad[3]) \
                         + (vgrad[2] + vgrad[6]) * (b_vgrad[2] + b_vgrad[6]) \
                         + (vgrad[5] + vgrad[7]) * (b_vgrad[5] + b_vgrad[7]) \
                         + ((2.0*vgrad[0]   - vgrad[4]       - vgrad[8]      ) * \
                            (2.0*b_vgrad[0] - b_vgrad[4]     - b_vgrad[8]    )   \
                          + (-vgrad[0]      + 2.0*vgrad[4]   - vgrad[8]      ) * \
                            (-b_vgrad[0]    + 2.0*b_vgrad[4] - b_vgrad[8]    )   \
                          + (-vgrad[0]      - vgrad[4]       + 2.0*vgrad[8]  ) * \
                            (-b_vgrad[0]    - b_vgrad[4]     + 2.0*b_vgrad[8])) * (4.0/9.0)
        comment()

        print hess_rho_kl[ijk_ispec] === hess_rho_kl[ijk_ispec] + deltat * hess_rhol
        print hess_kappa_kl[ijk_ispec] === hess_kappa_kl[ijk_ispec] + deltat * hess_kappal
        print hess_mu_kl[ijk_ispec] === hess_mu_kl[ijk_ispec] + deltat * hess_mul
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
