module BOAST
  def BOAST::compute_element_att_memory(type, n_gll3 = 125, n_gll3_padded = 128, n_sls = 3 )
    if type == :inner_core then
      function_name = "compute_element_ic_att_memory"
    elsif type == :crust_mantle then
      function_name = "compute_element_cm_att_memory"
    else
      raise "Unsupported type : #{type}!"
    end

    v = []
    v.push tx                = Int( "tx",                :dir => :in)
    v.push working_element   = Int( "working_element",   :dir => :in)
    v.push d_muv             = Real("d_muv",             :dir => :in, :dim => [Dim()] )
    v.push factor_common     = Real("factor_common",     :dir => :in, :dim => [Dim()] )
    v.push alphaval          = Real("alphaval",          :dir => :in, :dim => [Dim()] )
    v.push betaval           = Real("betaval",           :dir => :in, :dim => [Dim()] )
    v.push gammaval          = Real("gammaval",          :dir => :in, :dim => [Dim()] )
    v.push r_xx              = Real("R_xx",              :dir => :inout, :dim => [Dim()] )
    v.push r_yy              = Real("R_yy",              :dir => :inout, :dim => [Dim()] )
    v.push r_xy              = Real("R_xy",              :dir => :inout, :dim => [Dim()] )
    v.push r_xz              = Real("R_xz",              :dir => :inout, :dim => [Dim()] )
    v.push r_yz              = Real("R_yz",              :dir => :inout, :dim => [Dim()] )
    v.push epsilondev_xx     = Real("epsilondev_xx",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yy     = Real("epsilondev_yy",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xy     = Real("epsilondev_xy",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xz     = Real("epsilondev_xz",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yz     = Real("epsilondev_yz",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xx_loc = Real("epsilondev_xx_loc", :dir => :in)
    v.push epsilondev_yy_loc = Real("epsilondev_yy_loc", :dir => :in)
    v.push epsilondev_xy_loc = Real("epsilondev_xy_loc", :dir => :in)
    v.push epsilondev_xz_loc = Real("epsilondev_xz_loc", :dir => :in)
    v.push epsilondev_yz_loc = Real("epsilondev_yz_loc", :dir => :in)
    if type == :crust_mantle then
      v.push d_c44store      = Real("d_c44store",        :dir => :in, :dim => [Dim()])
      v.push anisotropy      = Int( "ANISOTROPY",        :dir => :in)
    end
    v.push use_3d_attenuation_arrays = Int( "USE_3D_ATTENUATION_ARRAYS",    :dir => :in)

    ngll3 = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)
    nsls  = Int("N_SLS", :const => n_sls)

    p = Procedure(function_name, v, [], :local => true) {
      decl offset = Int("offset")
      decl i_sls  = Int("i_sls")
      decl mul = Real("mul")
      decl alphaval_loc = Real("alphaval_loc")
      decl betaval_loc = Real("betaval_loc")
      decl gammaval_loc = Real("gammaval_loc")
      decl factor_loc = Real("factor_loc")
      decl sn = Real("sn")
      decl snp1 = Real("snp1")
      if type == :crust_mantle then
        print If( anisotropy, lambda {
          print mul === d_c44store[tx + ngll3_padded*working_element]
        }, lambda {
          print mul === d_muv[tx + ngll3_padded*working_element]
        })
      else
        print mul === d_muv[tx + ngll3_padded*working_element]
      end
      print For( i_sls, 0, nsls - 1 ) {
        print offset === i_sls + nsls*(tx + ngll3*working_element)
        print If(use_3d_attenuation_arrays, lambda {
            print factor_loc  === mul * factor_common[offset]
        }, lambda {
            print factor_loc  === mul * factor_common[i_sls + nsls*working_element ]
        })
        print alphaval_loc === alphaval[i_sls]
        print  betaval_loc ===  betaval[i_sls]
        print gammaval_loc === gammaval[i_sls]

        [[r_xx, epsilondev_xx, epsilondev_xx_loc],
         [r_yy, epsilondev_yy, epsilondev_yy_loc],
         [r_xy, epsilondev_xy, epsilondev_xy_loc],
         [r_xz, epsilondev_xz, epsilondev_xz_loc],
         [r_yz, epsilondev_yz, epsilondev_yz_loc]].each { |r, epsilondev, epsilondev_loc|
          print sn   === factor_loc * epsilondev[tx + ngll3 * working_element]
          print snp1 === factor_loc * epsilondev_loc
          print r[offset] === alphaval_loc*r[offset] + betaval_loc*sn + gammaval_loc*snp1
        }
      }
    }
    return p
  end
end
