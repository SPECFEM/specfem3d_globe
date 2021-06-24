module BOAST
  def BOAST::compute_element_att_memory(type, n_gll3 = 125, n_gll3_padded = 128, n_sls = 3 )
    if type == :inner_core then
      function_name = "compute_element_ic_att_memory"
    elsif type == :inner_core_lddrk then
      function_name = "compute_element_ic_att_memory_lddrk"
    elsif type == :crust_mantle then
      function_name = "compute_element_cm_att_memory"
    elsif type == :crust_mantle_lddrk then
      function_name = "compute_element_cm_att_memory_lddrk"
    else
      raise "Unsupported type : #{type}!"
    end

    v = []
    v.push tx                = Int( "tx",                :dir => :in)
    v.push working_element   = Int( "working_element",   :dir => :in)
    v.push d_muvstore        = Real("d_muvstore",        :dir => :in, :dim => [Dim()] )
    v.push factor_common     = Real("factor_common",     :dir => :in, :dim => [Dim()] )
    case type
    when :inner_core, :crust_mantle
      v.push alphaval          = Real("alphaval",          :dir => :in, :dim => [Dim()] )
      v.push betaval           = Real("betaval",           :dir => :in, :dim => [Dim()] )
      v.push gammaval          = Real("gammaval",          :dir => :in, :dim => [Dim()] )
    when :inner_core_lddrk, :crust_mantle_lddrk
      v.push tau_sigmainvval   = Real("tau_sigmainvval",   :dir => :in, :dim => [Dim()] )
    end
    v.push r_xx              = Real("R_xx",              :dir => :inout, :dim => [Dim()] )
    v.push r_yy              = Real("R_yy",              :dir => :inout, :dim => [Dim()] )
    v.push r_xy              = Real("R_xy",              :dir => :inout, :dim => [Dim()] )
    v.push r_xz              = Real("R_xz",              :dir => :inout, :dim => [Dim()] )
    v.push r_yz              = Real("R_yz",              :dir => :inout, :dim => [Dim()] )
    case type
    when :inner_core, :crust_mantle
      v.push epsilondev_xx     = Real("epsilondev_xx",     :dir => :in, :dim => [Dim()] )
      v.push epsilondev_yy     = Real("epsilondev_yy",     :dir => :in, :dim => [Dim()] )
      v.push epsilondev_xy     = Real("epsilondev_xy",     :dir => :in, :dim => [Dim()] )
      v.push epsilondev_xz     = Real("epsilondev_xz",     :dir => :in, :dim => [Dim()] )
      v.push epsilondev_yz     = Real("epsilondev_yz",     :dir => :in, :dim => [Dim()] )
    when :inner_core_lddrk, :crust_mantle_lddrk
      v.push r_xx_lddrk      = Real("R_xx_lddrk",        :dir => :inout, :dim => [Dim()] )
      v.push r_yy_lddrk      = Real("R_yy_lddrk",        :dir => :inout, :dim => [Dim()] )
      v.push r_xy_lddrk      = Real("R_xy_lddrk",        :dir => :inout, :dim => [Dim()] )
      v.push r_xz_lddrk      = Real("R_xz_lddrk",        :dir => :inout, :dim => [Dim()] )
      v.push r_yz_lddrk      = Real("R_yz_lddrk",        :dir => :inout, :dim => [Dim()] )
      v.push alpha_lddrk     = Real("alpha_lddrk",       :dir => :in )
      v.push beta_lddrk      = Real("beta_lddrk",        :dir => :in )
      v.push deltat          = Real("deltat",            :dir => :in )
    end
    v.push epsilondev_xx_loc = Real("epsilondev_xx_loc", :dir => :in)
    v.push epsilondev_yy_loc = Real("epsilondev_yy_loc", :dir => :in)
    v.push epsilondev_xy_loc = Real("epsilondev_xy_loc", :dir => :in)
    v.push epsilondev_xz_loc = Real("epsilondev_xz_loc", :dir => :in)
    v.push epsilondev_yz_loc = Real("epsilondev_yz_loc", :dir => :in)
    v.push use_3d_attenuation_arrays = Int( "USE_3D_ATTENUATION_ARRAYS",    :dir => :in)

    ngll3 = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)
    nsls  = Int("N_SLS", :const => n_sls)

    p = Procedure(function_name, v, :local => true) {
      decl offset = Int("offset")
      decl i_sls  = Int("i_sls")
      decl mul = Real("mul")
      decl factor_loc = Real("factor_loc")
      decl sn = Real("sn")
      decl snp1 = Real("snp1")
      case type
      when :inner_core, :crust_mantle
        decl alphaval_loc = Real("alphaval_loc")
        decl betaval_loc = Real("betaval_loc")
        decl gammaval_loc = Real("gammaval_loc")
      when :inner_core_lddrk, :crust_mantle_lddrk
        decl tau_sigmainv_loc = Real("tau_sigmainv_loc")
      end
      comment()

      print mul === d_muvstore[tx + ngll3_padded*working_element]
      comment()

      print For( i_sls, 0, nsls - 1 ) {
        # indices
        # note: index for R_xx,... here is (i,j,k,i_sls,ispec) and not (i,j,k,ispec,i_sls) as in local version
        #
        # index:
        # (i,j,k,i_sls,ispec) -> offset_sls = tx + NGLL3*(i_sls + N_SLS*working_element)
        print offset === tx + ngll3*(i_sls + nsls*working_element)
        print If(use_3d_attenuation_arrays => lambda {
            print factor_loc  === mul * factor_common[offset]
        }, :else => lambda {
            print factor_loc  === mul * factor_common[i_sls + nsls*working_element ]
        })

        case type
        when :inner_core, :crust_mantle
          print alphaval_loc === alphaval[i_sls]
          print  betaval_loc ===  betaval[i_sls]
          print gammaval_loc === gammaval[i_sls]
          comment()

          # non-LDDRK update
          [[r_xx, epsilondev_xx, epsilondev_xx_loc],
           [r_yy, epsilondev_yy, epsilondev_yy_loc],
           [r_xy, epsilondev_xy, epsilondev_xy_loc],
           [r_xz, epsilondev_xz, epsilondev_xz_loc],
           [r_yz, epsilondev_yz, epsilondev_yz_loc]].each { |r, epsilondev, epsilondev_loc|
            print sn   === factor_loc * epsilondev[tx + ngll3 * working_element]
            print snp1 === factor_loc * epsilondev_loc
            print r[offset] === alphaval_loc*r[offset] + betaval_loc*sn + gammaval_loc*snp1
          }
        when :inner_core_lddrk, :crust_mantle_lddrk
          # LDDRK update
          print tau_sigmainv_loc === tau_sigmainvval[i_sls]
          comment()

          [[r_xx, r_xx_lddrk, epsilondev_xx_loc],
           [r_yy, r_yy_lddrk, epsilondev_yy_loc],
           [r_xy, r_xy_lddrk, epsilondev_xy_loc],
           [r_xz, r_xz_lddrk, epsilondev_xz_loc],
           [r_yz, r_yz_lddrk, epsilondev_yz_loc]].each { |r, r_lddrk, epsilondev_loc|
            # see compute_element_att_memory_ic_lddrk
            #R_xx_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_xx_lddrk(INDEX_IJK,i_SLS,ispec) &
            #  + deltat * ( factor_common_use(INDEX_IJK)*epsilondev_loc(INDEX_IJK,1) &
            #               - R_xx(INDEX_IJK,i_SLS,ispec)*tau_sigmainv_CUSTOM_REAL(i_SLS) )
            #
            #R_xx(INDEX_IJK,i_SLS,ispec) = R_xx(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_xx_lddrk(INDEX_IJK,i_SLS,ispec)

            print sn   === tau_sigmainv_loc * r[offset]
            print snp1 === factor_loc * epsilondev_loc

            print r_lddrk[offset] === alpha_lddrk * r_lddrk[offset] + deltat * (snp1 - sn)
            print r[offset] === r[offset] + beta_lddrk * r_lddrk[offset]
          }
        end
      }
    }
    return p
  end
end
