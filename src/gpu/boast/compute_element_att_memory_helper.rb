module BOAST
  def BOAST::compute_element_att_memory(type, n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, use_cuda_shared_async = false )
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
    if use_cuda_shared_async then
      v.push sh_mul            = Real("sh_mul",            :dir => :in, :dim => [Dim()] )
      v.push sh_factor_common  = Real("sh_factor_common",  :dir => :in, :dim => [Dim()] )
    else
      v.push d_muvstore        = Real("d_muvstore",        :dir => :in, :dim => [Dim()] )
      v.push factor_common     = Real("factor_common",     :dir => :in, :dim => [Dim()] )
    end
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
    if use_cuda_shared_async then
      v.push sh_r_xx              = Real("sh_R_xx",         :dir => :inout, :dim => [Dim()] )
      v.push sh_r_yy              = Real("sh_R_yy",         :dir => :inout, :dim => [Dim()] )
      v.push sh_r_xy              = Real("sh_R_xy",         :dir => :inout, :dim => [Dim()] )
      v.push sh_r_xz              = Real("sh_R_xz",         :dir => :inout, :dim => [Dim()] )
      v.push sh_r_yz              = Real("sh_R_yz",         :dir => :inout, :dim => [Dim()] )
    end
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

    # cuda asynchronuous copies modification
    # to allow the two expressions as function calls below
    # not used yet, as it leads to problems with the output in subsequent kernel files
    #if use_cuda_shared_async then
    #  register_funccall("compute_offset_sh")
    #  register_funccall("compute_offset")
    #end

    p = Procedure(function_name, v, :local => true) {
      decl offset = Int("offset")
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

      if use_cuda_shared_async then
        # not needed in this routine here, will put #ifdef-statement around whole routine
        # because function argument list changes depending on cuda_shared_async case
        # and BOAST doesn't allow to put #ifdef-statements inside the function arguments list
        # get_output.puts "#ifdef CUDA_SHARED_ASYNC"
        decl offset_sh = Int("offset_sh")
        comment()
        print mul === sh_mul[tx];
      else
        print mul === d_muvstore[tx + ngll3_padded*working_element]
        if (type == :inner_core or type == :crust_mantle) then
          decl offset_eps = Int("offset_eps", :const => tx + ngll3 * working_element)
        end
      end
      comment()

      i_sls  = Int("i_sls")
      print For( i_sls, 0, nsls, :operator => "<", :declit => true ) {
        # indices
        # note: index for R_xx,... here is (i,j,k,i_sls,ispec) and not (i,j,k,ispec,i_sls) as in local version
        #
        # index:
        # (i,j,k,i_sls,ispec) -> offset_sls = tx + NGLL3*(i_sls + N_SLS*working_element)
        if use_cuda_shared_async then
          #print offset_sh === compute_offset_sh(tx, i_sls)
          get_output.puts "    offset_sh = compute_offset_sh(tx, i_sls);"
          print factor_loc === mul * sh_factor_common[offset_sh]
          #print offset === compute_offset(tx, i_sls, working_element)
          get_output.puts "    offset = compute_offset(tx, i_sls, working_element);"
        else
          print offset === tx + ngll3*(i_sls + nsls*working_element)
          print If(use_3d_attenuation_arrays => lambda {
              print factor_loc  === mul * factor_common[offset]
          }, :else => lambda {
              print factor_loc  === mul * factor_common[i_sls + nsls*working_element ]
          })
        end
        comment()

        case type
        when :inner_core, :crust_mantle
          print alphaval_loc === alphaval[i_sls]
          print  betaval_loc ===  betaval[i_sls]
          print gammaval_loc === gammaval[i_sls]
          comment()

          # non-LDDRK update
          if use_cuda_shared_async then
            [[r_xx, sh_r_xx, epsilondev_xx, epsilondev_xx_loc],
             [r_yy, sh_r_yy, epsilondev_yy, epsilondev_yy_loc],
             [r_xy, sh_r_xy, epsilondev_xy, epsilondev_xy_loc],
             [r_xz, sh_r_xz, epsilondev_xz, epsilondev_xz_loc],
             [r_yz, sh_r_yz, epsilondev_yz, epsilondev_yz_loc]].each { |r, sh_r, epsilondev, epsilondev_loc|
              print sn   === factor_loc * epsilondev[tx]
              print snp1 === factor_loc * epsilondev_loc
              print r[offset] === alphaval_loc*sh_r[offset_sh] + betaval_loc*sn + gammaval_loc*snp1
            }
          else
            [[r_xx, epsilondev_xx, epsilondev_xx_loc],
             [r_yy, epsilondev_yy, epsilondev_yy_loc],
             [r_xy, epsilondev_xy, epsilondev_xy_loc],
             [r_xz, epsilondev_xz, epsilondev_xz_loc],
             [r_yz, epsilondev_yz, epsilondev_yz_loc]].each { |r, epsilondev, epsilondev_loc|
              print sn   === factor_loc * epsilondev[offset_eps]
              print snp1 === factor_loc * epsilondev_loc
              print r[offset] === alphaval_loc*r[offset] + betaval_loc*sn + gammaval_loc*snp1
            }
          end
        when :inner_core_lddrk, :crust_mantle_lddrk
          # LDDRK update
          print tau_sigmainv_loc === tau_sigmainvval[i_sls]
          comment()

          if use_cuda_shared_async then
            [[r_xx, sh_r_xx, r_xx_lddrk, epsilondev_xx_loc],
             [r_yy, sh_r_yy, r_yy_lddrk, epsilondev_yy_loc],
             [r_xy, sh_r_xy, r_xy_lddrk, epsilondev_xy_loc],
             [r_xz, sh_r_xz, r_xz_lddrk, epsilondev_xz_loc],
             [r_yz, sh_r_yz, r_yz_lddrk, epsilondev_yz_loc]].each { |r, sh_r, r_lddrk, epsilondev_loc|
              # see compute_element_att_memory_ic_lddrk
              #R_xx_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_xx_lddrk(INDEX_IJK,i_SLS,ispec) &
              #  + deltat * ( factor_common_use(INDEX_IJK)*epsilondev_loc(INDEX_IJK,1) &
              #               - R_xx(INDEX_IJK,i_SLS,ispec)*tau_sigmainv_CUSTOM_REAL(i_SLS) )
              #
              #R_xx(INDEX_IJK,i_SLS,ispec) = R_xx(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_xx_lddrk(INDEX_IJK,i_SLS,ispec)

              print sn   === tau_sigmainv_loc * sh_r[offset_sh]
              print snp1 === factor_loc * epsilondev_loc

              print r_lddrk[offset] === alpha_lddrk * r_lddrk[offset] + deltat * (snp1 - sn)
              print r[offset] === sh_r[offset_sh] + beta_lddrk * r_lddrk[offset]
            }
          else
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
        end
      }
    }
    return p
  end
end
