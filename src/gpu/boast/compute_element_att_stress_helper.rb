module BOAST
  def BOAST::compute_element_att_stress( type, n_gll3 = 125, n_sls = 3, use_cuda_shared_sync = false )
    if type == :inner_core
      function_name = "compute_element_ic_att_stress"
    elsif type == :crust_mantle
      function_name = "compute_element_cm_att_stress"
    else
      raise "Unsupported interface type : #{type}!"
    end
    v = []
    v.push tx              = Int( "tx",              :dir => :in)
    v.push working_element = Int( "working_element", :dir => :in)
    v.push r_xx            = Real("R_xx",            :dir => :in, :dim => [Dim()] )
    v.push r_yy            = Real("R_yy",            :dir => :in, :dim => [Dim()] )
    v.push r_xy            = Real("R_xy",            :dir => :in, :dim => [Dim()] )
    v.push r_xz            = Real("R_xz",            :dir => :in, :dim => [Dim()] )
    v.push r_yz            = Real("R_yz",            :dir => :in, :dim => [Dim()] )
    v.push sigma_xx        = Real("sigma_xx",        :dir => :inout, :dim => [Dim()], :private => true )
    v.push sigma_yy        = Real("sigma_yy",        :dir => :inout, :dim => [Dim()], :private => true )
    v.push sigma_zz        = Real("sigma_zz",        :dir => :inout, :dim => [Dim()], :private => true )
    v.push sigma_xy        = Real("sigma_xy",        :dir => :inout, :dim => [Dim()], :private => true )
    v.push sigma_xz        = Real("sigma_xz",        :dir => :inout, :dim => [Dim()], :private => true )
    v.push sigma_yz        = Real("sigma_yz",        :dir => :inout, :dim => [Dim()], :private => true )

    ngll3 = Int("NGLL3", :const => n_gll3)
    nsls  = Int("N_SLS", :const => n_sls)

    # cuda asynchronuous copies modification
    # to allow the two expressions as function calls below
    # not used yet, as it leads to problems with the output in subsequent kernel files
    #if use_cuda_shared_sync then
    #  register_funccall("compute_offset_sh")
    #  register_funccall("compute_offset")
    #end

    p = Procedure(function_name, v, :local => true) {
      decl offset = Int("offset")
      decl r_xx_val = Real("R_xx_val")
      decl r_yy_val = Real("R_yy_val")
      comment()

      i_sls  = Int("i_sls")
      print For( i_sls, 0, nsls, :operator => "<", :declit => true) {
        # index
        # note: index for R_xx,.. here is (i,j,k,i_sls,ispec) and not (i,j,k,ispec,i_sls) as in local version
        #       see local version: offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);
        # indexing examples:
        #   (i,j,k,ispec,i_sls) -> offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls)
        #   (i_sls,i,j,k,ispec) -> offset_sls = i_sls + N_SLS*(tx + NGLL3*working_element)
        #   (i,j,k,i_sls,ispec) -> offset_sls = tx + NGLL3*(i_sls + N_SLS*working_element)
        comment()
        # cuda asynchronuous copies modification
        if use_cuda_shared_sync then
          get_output.puts "#ifdef CUDA_SHARED_ASYNC"
          #print offset === compute_offset_sh(tx, i_sls)
          get_output.puts "    offset = compute_offset_sh(tx, i_sls);"
          get_output.puts "#else"
          #print offset === compute_offset(tx, i_sls, working_element)
          #get_output.puts "    offset = compute_offset(tx, i_sls, working_element);"
          print offset === tx + ngll3*(i_sls + nsls*working_element)
          get_output.puts "#endif"
        else
          print offset === tx + ngll3*(i_sls + nsls*working_element)
        end
        comment()
        print r_xx_val === r_xx[offset]
        print r_yy_val === r_yy[offset]
        comment()

        print sigma_xx[0] === sigma_xx[0] - r_xx_val
        print sigma_yy[0] === sigma_yy[0] - r_yy_val
        print sigma_zz[0] === sigma_zz[0] + r_xx_val + r_yy_val
        print sigma_xy[0] === sigma_xy[0] - r_xy[offset]
        print sigma_xz[0] === sigma_xz[0] - r_xz[offset]
        print sigma_yz[0] === sigma_yz[0] - r_yz[offset]
      }
    }
    return p
  end
end
