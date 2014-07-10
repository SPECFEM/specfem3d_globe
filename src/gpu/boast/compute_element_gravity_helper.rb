module BOAST
  def BOAST::compute_element_gravity( type, n_gll3 = 125, r_earth_km = 6371.0 )
    if type == :inner_core then
      function_name = "compute_element_ic_gravity"
    elsif type == :crust_mantle then
      function_name = "compute_element_cm_gravity"
    else
      raise "Unsupported type : #{type}!"
    end
    ngll3 = Int("NGLL3", :const => n_gll3)
    v = []
    v.push tx                = Int( "tx",                :dir => :in)
    v.push iglob             = Int( "iglob",             :dir => :in)
    #v.push working_element   = Int( "working_element",   :dir => :in)
    #v.push d_ibool                 = Int("d_ibool",                  :dir => :in, :dim => [Dim()] )
    v.push *d_store = [d_xstore    = Real("d_xstore",                :dir => :in, :restrict => true, :dim => [Dim()] ), d_ystore = Real("d_ystore",:dir => :in, :restrict => true, :dim => [Dim()] ), d_zstore = Real("d_zstore",:dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push d_minus_gravity_table   = Real("d_minus_gravity_table",   :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_minus_deriv_gravity_table = Real("d_minus_deriv_gravity_table", :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_density_table         = Real("d_density_table",         :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push wgll_cube               = Real("wgll_cube",               :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push jacobianl               = Real("jacobianl", :dir => :in)
    v.push *s_dummy_loc = ["x", "y", "z"].collect { |a|
                                     Real("s_dummy#{a}_loc", :dir => :in, :dim => [Dim(ngll3)], :local => true )
    }
    sigma = ["x", "y", "z"].collect { |a1|
        ["x", "y", "z"].collect { |a2|
                                     Real("sigma_#{a1}#{a2}", :dir => :inout, :dim => [Dim()], :private => true )
        }
    }
    v.push sigma[0][0], sigma[1][1], sigma[2][2],\
           sigma[0][1], sigma[1][0], sigma[0][2],\
           sigma[2][0], sigma[1][2], sigma[2][1]
    v.push *rho_s_H = [1,2,3].collect {|n|
                                     Real("rho_s_H#{n}", :dir => :inout, :dim => [Dim()], :private => true )
    }

    p = Procedure(function_name, v, [], :local => true) {
      decl radius = Real("radius"), theta = Real("theta"), phi = Real("phi")
      decl cos_theta = Real("cos_theta"), sin_theta = Real("sin_theta"), cos_phi = Real("cos_phi"), sin_phi = Real("sin_phi")
      decl cos_theta_sq = Real("cos_theta_sq"), sin_theta_sq = Real("sin_theta_sq"), cos_phi_sq = Real("cos_phi_sq"), sin_phi_sq = Real("sin_phi_sq")
      decl minus_g = Real("minus_g"), minus_dg = Real("minus_dg")
      decl rho = Real("rho")
      decl *gl = [ Real("gxl"), Real("gyl"), Real("gzl") ]
      decl minus_g_over_radius = Real("minus_g_over_radius"), minus_dg_plus_g_over_radius = Real("minus_dg_plus_g_over_radius")
      decl hxxl = Real("Hxxl"), hyyl = Real("Hyyl"), hzzl = Real("Hzzl"), hxyl = Real("Hxyl"), hxzl = Real("Hxzl"), hyzl = Real("Hyzl")
      decl *s_l = [ Real("sx_l"), Real("sy_l"), Real("sz_l") ]
      decl factor = Real("factor")
      decl int_radius = Int("int_radius")

      print radius === d_xstore[iglob]
      print If(radius < ( 100.0 / (r_earth_km*1000.0))) {
        print radius ===  100.0 / (r_earth_km*1000.0)
      }
      print theta === d_ystore[iglob]
      print phi === d_zstore[iglob]
      if(get_lang == CL) then
        print sin_theta === sincos(theta, cos_theta.address)
        print sin_phi   === sincos(phi,   cos_phi.address)
      else
        if(get_default_real_size == 4) then
          print sincosf(theta, sin_theta.address, cos_theta.address)
          print sincosf(phi,   sin_phi.address,   cos_phi.address)
        else
          print cos_theta === cos(theta)
          print sin_theta === sin(theta)
          print cos_phi   === cos(phi)
          print sin_phi   === sin(phi)
        end
      end
      print int_radius === rint(radius * r_earth_km * 10.0 ) - 1
      print If( int_radius < 0 ) {
        print int_radius === 0
      }
      print minus_g  === d_minus_gravity_table[int_radius]
      print minus_dg === d_minus_deriv_gravity_table[int_radius]
      print rho      === d_density_table[int_radius]

      print gl[0] === minus_g*sin_theta*cos_phi
      print gl[1] === minus_g*sin_theta*sin_phi
      print gl[2] === minus_g*cos_theta

      print minus_g_over_radius === minus_g / radius
      print minus_dg_plus_g_over_radius === minus_dg - minus_g_over_radius

      print cos_theta_sq === cos_theta*cos_theta
      print sin_theta_sq === sin_theta*sin_theta
      print cos_phi_sq   === cos_phi*cos_phi
      print sin_phi_sq   === sin_phi*sin_phi

      print hxxl === minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq
      print hyyl === minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq
      print hzzl === cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq
      print hxyl === cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq
      print hxzl === cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta
      print hyzl === cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta

      (0..2).each { |indx|
        print s_l[indx] === rho * s_dummy_loc[indx][tx]
      }
      print sigma[0][0].dereference === sigma[0][0].dereference + s_l[1]*gl[1] + s_l[2]*gl[2];
      print sigma[1][1].dereference === sigma[1][1].dereference + s_l[0]*gl[0] + s_l[2]*gl[2];
      print sigma[2][2].dereference === sigma[2][2].dereference + s_l[0]*gl[0] + s_l[1]*gl[1];

      print sigma[0][1].dereference === sigma[0][1].dereference - s_l[0] * gl[1];
      print sigma[1][0].dereference === sigma[1][0].dereference - s_l[1] * gl[0];

      print sigma[0][2].dereference === sigma[0][2].dereference - s_l[0] * gl[2];
      print sigma[2][0].dereference === sigma[2][0].dereference - s_l[2] * gl[0];

      print sigma[1][2].dereference === sigma[1][2].dereference - s_l[1] * gl[2];
      print sigma[2][1].dereference === sigma[2][1].dereference - s_l[2] * gl[1];

      print factor === jacobianl * wgll_cube[tx]
      print rho_s_H[0][0] === factor * (s_l[0]*hxxl + s_l[1]*hxyl + s_l[2]*hxzl)
      print rho_s_H[1][0] === factor * (s_l[0]*hxyl + s_l[1]*hyyl + s_l[2]*hyzl)
      print rho_s_H[2][0] === factor * (s_l[0]*hxzl + s_l[1]*hyzl + s_l[2]*hzzl)
    }
    return p
  end
end
