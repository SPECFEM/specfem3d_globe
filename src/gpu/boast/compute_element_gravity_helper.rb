module BOAST
  def BOAST::compute_element_gravity( type, n_gll3 = 125, use_cuda_shared_async = false )
    if type == :inner_core then
      function_name = "compute_element_ic_gravity"
    elsif type == :crust_mantle then
      function_name = "compute_element_cm_gravity"
    else
      raise "Unsupported type : #{type}!"
    end
    ngll3 = Int("NGLL3", :const => n_gll3)
    v = []
    v.push tx                      = Int( "tx",                :dir => :in)
    v.push iglob                   = Int( "iglob",             :dir => :in)
    v.push gravity_pre_store       = Real("gravity_pre_store",     :dir => :in, :restrict => true, :dim => [Dim(3), Dim()] )
    v.push gravity_H               = Real("gravity_H",             :dir => :in, :restrict => true, :dim => [Dim(6), Dim()] )
    v.push wgll_cube               = Real("wgll_cube",             :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push jacobianl               = Real("jacobianl", :dir => :in)
    if use_cuda_shared_async then
      v.push s_dummy_loc             = Real("s_dummy_loc", :dir => :in, :dim => [Dim(3),Dim()], :local => true )
    else
      v.push *s_dummy_loc = ["x", "y", "z"].collect { |a|
                                     Real("s_dummy#{a}_loc", :dir => :in, :dim => [Dim(ngll3)], :local => true )
      }
    end

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

    p = Procedure(function_name, v, :local => true) {
      decl *gl = [ Real("gxl"), Real("gyl"), Real("gzl") ]
      decl hxxl = Real("Hxxl"), hyyl = Real("Hyyl"), hzzl = Real("Hzzl"), hxyl = Real("Hxyl"), hxzl = Real("Hxzl"), hyzl = Real("Hyzl")
      decl *s_l = [ Real("sx_l"), Real("sy_l"), Real("sz_l") ]
      decl factor = Real("factor")
      comment()

      # gravity
      # new with pre-calculated arrays
      # Cartesian components of the gravitational acceleration
      if use_cuda_shared_async then
        print gl[0] === gravity_pre_store[0,tx] # minus_g*sin_theta*cos_phi * rho
        print gl[1] === gravity_pre_store[1,tx] # minus_g*sin_theta*sin_phi * rho
        print gl[2] === gravity_pre_store[2,tx] # minus_g*cos_theta * rho
        comment()
        print hxxl === gravity_H[0,tx]
        print hyyl === gravity_H[1,tx]
        print hzzl === gravity_H[2,tx]
        print hxyl === gravity_H[3,tx]
        print hxzl === gravity_H[4,tx]
        print hyzl === gravity_H[5,tx]
      else
        print gl[0] === gravity_pre_store[0,iglob] # minus_g*sin_theta*cos_phi * rho
        print gl[1] === gravity_pre_store[1,iglob] # minus_g*sin_theta*sin_phi * rho
        print gl[2] === gravity_pre_store[2,iglob] # minus_g*cos_theta * rho
        comment()
        print hxxl === gravity_H[0,iglob]
        print hyyl === gravity_H[1,iglob]
        print hzzl === gravity_H[2,iglob]
        print hxyl === gravity_H[3,iglob]
        print hxzl === gravity_H[4,iglob]
        print hyzl === gravity_H[5,iglob]
      end
      comment()

      if use_cuda_shared_async then
        (0..2).each { |indx|
          print s_l[indx] === s_dummy_loc[indx,tx]
        }
      else
        (0..2).each { |indx|
          print s_l[indx] === s_dummy_loc[indx][tx]
        }
      end
      comment()

      print sigma[0][0].dereference === sigma[0][0].dereference + s_l[1]*gl[1] + s_l[2]*gl[2];
      print sigma[1][1].dereference === sigma[1][1].dereference + s_l[0]*gl[0] + s_l[2]*gl[2];
      print sigma[2][2].dereference === sigma[2][2].dereference + s_l[0]*gl[0] + s_l[1]*gl[1];

      print sigma[0][1].dereference === sigma[0][1].dereference - s_l[0] * gl[1];
      print sigma[1][0].dereference === sigma[1][0].dereference - s_l[1] * gl[0];

      print sigma[0][2].dereference === sigma[0][2].dereference - s_l[0] * gl[2];
      print sigma[2][0].dereference === sigma[2][0].dereference - s_l[2] * gl[0];

      print sigma[1][2].dereference === sigma[1][2].dereference - s_l[1] * gl[2];
      print sigma[2][1].dereference === sigma[2][1].dereference - s_l[2] * gl[1];
      comment()

      # precompute vector
      print factor === jacobianl * wgll_cube[tx]
      print rho_s_H[0].dereference === factor * (s_l[0]*hxxl + s_l[1]*hxyl + s_l[2]*hxzl)
      print rho_s_H[1].dereference === factor * (s_l[0]*hxyl + s_l[1]*hyyl + s_l[2]*hyzl)
      print rho_s_H[2].dereference === factor * (s_l[0]*hxzl + s_l[1]*hyzl + s_l[2]*hzzl)
    }
    return p
  end
end
