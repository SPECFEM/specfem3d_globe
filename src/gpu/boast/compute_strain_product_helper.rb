module BOAST

  def BOAST::compute_strain_product
    function_name = "compute_strain_product"

    prod               = Real("prod",              :dir => :out,:dim => [Dim(21)], :register => true )
    eps_trace_over_3   = Real("eps_trace_over_3",  :dir => :in)
    epsdev             = Real("epsdev",            :dir => :in, :dim => [Dim(5)], :register => true )
    b_eps_trace_over_3 = Real("b_eps_trace_over_3",:dir => :in)
    b_epsdev           = Real("b_epsdev",          :dir => :in, :dim => [Dim(5)], :register => true )

    sub = Procedure(function_name, [prod, eps_trace_over_3, epsdev, b_eps_trace_over_3, b_epsdev], :local => true) {
      decl eps   = Real("eps", :dim => [Dim(6)], :allocate => true )
      decl b_eps = Real("b_eps", :dim => [Dim(6)], :allocate => true )
      comment()

      print eps[0] === epsdev[0] + eps_trace_over_3
      print eps[1] === epsdev[1] + eps_trace_over_3
      print eps[2] === -(eps[0] + eps[1]) + eps_trace_over_3*3.0
      print eps[3] === epsdev[4]
      print eps[4] === epsdev[3]
      print eps[5] === epsdev[2]
      comment()

      print b_eps[0] === b_epsdev[0] + b_eps_trace_over_3
      print b_eps[1] === b_epsdev[1] + b_eps_trace_over_3
      print b_eps[2] === -(b_eps[0] + b_eps[1]) + b_eps_trace_over_3*3.0
      print b_eps[3] === b_epsdev[4]
      print b_eps[4] === b_epsdev[3]
      print b_eps[5] === b_epsdev[2]
      comment()

      # Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
      print prod[0] === (eps[0]) * (b_eps[0]);
      print prod[1] === (eps[0]) * (b_eps[1]);
      print prod[1] === prod[1] + (eps[1]) * (b_eps[0]);
      print prod[2] === (eps[0]) * (b_eps[2]);
      print prod[2] === prod[2] + (eps[2]) * (b_eps[0]);
      print prod[3] === (eps[0]) * (b_eps[3]);
      print prod[3] === prod[3] + (eps[3]) * (b_eps[0]);
      print prod[3] === (prod[3]) * (2.0);
      print prod[4] === (eps[0]) * (b_eps[4]);
      print prod[4] === prod[4] + (eps[4]) * (b_eps[0]);
      print prod[4] === (prod[4]) * (2.0);
      print prod[5] === (eps[0]) * (b_eps[5]);
      print prod[5] === prod[5] + (eps[5]) * (b_eps[0]);
      print prod[5] === (prod[5]) * (2.0);
      print prod[6] === (eps[1]) * (b_eps[1]);
      print prod[7] === (eps[1]) * (b_eps[2]);
      print prod[7] === prod[7] + (eps[2]) * (b_eps[1]);
      print prod[8] === (eps[1]) * (b_eps[3]);
      print prod[8] === prod[8] + (eps[3]) * (b_eps[1]);
      print prod[8] === (prod[8]) * (2.0);
      print prod[9] === (eps[1]) * (b_eps[4]);
      print prod[9] === prod[9] + (eps[4]) * (b_eps[1]);
      print prod[9] === (prod[9]) * (2.0);
      print prod[10] === (eps[1]) * (b_eps[5]);
      print prod[10] === prod[10] + (eps[5]) * (b_eps[1]);
      print prod[10] === (prod[10]) * (2.0);
      print prod[11] === (eps[2]) * (b_eps[2]);
      print prod[12] === (eps[2]) * (b_eps[3]);
      print prod[12] === prod[12] + (eps[3]) * (b_eps[2]);
      print prod[12] === (prod[12]) * (2.0);
      print prod[13] === (eps[2]) * (b_eps[4]);
      print prod[13] === prod[13] + (eps[4]) * (b_eps[2]);
      print prod[13] === (prod[13]) * (2.0);
      print prod[14] === (eps[2]) * (b_eps[5]);
      print prod[14] === prod[14] + (eps[5]) * (b_eps[2]);
      print prod[14] === (prod[14]) * (2.0);
      print prod[15] === (eps[3]) * (b_eps[3]);
      print prod[15] === (prod[15]) * (4.0);
      print prod[16] === (eps[3]) * (b_eps[4]);
      print prod[16] === prod[16] + (eps[4]) * (b_eps[3]);
      print prod[16] === (prod[16]) * (4.0);
      print prod[17] === (eps[3]) * (b_eps[5]);
      print prod[17] === prod[17] + (eps[5]) * (b_eps[3]);
      print prod[17] === (prod[17]) * (4.0);
      print prod[18] === (eps[4]) * (b_eps[4]);
      print prod[18] === (prod[18]) * (4.0);
      print prod[19] === (eps[4]) * (b_eps[5]);
      print prod[19] === prod[19] + (eps[5]) * (b_eps[4]);
      print prod[19] === (prod[19]) * (4.0);
      print prod[20] === (eps[5]) * (b_eps[5]);
      print prod[20] === (prod[20]) * (4.0);

    }
    return sub
  end

end
