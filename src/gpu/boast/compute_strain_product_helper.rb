module BOAST

  def BOAST::compute_strain_product
    function_name = "compute_strain_product"

    prod               = Real("prod",              :dir => :out,:dim => [Dim(21)], :register => true )
    eps_trace_over_3   = Real("eps_trace_over_3",  :dir => :in)
    epsdev             = Real("epsdev",            :dir => :in, :dim => [Dim(5)], :register => true )
    b_eps_trace_over_3 = Real("b_eps_trace_over_3",:dir => :in)
    b_epsdev           = Real("b_epsdev",          :dir => :in, :dim => [Dim(5)], :register => true )

    sub = Procedure(function_name, [prod, eps_trace_over_3, epsdev, b_eps_trace_over_3, b_epsdev], [], :local => true) {
      decl eps   = Real("eps", :dim => [Dim(6)], :allocate => true )
      decl b_eps = Real("b_eps", :dim => [Dim(6)], :allocate => true )
      decl p = Int("p")
      decl i = Int("i")
      decl j = Int("j")

      print eps[0] === epsdev[0] + eps_trace_over_3
      print eps[1] === epsdev[1] + eps_trace_over_3
      print eps[2] === -(eps[0] + eps[1]) + eps_trace_over_3*3.0
      print eps[3] === epsdev[4]
      print eps[4] === epsdev[3]
      print eps[5] === epsdev[2]

      print b_eps[0] === b_epsdev[0] + b_eps_trace_over_3
      print b_eps[1] === b_epsdev[1] + b_eps_trace_over_3
      print b_eps[2] === -(b_eps[0] + b_eps[1]) + b_eps_trace_over_3*3.0
      print b_eps[3] === b_epsdev[4]
      print b_eps[4] === b_epsdev[3]
      print b_eps[5] === b_epsdev[2]

      print p === 0
      print For(i, 0, 5) {
        print For(j, 0, 5) {
          print prod[p] === eps[i]*b_eps[j]
          print If(j>i) {
            print prod[p] === prod[p] + eps[j]*b_eps[i]
            print If( Expression("&&", j>2, i<3) ) {
              print prod[p] === prod[p]*2.0
            }
          }
        }
      }
    }
    return sub
  end

end
