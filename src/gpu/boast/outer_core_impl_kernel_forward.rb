module BOAST
  def BOAST::compute_element_oc_rotation( n_gll3 = 125)
    function_name = "compute_element_oc_rotation"
    v = []
    v.push tx                    = Int( "tx",                    :dir => :in)
    v.push working_element       = Int( "working_element",       :dir => :in)
    v.push time                  = Real("time",                  :dir => :in)
    v.push two_omega_earth       = Real("two_omega_earth",       :dir => :in)
    v.push deltat                = Real("deltat",                :dir => :in)
    v.push d_A_array_rotation    = Real("d_A_array_rotation",    :dir => :inout, :dim => [Dim()] )
    v.push d_B_array_rotation    = Real("d_B_array_rotation",    :dir => :inout, :dim => [Dim()] )
    v.push d_A_array_rotation_lddrk = Real("d_A_array_rotation_lddrk",    :dir => :inout, :dim => [Dim()] )
    v.push d_B_array_rotation_lddrk = Real("d_B_array_rotation_lddrk",    :dir => :inout, :dim => [Dim()] )
    v.push alpha_lddrk           = Real("alpha_lddrk",           :dir => :in)
    v.push beta_lddrk            = Real("beta_lddrk",            :dir => :in)
    v.push use_lddrk             = Int( "use_lddrk",             :dir => :in)
    v.push dpotentialdxl         = Real("dpotentialdxl",         :dir => :in)
    v.push dpotentialdyl         = Real("dpotentialdyl",         :dir => :in)
    v.push dpotentialdx_with_rot = Real("dpotentialdx_with_rot", :dir => :out, :dim => [Dim(1)], :private => true)
    v.push dpotentialdy_with_rot = Real("dpotentialdy_with_rot", :dir => :out, :dim => [Dim(1)], :private => true)

    ngll3 = Int("NGLL3", :const => n_gll3)

    p = Procedure(function_name, v, :local => true) {
      decl two_omega_deltat = Real("two_omega_deltat")
      decl cos_two_omega_t  = Real("cos_two_omega_t")
      decl sin_two_omega_t  = Real("sin_two_omega_t")
      decl a_rotation = Real("A_rotation")
      decl b_rotation = Real("B_rotation")
      decl source_euler_A = Real("source_euler_A")
      decl source_euler_B = Real("source_euler_B")
      comment()

      comment("  // store the source for the Euler scheme for A_rotation and B_rotation")
      if (get_lang == CL) then
        print sin_two_omega_t === sincos(two_omega_earth*time, cos_two_omega_t.address)
      else
        if (get_default_real_size == 4) then
          print sincosf(two_omega_earth*time, sin_two_omega_t.address, cos_two_omega_t.address)
        else
          print cos_two_omega_t === cos(two_omega_earth*time)
          print sin_two_omega_t === sin(two_omega_earth*time)
        end
      end
      comment()

      comment("  // time step deltat of Euler scheme is included in the source")
      print two_omega_deltat === deltat * two_omega_earth

      print source_euler_A === two_omega_deltat * (cos_two_omega_t * dpotentialdyl + sin_two_omega_t * dpotentialdxl)
      print source_euler_B === two_omega_deltat * (sin_two_omega_t * dpotentialdyl - cos_two_omega_t * dpotentialdxl)

      print a_rotation === d_A_array_rotation[tx + working_element*ngll3]
      print b_rotation === d_B_array_rotation[tx + working_element*ngll3]

      print dpotentialdx_with_rot[0] === dpotentialdxl + ( a_rotation*cos_two_omega_t + b_rotation*sin_two_omega_t)
      print dpotentialdy_with_rot[0] === dpotentialdyl + (-a_rotation*sin_two_omega_t + b_rotation*cos_two_omega_t)
      comment()

      comment("  // updates rotation term with Euler scheme (non-padded offset)")
      print If(!use_lddrk => lambda {
        # non-LDDRK update
        print d_A_array_rotation[tx + working_element*ngll3] === d_A_array_rotation[tx + working_element*ngll3] + source_euler_A
        print d_B_array_rotation[tx + working_element*ngll3] === d_B_array_rotation[tx + working_element*ngll3] + source_euler_B
      }, :else => lambda {
        # LDDRK update
        print d_A_array_rotation_lddrk[tx + working_element*ngll3] === \
                alpha_lddrk * d_A_array_rotation_lddrk[tx + working_element*ngll3] + source_euler_A
        print d_B_array_rotation_lddrk[tx + working_element*ngll3] === \
                alpha_lddrk * d_B_array_rotation_lddrk[tx + working_element*ngll3] + source_euler_B
        print d_A_array_rotation[tx + working_element*ngll3] === \
                d_A_array_rotation[tx + working_element*ngll3] + beta_lddrk * d_A_array_rotation_lddrk[tx + working_element*ngll3]
        print d_B_array_rotation[tx + working_element*ngll3] === \
                d_B_array_rotation[tx + working_element*ngll3] + beta_lddrk * d_B_array_rotation_lddrk[tx + working_element*ngll3]
      })
    }
    return p
  end

  def BOAST::outer_core_impl_kernel_forward(ref = true, elem_per_thread = 1, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = false, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, coloring_min_nspec_outer_core = 1000)
    return BOAST::outer_core_impl_kernel(true, ref, elem_per_thread, mesh_coloring, textures_fields, textures_constants, unroll_loops, n_gllx, n_gll2, n_gll3, n_gll3_padded, coloring_min_nspec_outer_core)
  end

  def BOAST::outer_core_impl_kernel(forward, ref = true, elem_per_thread = 1, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = false, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, coloring_min_nspec_outer_core = 1000)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    v = []
    function_name = "outer_core_impl_kernel"
    if forward then
      function_name += "_forward"
    else
      function_name += "_adjoint"
    end
    v.push nb_blocks_to_compute    = Int("nb_blocks_to_compute",     :dir => :in)
    v.push d_ibool                 = Int("d_ibool",                  :dir => :in, :dim => [Dim()] )
    v.push d_phase_ispec_inner     = Int("d_phase_ispec_inner",      :dir => :in, :dim => [Dim()] )
    v.push num_phase_ispec         = Int("num_phase_ispec",          :dir => :in)
    v.push d_iphase                = Int("d_iphase",                 :dir => :in)
    v.push use_mesh_coloring_gpu   = Int("use_mesh_coloring_gpu",    :dir => :in)
    v.push d_potential             = Real("d_potential",             :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_potential_dot_dot     = Real("d_potential_dot_dot",     :dir => :inout, :dim => [Dim()] )
    v.push *d_xi = [d_xix          = Real("d_xix",                   :dir => :in, :restrict => true, :dim => [Dim()] ),
                    d_xiy          = Real("d_xiy",                   :dir => :in, :restrict => true, :dim => [Dim()] ),
                    d_xiz          = Real("d_xiz",                   :dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push *d_eta = [d_etax        = Real("d_etax",                  :dir => :in, :restrict => true, :dim => [Dim()] ),
                     d_etay        = Real("d_etay",                  :dir => :in, :restrict => true, :dim => [Dim()] ),
                     d_etaz        = Real("d_etaz",                  :dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax    = Real("d_gammax",                :dir => :in, :restrict => true, :dim => [Dim()] ),
                       d_gammay    = Real("d_gammay",                :dir => :in, :restrict => true, :dim => [Dim()] ),
                       d_gammaz    = Real("d_gammaz",                :dir => :in, :restrict => true, :dim => [Dim()] ) ]
    v.push d_hprime_xx             = Real("d_hprime_xx",             :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push d_hprimewgll_xx         = Real("d_hprimewgll_xx",         :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push wgllwgll_xy             = Real("wgllwgll_xy",             :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push wgllwgll_xz             = Real("wgllwgll_xz",             :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push wgllwgll_yz             = Real("wgllwgll_yz",             :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push gravity                 = Int( "GRAVITY",                 :dir => :in)
    v.push gravity_pre_store_outer_core = Real("d_gravity_pre_store_outer_core",:dir => :in, :restrict => true, :dim => [Dim(3), Dim()] )
    v.push wgll_cube               = Real("wgll_cube",               :dir => :in, :restrict => true, :dim => [Dim()] )
    v.push rotation                = Int( "ROTATION",                :dir => :in)
    v.push time                    = Real("time",                    :dir => :in)
    v.push two_omega_earth         = Real("two_omega_earth",         :dir => :in)
    v.push deltat                  = Real("deltat",                  :dir => :in)
    v.push d_A_array_rotation      = Real("d_A_array_rotation",      :dir => :inout, :dim => [Dim()] )
    v.push d_B_array_rotation      = Real("d_B_array_rotation",      :dir => :inout, :dim => [Dim()] )
    v.push d_A_array_rotation_lddrk = Real("d_A_array_rotation_lddrk",      :dir => :inout, :dim => [Dim()] )
    v.push d_B_array_rotation_lddrk = Real("d_B_array_rotation_lddrk",      :dir => :inout, :dim => [Dim()] )
    v.push alpha_lddrk             = Real("alpha_lddrk",             :dir => :in)
    v.push beta_lddrk              = Real("beta_lddrk",              :dir => :in)
    v.push use_lddrk               = Int( "use_lddrk",               :dir => :in)
    v.push nspec_outer_core        = Int( "NSPEC_OUTER_CORE",        :dir => :in)

    ngllx        = Int("NGLLX", :const => n_gllx)
    ngll2        = Int("NGLL2", :const => n_gll2)
    ngll3        = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    use_mesh_coloring       = Int("USE_MESH_COLORING_GPU",   :const => mesh_coloring)
    use_textures_constants  = Int("USE_TEXTURES_CONSTANTS",  :const => textures_constants)
    use_textures_fields     = Int("USE_TEXTURES_FIELDS",     :const => textures_fields)
    manually_unrolled_loops = Int("MANUALLY_UNROLLED_LOOPS", :const => unroll_loops)

    constants = [] #ngllx, ngll2, ngll3, ngll3_padded]
    textures_fields = []
    textures_constants = []

    d_displ_oc_tex = Real("d_#{forward ? "":"b_"}displ_oc_tex", :texture => true, :dir => :in, :dim => [Dim()] )
    d_accel_oc_tex = Real("d_#{forward ? "":"b_"}accel_oc_tex", :texture => true, :dir => :in, :dim => [Dim()] )

    if get_lang == CL then
      textures_fields.push(d_displ_oc_tex, d_accel_oc_tex)
    end
    d_hprime_xx_oc_tex = Real("d_hprime_xx_oc_tex", :texture => true, :dir => :in, :dim => [Dim()] )
    d_hprimewgll_xx_oc_tex = Real("d_hprimewgll_xx_oc_tex", :texture => true, :dir => :in, :dim => [Dim()] )
    if get_lang == CL then
      textures_constants.push(d_hprime_xx_oc_tex, d_hprimewgll_xx_oc_tex)
    end

    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu".gsub("_forward","").gsub("_adjoint",""))
    elsif(get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      # header
      make_specfem3d_header(:ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded,
                            :coloring_min_nspec_outer_core => coloring_min_nspec_outer_core)

      #DEACTIVATE USE TEXTURES CONSTANTS
      get_output.puts "#ifdef #{use_textures_constants}"
      get_output.puts "#undef #{use_textures_constants}"
      get_output.puts "#endif"

#      if get_lang == CUDA
#        get_output.puts "#ifdef #{use_textures_fields}"
#          decl d_displ_oc_tex
#          decl d_accel_oc_tex
#        get_output.puts "#endif"
#        get_output.puts "#ifdef #{use_textures_constants}"
#          decl d_hprime_xx_oc_tex
#        get_output.puts "#endif"
#      end

      comment()
      comment("// fluid rotation")
      comment()

      sub_kernel =  compute_element_oc_rotation(n_gll3)
      print sub_kernel

      comment()
      comment("// KERNEL 2")
      comment("//")
      comment("// for outer core ( acoustic domain )")
      comment()

      if get_lang == CL then
        get_output.puts "#ifdef #{use_textures_fields}"
          get_output.puts "#ifdef #{use_textures_constants}"
            p = Procedure(function_name, v+textures_fields+textures_constants, :constants => constants)
            open p
            set_indent_level(0)
          get_output.puts "#else"
            p = Procedure(function_name, v+textures_fields, :constants => constants)
            open p
            set_indent_level(0)
          get_output.puts "#endif"
        get_output.puts "#else"
          get_output.puts "#ifdef #{use_textures_constants}"
            p = Procedure(function_name, v+textures_constants, :constants => constants)
            open p
            set_indent_level(0)
          get_output.puts "#else"
            p = Procedure(function_name, v, :constants => constants)
            open p
          get_output.puts "#endif"
        get_output.puts "#endif"

        get_output.puts "#ifdef #{use_textures_fields}"
          decl d_displ_oc_tex.sampler
          decl d_accel_oc_tex.sampler
        get_output.puts "#endif"
        get_output.puts "#ifdef #{use_textures_constants}"
          decl d_hprime_xx_oc_tex.sampler
          decl d_hprimewgll_xx_oc_tex.sampler
        get_output.puts "#endif"
      else
        p = Procedure(function_name, v, :constants => constants)
        open p
      end
      decl bx = Int("bx")
      decl tx = Int("tx")
      decl k  = Int("K"), j = Int("J"), i = Int("I")

      #get_output.puts "#ifndef #{manually_unrolled_loops}"
      l = Int("l")
      #get_output.puts "#endif"

      active = (1..elem_per_thread).collect { |e_i| Int("active_#{e_i}", :size => 2, :signed => false) }
      decl *active

      decl offset = Int("offset")
      iglob = (1..elem_per_thread).collect { |e_i| Int("iglob_#{e_i}") }

      decl *iglob
      decl working_element = Int("working_element")
      decl *templ = [ Real("temp1l"), Real("temp2l"), Real("temp3l") ]
      decl *xil   = [ Real("xixl"),   Real("xiyl"),   Real("xizl")   ]
      decl *etal  = [ Real("etaxl"),  Real("etayl"),  Real("etazl")  ]
      decl *gammal= [ Real("gammaxl"),Real("gammayl"),Real("gammazl")]
      decl jacobianl= Real("jacobianl")
      decl *dpotentialdl = [ Real("dpotentialdxl"), Real("dpotentialdyl"), Real("dpotentialdzl") ]
      decl dpotentialdx_with_rot = Real("dpotentialdx_with_rot")
      decl dpotentialdy_with_rot = Real("dpotentialdy_with_rot")

      decl sum_terms = Real("sum_terms")
      gravity_term = (1..elem_per_thread).collect { |e_i| Real("gravity_term_#{e_i}") }
      decl *gravity_term

      # gravity
      decl vec_x = Real("vec_x"), vec_y = Real("vec_y"), vec_z = Real("vec_z")

      decl s_dummy_loc = Real("s_dummy_loc", :local => true, :dim => [Dim(ngll3)] )

      decl *s_temp = [ Real("s_temp1", :local => true, :dim => [Dim(ngll3)]),
                       Real("s_temp2", :local => true, :dim => [Dim(ngll3)]),
                       Real("s_temp3", :local => true, :dim => [Dim(ngll3)])]

      decl sh_hprime_xx     = Real("sh_hprime_xx",     :local => true, :dim => [Dim(ngll2)] )
      decl sh_hprimewgll_xx = Real("sh_hprimewgll_xx", :local => true, :dim => [Dim(ngll2)] )
      comment()

      print bx === get_group_id(1)*get_num_groups(0)+get_group_id(0)

      elem_per_thread.times { |elem_index|
        print tx === get_local_id(0) + ngll3_padded * elem_index / elem_per_thread

        comment()
        comment("  // use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,")
        comment("  // because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses")

        print active[elem_index] === Ternary( Expression("&&", tx < ngll3, bx < nb_blocks_to_compute), 1, 0)

        comment()
        comment("  // copy from global memory to shared memory")
        comment("  // each thread writes one of the NGLL^3 = 125 data points")

        print If(active[elem_index]) {

        if elem_index == 0 then
          get_output.puts "#ifdef #{use_mesh_coloring}"
            print working_element === bx
          get_output.puts "#else"
            print If(use_mesh_coloring_gpu => lambda {
              print working_element === bx
            }, :else => lambda {
              print working_element === d_phase_ispec_inner[bx + num_phase_ispec*(d_iphase-1)]-1
            })
          get_output.puts "#endif"
        end
          print iglob[elem_index] === d_ibool[working_element*ngll3 + tx]-1

          get_output.puts "#ifdef #{use_textures_fields}"
            print s_dummy_loc[tx] === d_displ_oc_tex[iglob[elem_index]]
          get_output.puts "#else"
            print s_dummy_loc[tx] === d_potential[iglob[elem_index]]
          get_output.puts "#endif"
        }
        print If(tx < ngll2) {
          get_output.puts "#ifdef #{use_textures_constants}"
            print sh_hprime_xx[tx] === d_hprime_xx_oc_tex[tx]
            print sh_hprimewgll_xx[tx] === d_hprimewgll_xx_oc_tex[tx]
          get_output.puts "#else"
            print sh_hprime_xx[tx] === d_hprime_xx[tx]
            print sh_hprimewgll_xx[tx] === d_hprimewgll_xx[tx]
          get_output.puts "#endif"

        }
      }
      print barrier(:local)
      comment()


      elem_per_thread.times { |elem_index|
        if elem_per_thread > 1 then
          print tx === get_local_id(0) + ngll3_padded * elem_index / elem_per_thread
        end
        print k === tx/ngll2
        print j === (tx-k*ngll2)/ngllx
        print i === tx - k*ngll2 - j*ngllx
        comment()

        print If(active[elem_index]) {
          (0..2).each { |indx| print templ[indx] === 0.0 }
          for_loop = For(l, 0, ngllx-1, :declit => true) {
             print templ[0] === templ[0] + s_dummy_loc[k*ngll2+j*ngllx+l]*sh_hprime_xx[l*ngllx+i]
             print templ[1] === templ[1] + s_dummy_loc[k*ngll2+l*ngllx+i]*sh_hprime_xx[l*ngllx+j]
             print templ[2] === templ[2] + s_dummy_loc[l*ngll2+j*ngllx+i]*sh_hprime_xx[l*ngllx+k]
          }
          get_output.puts "#ifdef #{manually_unrolled_loops}"
            print for_loop.unroll
          get_output.puts "#else"
            print for_loop
          get_output.puts "#endif"
          comment()

          print offset === working_element*ngll3_padded + tx
          (0..2).each { |indx|
            print xil[indx]    === d_xi[indx][offset]
            print etal[indx]   === d_eta[indx][offset]
            print gammal[indx] === d_gamma[indx][offset]
          }
          comment()

          print jacobianl === Expression("/", 1.0, xil[0]*(etal[1]*gammal[2] - etal[2]*gammal[1]) - xil[1]*(etal[0]*gammal[2] - etal[2]*gammal[0]) + xil[2]*(etal[0]*gammal[1] - etal[1]*gammal[0]))
          (0..2).each { |indx|
            print dpotentialdl[indx] === xil[indx]*templ[0] + etal[indx]*templ[1] + gammal[indx]*templ[2]
          }
          comment()

          print If(rotation => lambda {
            print sub_kernel.call(tx,working_element,\
                                  time,two_omega_earth,deltat,\
                                  d_A_array_rotation,\
                                  d_B_array_rotation,\
                                  d_A_array_rotation_lddrk,\
                                  d_B_array_rotation_lddrk,\
                                  alpha_lddrk,beta_lddrk,use_lddrk, \
                                  dpotentialdl[0],\
                                  dpotentialdl[1],\
                                  dpotentialdx_with_rot.address,\
                                  dpotentialdy_with_rot.address)
          }, :else => lambda {
            print dpotentialdx_with_rot === dpotentialdl[0]
            print dpotentialdy_with_rot === dpotentialdl[1]
          })
          comment()

          # gravity way using pre-calculated arrays
          print vec_x === gravity_pre_store_outer_core[0,iglob[elem_index]]
          print vec_y === gravity_pre_store_outer_core[1,iglob[elem_index]]
          print vec_z === gravity_pre_store_outer_core[2,iglob[elem_index]]

          comment()

          print If(!gravity => lambda {
            print dpotentialdx_with_rot === dpotentialdx_with_rot + s_dummy_loc[tx] * vec_x
            print dpotentialdy_with_rot === dpotentialdy_with_rot + s_dummy_loc[tx] * vec_y
            print dpotentialdl[2]       === dpotentialdl[2] +       s_dummy_loc[tx] * vec_z
          }, :else => lambda {
            print gravity_term[elem_index] === jacobianl * wgll_cube[tx] * \
                                   (dpotentialdx_with_rot * vec_x + dpotentialdy_with_rot * vec_y + dpotentialdl[2] * vec_z)
          })
          comment()

          print s_temp[0][tx] === jacobianl*(   xil[0]*dpotentialdx_with_rot +    xil[1]*dpotentialdy_with_rot +    xil[2]*dpotentialdl[2])
          print s_temp[1][tx] === jacobianl*(  etal[0]*dpotentialdx_with_rot +   etal[1]*dpotentialdy_with_rot +   etal[2]*dpotentialdl[2])
          print s_temp[2][tx] === jacobianl*(gammal[0]*dpotentialdx_with_rot + gammal[1]*dpotentialdy_with_rot + gammal[2]*dpotentialdl[2])
        }
      }
      print barrier(:local)
      comment()

      elem_per_thread.times { |elem_index|
        if elem_per_thread > 1 then
          print tx === get_local_id(0) + ngll3_padded * elem_index / elem_per_thread
          print k === tx/ngll2
          print j === (tx-k*ngll2)/ngllx
          print i === tx - k*ngll2 - j*ngllx
        end
        comment()

        print If(active[elem_index]) {
          (0..2).each { |indx| print templ[indx] === 0.0 }
          for_loop = For(l, 0, ngllx-1, :declit => true) {
             print templ[0] === templ[0] + s_temp[0][k*ngll2+j*ngllx+l]*sh_hprimewgll_xx[i*ngllx+l]
             print templ[1] === templ[1] + s_temp[1][k*ngll2+l*ngllx+i]*sh_hprimewgll_xx[j*ngllx+l]
             print templ[2] === templ[2] + s_temp[2][l*ngll2+j*ngllx+i]*sh_hprimewgll_xx[k*ngllx+l]
          }
          get_output.puts "#ifdef #{manually_unrolled_loops}"
            print for_loop.unroll
          get_output.puts "#else"
            print for_loop
          get_output.puts "#endif"
          comment()

          print sum_terms === -(wgllwgll_yz[k*ngllx+j]*templ[0] + wgllwgll_xz[k*ngllx+i]*templ[1] + wgllwgll_xy[j*ngllx+i]*templ[2])

          print If(gravity) {
            print  sum_terms === sum_terms + gravity_term[elem_index]
          }
          get_output.puts "#ifdef #{use_mesh_coloring}"
            get_output.puts "#ifdef #{use_textures_fields}"
              print d_potential_dot_dot[iglob[elem_index]] === d_accel_oc_tex[iglob[elem_index]] + sum_terms
            get_output.puts "#else"
              print d_potential_dot_dot[iglob[elem_index]] === d_potential_dot_dot[iglob[elem_index]] + sum_terms
            get_output.puts "#endif"
          get_output.puts "#else"
            print If(use_mesh_coloring_gpu => lambda {
              print If(nspec_outer_core > coloring_min_nspec_outer_core => lambda {
                get_output.puts "#ifdef #{use_textures_fields}"
                  print d_potential_dot_dot[iglob[elem_index]] === d_accel_oc_tex[iglob[elem_index]] + sum_terms
                get_output.puts "#else"
                  print d_potential_dot_dot[iglob[elem_index]] === d_potential_dot_dot[iglob[elem_index]] + sum_terms
                get_output.puts "#endif"
              }, :else => lambda{
                print atomicAdd(d_potential_dot_dot+iglob[elem_index],sum_terms)
              })
            }, :else => lambda {
              print atomicAdd(d_potential_dot_dot+iglob[elem_index],sum_terms)
            })
          get_output.puts "#endif"
        }
      }
      close p
    else
      raise "Unsupported language!"
    end
    pop_env(:array_start)
    kernel.procedure = p
    return kernel
  end
end
