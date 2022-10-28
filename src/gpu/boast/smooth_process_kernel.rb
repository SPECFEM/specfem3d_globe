module BOAST
  def BOAST::smooth_process_kernel(ref = true, n_gll3 = 125)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    function_name = "smooth_process_kernel"

    v = []
    v.push xstore_me    = Real("xstore_me", :dir => :in, :dim => [ Dim()])
    v.push ystore_me    = Real("ystore_me", :dir => :in, :dim => [ Dim()])
    v.push zstore_me    = Real("zstore_me", :dir => :in, :dim => [ Dim()])
    v.push xstore_other = Real("xstore_other", :dir => :in, :dim => [ Dim()])
    v.push ystore_other = Real("ystore_other", :dir => :in, :dim => [ Dim()])
    v.push zstore_other = Real("zstore_other", :dir => :in, :dim => [ Dim()])
    v.push data_other   = Real("data_other", :dir => :in, :dim => [ Dim()])
    v.push sigma_h2_inv = Real("sigma_h2_inv", :dir => :in)
    v.push sigma_v2_inv = Real("sigma_v2_inv", :dir => :in)
    v.push iker         = Int("iker", :dir => :in)
    v.push nspec_me     = Int("nspec_me", :dir => :in)
    v.push nspec_other  = Int("nspec_other", :dir => :in)
    v.push v_criterion  = Real("v_criterion", :dir => :in)
    v.push h_criterion  = Real("h_criterion", :dir => :in)
    v.push integ_factor = Real("integ_factor", :dir => :in, :dim => [ Dim()])

    v.push data_smooth  = Real("data_smooth", :dir => :out, :dim => [ Dim()])
    v.push normalisation = Real("normalisation", :dir => :out, :dim => [ Dim()])
    v.push use_vector_distance = Int("use_vector_distance", :dir => :in)

    ngll3 = Int("NGLL3", :const => n_gll3)

    p = Procedure(function_name, v)

    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      make_specfem3d_header( :ngll3 => n_gll3 )

      open p
      decl ispec = Int("ispec"), igll = Int("igll")
      decl gll_other = Int("gll_other")
      decl n_loop = Int("n_loop")
      decl ispec_other = Int("ispec_other"), ispec_test = Int("ispec_test")

      decl x_me = Real("x_me"), y_me = Real("y_me"), z_me = Real("z_me")
      decl x_other = Real("x_other"), y_other = Real("y_other"), z_other = Real("z_other")
      decl center_x = Real("center_x"), center_y = Real("center_y"), center_z = Real("center_z")
      decl vx = Real("vx"), vy = Real("vy"), vz = Real("vz")
      decl alpha = Real("alpha"), ratio = Real("ratio"), theta = Real("theta")
      decl r0 = Real("r0"), r1 = Real("r1")
      decl r0_squared = Real("r0_squared"), r1_squared = Real("r1_squared")
      decl dist_h = Real("dist_h"), dist_v = Real("dist_v")
      decl val = Real("val"), val_gaussian = Real("val_gaussian")
      decl coef = Real("coef"), normalisation_slice = Real("normalisation_slice")
      decl dat = Real("dat")
      comment()

      decl sh_test     = Int("sh_test",     :local => true, :dim => [Dim(ngll3)] )
      decl sh_x_other  = Real("sh_x_other", :local => true, :dim => [Dim(ngll3)] )
      decl sh_y_other  = Real("sh_y_other", :local => true, :dim => [Dim(ngll3)] )
      decl sh_z_other  = Real("sh_z_other", :local => true, :dim => [Dim(ngll3)] )
      decl sh_integ_factor  = Real("sh_integ_factor", :local => true, :dim => [Dim(ngll3)] )
      decl sh_data          = Real("sh_data",         :local => true, :dim => [Dim(ngll3)] )
      comment()
      # PI squared
      decl pi2 = Real("PI2", :const => "9.869604401089358f")
      comment()

      print ispec === get_group_id(0)+get_group_id(1)*get_num_groups(0)
      print igll === get_local_id(0)
      comment()

      # for each reference GLL point, we can check a block of 125 neighbor elements
      # by convenience, the block size is set to the number of threads 125 of this kernel
      print n_loop === nspec_other / ngll3 + 1;

      # reference GLL point position
      print x_me === xstore_me[ngll3 * ispec + igll]
      print y_me === ystore_me[ngll3 * ispec + igll]
      print z_me === zstore_me[ngll3 * ispec + igll]

      #__syncthreads();
      print barrier(:local)
      comment()

      print dat === 0.0
      print normalisation_slice === 0.0
      comment()

      # We test 125 spectral elements at a time
      i = Int("i")
      print For(i, 0, n_loop, :operator => "<", :declit => true){
        # __syncthreads();
        print barrier(:local)
        comment()

        # each thread helps to test a different element in the other slice (using the center position)
        print ispec_other === ngll3 * i + igll

        print If(ispec_other < nspec_other){
          # center position
          print center_x === (xstore_other[ispec_other * ngll3] + xstore_other[ispec_other * ngll3 + (ngll3 - 1)]) * 0.5
          print center_y === (ystore_other[ispec_other * ngll3] + ystore_other[ispec_other * ngll3 + (ngll3 - 1)]) * 0.5
          print center_z === (zstore_other[ispec_other * ngll3] + zstore_other[ispec_other * ngll3 + (ngll3 - 1)]) * 0.5
          comment()

          # radius (squared)
          print r0_squared === x_me*x_me + y_me*y_me + z_me*z_me
          print r1_squared === center_x*center_x + center_y*center_y + center_z*center_z
          comment()

          print If(use_vector_distance => lambda {
            # vector approximation (fast computation): neglects curvature
            # radius
            print r0 === sqrt( r0_squared )
            print r1 === sqrt( r1_squared )

            # vertical distance (squared)
            print dist_v === (r1 - r0)*(r1 - r0)
            comment()

            # horizontal distance
            print alpha === r1 / r0
            print vx === alpha * x_me
            print vy === alpha * y_me
            print vz === alpha * z_me

            # vector in horizontal between new r0 and r1
            print vx === center_x - vx
            print vy === center_y - vy
            print vz === center_z - vz

            # distance is vector length
            #dist_h = sqrt( vx*vx + vy*vy + vz*vz )
            # or squared:
            print dist_h === vx*vx + vy*vy + vz*vz
          }, :else => lambda {
            # w/ exact epicentral distance calculation
            # vertical distance (squared)
            print alpha === sqrt( r0_squared * r1_squared )
            print dist_v === r1_squared + r0_squared - 2.0 * alpha
            comment()

            # epicentral distance
            print If(alpha > 0.0 => lambda {
              print ratio === (x_me*center_x + y_me*center_y + z_me*center_z) / alpha
            }, :else => lambda {
              print ratio === 1.0
            })

            # checks boundaries of ratio (due to numerical inaccuracies)
            print If(ratio >= 1.0 => lambda {
              print dist_h === 0.0
            }, (ratio <= -1.0) => lambda {
              print dist_h === r1_squared * pi2
            }, :else => lambda {
              print theta === acos( ratio )
              print dist_h === r1_squared * (theta*theta)
            })
          })
        }
        comment()

        # tests if element is too far away
        #sh_test[igll] = ( ispec_other >= nspec_other
        #                  || dist_h > h_criterion
        #                  || dist_v > v_criterion ) ? 1 : 0 ;
        print sh_test[igll] === Ternary(ispec_other >= nspec_other, 1, 0)
        print sh_test[igll] === Ternary(Or(dist_h > h_criterion,sh_test[igll]), 1, 0)
        print sh_test[igll] === Ternary(Or(dist_v > v_criterion,sh_test[igll]), 1, 0)

        #__syncthreads();
        print barrier(:local)
        comment()

        # loops over each spectral element tested
        k = Int("k")
        print For(k, 0, ngll3, :operator => "<", :declit => true){
          #__syncthreads();
          print barrier(:local)
          comment()

          # skips element if test was true (too far away)
          print If(sh_test[k]){
            BOAST::get_output.puts "        continue;"
          }
          comment()

          # loads data from other slice to shared memory
          print ispec_test === i * ngll3 + k
          print sh_x_other[igll] === xstore_other[ispec_test * ngll3 + igll]
          print sh_y_other[igll] === ystore_other[ispec_test * ngll3 + igll]
          print sh_z_other[igll] === zstore_other[ispec_test * ngll3 + igll]
          comment()

          print sh_data[igll] === data_other[ispec_test * ngll3 + igll]
          print sh_integ_factor[igll] === integ_factor[ispec_test * ngll3 + igll]
          comment()

          #__syncthreads();
          print barrier(:local)
          comment()

          # loops over gll points
          j = Int("j")
          print For(j, 0, ngll3, :operator => "<", :declit => true){
            print gll_other === Modulo(igll + j, ngll3)
            comment()

            print x_other === sh_x_other[gll_other]
            print y_other === sh_y_other[gll_other]
            print z_other === sh_z_other[gll_other]
            comment()

            # radius (squared)
            print r0_squared === x_me*x_me + y_me*y_me + z_me*z_me
            print r1_squared === x_other*x_other + y_other*y_other + z_other*z_other
            comment()

            print If(use_vector_distance => lambda {
              # vector approximation (fast computation): neglects curvature
              # radius
              print r0 === sqrt( r0_squared )
              print r1 === sqrt( r1_squared )

              # vertical distance (squared)
              print dist_v === (r1 - r0)*(r1 - r0)
              comment()

              # horizontal distance
              print alpha === r1 / r0
              print vx === alpha * x_me
              print vy === alpha * y_me
              print vz === alpha * z_me

              # vector in horizontal between new r0 and r1
              print vx === x_other - vx
              print vy === y_other - vy
              print vz === z_other - vz

              # distance is vector length
              #dist_h = sqrt( vx*vx + vy*vy + vz*vz )
              # or squared:
              print dist_h === vx*vx + vy*vy + vz*vz
            }, :else => lambda {
              # w/ exact epicentral distance calculation
              # vertical distance (squared)
              print alpha === sqrt( r0_squared * r1_squared )
              print dist_v === r1_squared + r0_squared - 2.0 * alpha
              comment()

              # epicentral distance
              print If(alpha > 0.0 => lambda {
                print ratio === (x_me*x_other + y_me*y_other + z_me*z_other) / alpha
              }, :else => lambda {
                print ratio === 1.0
              })

              # checks boundaries of ratio (due to numerical inaccuracies)
              print If(ratio >= 1.0 => lambda {
                print dist_h === 0.0
              }, (ratio <= -1.0) => lambda {
                print dist_h === r1_squared * pi2
              }, :else => lambda {
                print theta === acos(ratio)
                print dist_h === r1_squared * (theta*theta)
              })
            })
            comment()

            # Gaussian function
            print val === - dist_h*sigma_h2_inv - dist_v*sigma_v2_inv

            # limits to single precision
            print If(val < - 86.0 => lambda {
              # smaller than numerical precision: exp(-86) < 1.e-37
              print val_gaussian === 0.0
            }, :else => lambda {
              print val_gaussian === exp(val)
            })

            print coef === val_gaussian * sh_integ_factor[gll_other]

            print normalisation_slice === normalisation_slice + coef
            print dat === dat + sh_data[gll_other] * coef
          }
        }
      }

      # data_smooth[NGLL3*nspec_me*iker + NGLL3*ispec + igll] += dat;
      print val === data_smooth[ngll3*nspec_me*iker + ngll3*ispec + igll] + dat
      print data_smooth[ngll3*nspec_me*iker + ngll3*ispec + igll] === val

      # note: normalization coefficient is added nker times
      #normalisation[NGLL3*ispec + igll] += normalisation_slice;
      print val === normalisation[ngll3*ispec + igll] + normalisation_slice
      print normalisation[ngll3*ispec + igll] === val

      close p
    else
      raise "Unsupported language!"
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end
end
