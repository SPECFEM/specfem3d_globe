module BOAST

  def BOAST::resort_array(ref = true, n_gll3 = 125)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    function_name = "resort_array"

    v = []
    v.push old_array  = Real("old_array",  :dir => :inout,:dim => [Dim()] )
    v.push nspec      = Int("NSPEC",       :dir => :in)

    ngll3 = Int("NGLL3", :const => n_gll3)

    p = Procedure(function_name, v)

    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      make_specfem3d_header( :ngll3 => n_gll3 )

      open p
      decl ispec =  Int("ispec")  # using default int (signed), otherwise we would compare ispec (unsigned) < NSPEC (signed)
      decl id =     Int("id",       :signed => false)
      decl idx =    Int("idx",      :signed => false)
      decl t_idx =  Int("t_idx",    :signed => false)
      decl tx =     Int("tx",       :signed => false)
      decl offset = Int("offset",   :signed => false)
      decl sh_tmp = Real("sh_tmp",  :local => true, :dim => [Dim(21*n_gll3)] )
      i = Int("i")
      comment()

      print ispec === get_group_id(0) + get_group_id(1)*get_num_groups(0)
      comment()

      print If(ispec < nspec) {
        print tx === get_local_id(0)
        print offset === ispec*ngll3*21+tx
        print For(i, 0, 21, :operator => "<", :declit => true) {
            print sh_tmp[i*ngll3+tx] === old_array[i*ngll3+offset]
        }
      }
      # synchronizes threads
      print barrier(:local)
      comment()

      print If(ispec < nspec) {
        print For(i, 0, 21, :operator => "<", :declit => true) {
          print id === (i*ngll3+tx)
          print idx === id / 21
          print t_idx === Modulo(id, 21)
          print old_array[i*ngll3+offset] === sh_tmp[idx + t_idx*ngll3 ]
        }
      }

      close p
    else
      raise "Unsupported language!"
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end
end
