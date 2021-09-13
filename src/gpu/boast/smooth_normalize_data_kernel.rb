module BOAST
  def BOAST::smooth_normalize_data_kernel(ref = true, n_gll3 = 125)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    function_name = "smooth_normalize_data_kernel"

    v = []
    v.push data_smooth   = Real("data_smooth", :dir => :out, :dim => [ Dim()])
    v.push normalisation = Real("normalisation", :dir => :in, :dim => [ Dim()])
    v.push nker = Int("nker", :dir => :in)
    v.push nspec_me = Int("nspec_me", :dir => :in)

    ngll3 = Int("NGLL3", :const => n_gll3)

    p = Procedure(function_name, v)

    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA or get_lang == HIP) then
      make_specfem3d_header( :ngll3 => n_gll3 )

      open p
      decl ispec =   Int("ispec")
      decl igll  =   Int("igll")
      decl norm  =   Real("norm")
      decl val   =   Real("val")

      decl tol = Real("TOL", :const => "1.e-24")

      iker = Int("iker")
      comment()

      print ispec === get_group_id(0)+get_group_id(1)*get_num_groups(0)
      print igll === get_local_id(0)
      comment()

      # note: normalization coefficient is added nker times, thus divide by nker
      print norm === normalisation[ngll3 * ispec + igll] / nker

      # avoids division by zero
      print If(norm < tol) {
        print norm === 1.0
      }
      comment()

      # normalizes smoothed kernel values
      print For(iker, 0, nker, :operator => "<", :declit => true) {
        # BOAST seems not to have an /= operator:
        #   data_smooth[ngll3*nspec_me*iker + ngll3*ispec + igll] /= norm
        # working around it..
        print val === data_smooth[ngll3*nspec_me*iker + ngll3*ispec + igll] / norm
        print data_smooth[ngll3*nspec_me*iker + ngll3*ispec + igll] === val
      }
      comment()

      close p
    else
      raise "Unsupported language!"
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end
end
