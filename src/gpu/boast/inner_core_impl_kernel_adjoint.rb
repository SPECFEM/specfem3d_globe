module BOAST
  def BOAST::inner_core_impl_kernel_adjoint(ref = true, elem_per_thread = 1, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, n_sls = 3, coloring_min_nspec_inner_core = 1000, i_flag_in_fictitious_cube = 11)

    # setting default values
    launch_bounds = false
    min_blocks = 7

    # adjoint iso/tiso kernel (default)
    kernel = BOAST::impl_kernel(:inner_core, false, ref, elem_per_thread, mesh_coloring, textures_fields, textures_constants, unroll_loops, n_gllx, n_gll2, n_gll3, n_gll3_padded, n_sls, coloring_min_nspec_inner_core, i_flag_in_fictitious_cube, launch_bounds, min_blocks, false)

    # adjoint aniso kernel
    kernel_aniso = BOAST::impl_kernel(:inner_core, false, ref, elem_per_thread, mesh_coloring, textures_fields, textures_constants, unroll_loops, n_gllx, n_gll2, n_gll3, n_gll3_padded, n_sls, coloring_min_nspec_inner_core, i_flag_in_fictitious_cube, launch_bounds, min_blocks, true)

    # output both kernels
    if (get_lang == CL) then
      # OpenCL will need for each kernel a full output of subroutines to define a (const char*) variable for each kernel
      kernel_total = [kernel,kernel_aniso]
    else
      # CUDA / HIP
      # to create a single kernel with both procedure's codes
      # (since CUDA async memcpy headers can only appear in a single file)
      str = StringIO::new
      s1 = "" + kernel.to_s + "\n"
      s2 = "" + kernel_aniso.to_s + "\n"
      str << s1 + s2

      kernel_total = CKernel::new(code: str)
      # will need to set a procedure for file outputs
      kernel_total.procedure = [kernel.procedure,kernel_aniso.procedure]
    end

    return kernel_total
  end
end

