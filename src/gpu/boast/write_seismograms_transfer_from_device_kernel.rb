module BOAST

  def BOAST::write_seismograms_transfer_from_device_kernel(ref = true, n_gll3 = 125)
    BOAST::write_seismograms_from_device_kernel(:transfer, ref, n_gll3)
  end

  def BOAST::write_seismograms_from_device_kernel(type, ref = true, n_gll3 = 125)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    if type == :transfer then
      function_name = "write_seismograms_transfer_from_device_kernel"
      station_field =   Real("station_seismo_field",   :dir => :out,:dim => [Dim()] )
    elsif type == :transfer_strain then
      function_name = "write_seismograms_transfer_strain_from_device_kernel"
      station_field =   Real("station_strain_field",   :dir => :out,:dim => [Dim()] )
    else
      raise "Unsupported type : #{type}!"
    end

    number_receiver_global = Int("number_receiver_global",  :dir => :in, :dim => [Dim()] )
    ispec_selected_rec     = Int("ispec_selected_rec",      :dir => :in, :dim => [Dim()] )
    ibool                  = Int("ibool",                   :dir => :in, :dim => [Dim()] )
    d_field                = Real("d_field",                :dir => :in, :dim => [Dim()] )
    nrec_local             = Int("nrec_local",              :dir => :in)

    ngll3                  = Int("NGLL3",                   :const => n_gll3)

    p = Procedure(function_name, [number_receiver_global,ispec_selected_rec,ibool,station_field,d_field,nrec_local])
    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ngll3 => n_gll3 )
      open p
      blockID    = Int("blockID")
      tx         = Int("tx")
      irec       = Int("irec")
      ispec      = Int("ispec")
      if type == :transfer then
        iglob      = Int("iglob")
        decl iglob
      end
      decl tx, irec, ispec
      decl blockID
      comment()

      print blockID === get_group_id(0)+get_group_id(1)*get_num_groups(0)
      print tx      === get_local_id(0)
      comment()

      print If(blockID<nrec_local) {
        print irec === number_receiver_global[blockID] - 1
        print ispec === ispec_selected_rec[irec] - 1
        comment()

        # copy element values to output array
        if type == :transfer then
          # from displ array (NDIM,NGLOB)
          print iglob === ibool[tx + ngll3*ispec] - 1
          (0..2).each { |i|
            print station_field[ngll3*3*blockID + tx*3 + i] === d_field[iglob*3+i]
          }
        else
          # from strain epsilon_xx arrays (NGLLX,NGLLY,NGLLZ,nspec)
          print station_field[ngll3*blockID + tx] === d_field[ngll3*ispec + tx];
        end
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

