!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine finalize_simulation()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie
  implicit none

  ! synchronize all processes, waits until all processes have written their seismograms
  call synchronize_all()

  ! closes Stacey absorbing boundary snapshots
  if (ABSORBING_CONDITIONS) then
    ! crust mantle
    if (nspec2D_xmin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
      .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
      call close_file_abs(0)
    endif

    if (nspec2D_xmax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
      .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
      call close_file_abs(1)
    endif

    if (nspec2D_ymin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
      .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
      call close_file_abs(2)
    endif

    if (nspec2D_ymax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
      .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
      call close_file_abs(3)
    endif

    ! outer core
    if (nspec2D_xmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
      .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
      call close_file_abs(4)
    endif

    if (nspec2D_xmax_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
      .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
      call close_file_abs(5)
    endif

    if (nspec2D_ymin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
      .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
      call close_file_abs(6)
    endif

    if (nspec2D_ymax_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
      .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
      call close_file_abs(7)
    endif

    if (nspec2D_zmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
      .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
      call close_file_abs(8)
    endif

    ! frees memory
    deallocate(absorb_xmin_crust_mantle, &
               absorb_xmax_crust_mantle, &
               absorb_ymin_crust_mantle, &
               absorb_ymax_crust_mantle, &
               absorb_xmin_outer_core, &
               absorb_xmax_outer_core, &
               absorb_ymin_outer_core, &
               absorb_ymax_outer_core, &
               absorb_zmin_outer_core)
  endif

  ! save/read the surface movie using the same c routine as we do for absorbing boundaries (file ID is 9)
  if (NOISE_TOMOGRAPHY /= 0) then
    call close_file_abs(9)
  endif

  ! save files to local disk or tape system if restart file
  call save_forward_arrays()

  ! dump kernel arrays
  call save_kernels()

  ! VTK visualization
  if (VTK_MODE) then
    ! closes VTK window
    if (myrank == 0 ) call finish_vtkwindow()
  endif

  ! adios finalizes
  if (ADIOS_ENABLED .or. OUTPUT_SEISMOS_ASDF) then
    call adios_cleanup()
  endif

  ! frees dynamically allocated memory
  call finalize_simulation_cleanup()

  ! closes the main output file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'End of the simulation'
    write(IMAIN,*)
    call flush_IMAIN()
    close(IMAIN)
  endif

  ! synchronizes all the processes to make sure everybody has finished
  call synchronize_all()

  end subroutine finalize_simulation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_simulation_cleanup()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie
  implicit none

  ! mass matrices
  ! crust/mantle
  if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION)) then
    deallocate(rmassx_crust_mantle,rmassy_crust_mantle)
  else
    nullify(rmassx_crust_mantle,rmassy_crust_mantle)
  endif
  if (SIMULATION_TYPE == 3) then
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
      deallocate(b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    else
      nullify(b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    endif
    nullify(b_rmassz_crust_mantle)
  endif

  ! outer core
  if (SIMULATION_TYPE == 3 ) nullify(b_rmass_outer_core)

  ! inner core
  if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
    deallocate(rmassx_inner_core,rmassy_inner_core)
  else
    nullify(rmassx_inner_core,rmassy_inner_core)
  endif
  if (SIMULATION_TYPE == 3) then
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
      deallocate(b_rmassx_inner_core,b_rmassy_inner_core)
    else
      nullify(b_rmassx_inner_core,b_rmassy_inner_core)
    endif
    nullify(b_rmassz_inner_core)
  endif

  ! MPI buffers
  deallocate(buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle, &
             request_send_vector_cm,request_recv_vector_cm)
  deallocate(buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core, &
             request_send_scalar_oc,request_recv_scalar_oc)
  deallocate(buffer_send_vector_inner_core,buffer_recv_vector_inner_core, &
             request_send_vector_ic,request_recv_vector_ic)

  if (SIMULATION_TYPE == 3) then
    deallocate(b_buffer_send_vector_cm,b_buffer_recv_vector_cm, &
               b_request_send_vector_cm,b_request_recv_vector_cm)
    deallocate(b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core, &
               b_request_send_scalar_oc,b_request_recv_scalar_oc)
    deallocate(b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core, &
               b_request_send_vector_ic,b_request_recv_vector_ic)
  endif

  deallocate(my_neighbours_crust_mantle,nibool_interfaces_crust_mantle)
  deallocate(ibool_interfaces_crust_mantle)
  deallocate(my_neighbours_outer_core,nibool_interfaces_outer_core)
  deallocate(ibool_interfaces_outer_core)
  deallocate(my_neighbours_inner_core,nibool_interfaces_inner_core)
  deallocate(ibool_interfaces_inner_core)

  ! inner/outer elements
  deallocate(phase_ispec_inner_crust_mantle)
  deallocate(phase_ispec_inner_outer_core)
  deallocate(phase_ispec_inner_inner_core)

  ! coloring
  deallocate(num_elem_colors_crust_mantle)
  deallocate(num_elem_colors_outer_core)
  deallocate(num_elem_colors_inner_core)

  ! sources
  deallocate(islice_selected_source, &
             ispec_selected_source, &
             Mxx,Myy,Mzz,Mxy,Mxz,Myz)
  deallocate(xi_source,eta_source,gamma_source)
  deallocate(tshift_cmt,hdur,hdur_gaussian)
  deallocate(nu_source)

  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) deallocate(sourcearrays)
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    deallocate(iadj_vec)
    if (nadj_rec_local > 0) then
      deallocate(adj_sourcearrays)
      if (IO_ASYNC_COPY .and. NSTEP_SUB_ADJ > 1 ) deallocate(buffer_sourcearrays)
      deallocate(iadjsrc,iadjsrc_len)
    endif
  endif

  ! receivers
  deallocate(islice_selected_rec,ispec_selected_rec, &
             xi_receiver,eta_receiver,gamma_receiver)
  deallocate(station_name,network_name, &
             stlat,stlon,stele,stbur)
  deallocate(nu,number_receiver_global)
  if (nrec_local > 0) then
    deallocate(hxir_store, &
               hetar_store, &
               hgammar_store)
    deallocate(hlagrange_store)
    if (SIMULATION_TYPE == 2) then
      deallocate(moment_der,stshift_der)
    endif
  endif
  deallocate(seismograms)

  ! kernels
  if (SIMULATION_TYPE == 3) then
    if (APPROXIMATE_HESS_KL) then
      deallocate(hess_kl_crust_mantle)
    endif
    deallocate(beta_kl_outer_core)
  endif

  ! movies
  if (MOVIE_SURFACE) then
    deallocate(store_val_ux,store_val_uy,store_val_uz)
    deallocate(store_val_ux_all,store_val_uy_all,store_val_uz_all)
  endif
  if (MOVIE_VOLUME) then
    deallocate(nu_3dmovie)
  endif

  ! noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then
    deallocate(noise_sourcearray, &
              normal_x_noise,normal_y_noise,normal_z_noise, &
              mask_noise,noise_surface_movie)
  endif

  ! VTK visualization
  if (VTK_MODE) then
    ! frees memory
    deallocate(vtkdata,vtkmask)
    if (NPROCTOT_VAL > 1) then
      deallocate(vtkdata_points_all,vtkdata_offset_all)
      if (myrank == 0 ) deallocate(vtkdata_all)
    endif
  endif

  end subroutine finalize_simulation_cleanup
