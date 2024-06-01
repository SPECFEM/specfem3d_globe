!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

  use manager_adios

#ifdef USE_XSMM
  use my_libxsmm
#endif

  implicit none

  ! synchronize all processes, waits until all processes have written their seismograms
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'finalizing simulation'
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! closes Stacey absorbing boundary snapshots
  if (ABSORBING_CONDITIONS) then
    ! crust mantle
    if (num_abs_boundary_faces_crust_mantle > 0 .and. SAVE_STACEY) call close_file_abs(0)

    ! outer core
    if (num_abs_boundary_faces_outer_core > 0 .and. SAVE_STACEY) call close_file_abs(4)

    ! frees memory
    deallocate(absorb_buffer_crust_mantle)
    if (NSPEC_OUTER_CORE > 0) deallocate(absorb_buffer_outer_core)
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

  ! finalize full gravity
  if (FULL_GRAVITY) then
    call SIEM_finalize()
  endif

  ! adios finalizes
  if (ADIOS_ENABLED) then
    call finalize_adios()
  endif

  ! asdf finalizes
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3 .and. READ_ADJSRC_ASDF) then
    call asdf_cleanup()
  endif

#ifdef USE_XSMM
  ! finalizes LIBXSMM
  call libxsmm_finalize()
#endif

  ! synchronize all
  call synchronize_all()

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
  use specfem_par_noise
  use specfem_par_movie
  implicit none

  ! from here on, no gpu data is needed anymore
  ! frees allocated memory on GPU
  if (GPU_MODE) call prepare_cleanup_device(Mesh_pointer,NCHUNKS_VAL)

  ! mass matrices
  ! crust/mantle
  if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL)) then
    deallocate(rmassx_crust_mantle,rmassy_crust_mantle)
  else
    nullify(rmassx_crust_mantle,rmassy_crust_mantle)
  endif
  if (SIMULATION_TYPE == 3) then
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
      deallocate(b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    else
      nullify(b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    endif
    nullify(b_rmassz_crust_mantle)
  endif
  deallocate(rstore_crust_mantle)

  ! optimized arrays
  if (USE_DEVILLE_PRODUCTS_VAL) then
    deallocate(sum_terms_crust_mantle,sum_terms_inner_core,sum_terms_outer_core)
    deallocate(deriv_mapping_crust_mantle,deriv_mapping_inner_core,deriv_mapping_outer_core)
  endif
  if (use_inversed_arrays) then
    deallocate(ibool_inv_tbl_crust_mantle,ibool_inv_tbl_inner_core,ibool_inv_tbl_outer_core)
    deallocate(ibool_inv_st_crust_mantle,ibool_inv_st_inner_core,ibool_inv_st_outer_core)
    deallocate(phase_iglob_crust_mantle,phase_iglob_inner_core,phase_iglob_outer_core)
  endif

  ! outer core
  if (SIMULATION_TYPE == 3 ) nullify(b_rmass_outer_core)
  deallocate(rstore_outer_core)
  deallocate(gravity_pre_store_outer_core)

  ! inner core
  if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
    deallocate(rmassx_inner_core,rmassy_inner_core)
  else
    nullify(rmassx_inner_core,rmassy_inner_core)
  endif

  if (SIMULATION_TYPE == 3) then
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
      deallocate(b_rmassx_inner_core,b_rmassy_inner_core)
    else
      nullify(b_rmassx_inner_core,b_rmassy_inner_core)
    endif
    nullify(b_rmassz_inner_core)
  endif
  deallocate(rstore_inner_core)

  ! MPI buffers
  if (.not. USE_CUDA_AWARE_MPI) then
    deallocate(buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle)
    deallocate(buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core)
    deallocate(buffer_send_vector_inner_core,buffer_recv_vector_inner_core)
  endif
  deallocate(request_send_vector_cm,request_recv_vector_cm)
  deallocate(request_send_scalar_oc,request_recv_scalar_oc)
  deallocate(request_send_vector_ic,request_recv_vector_ic)

  if (SIMULATION_TYPE == 3) then
    if (.not. USE_CUDA_AWARE_MPI) then
      deallocate(b_buffer_send_vector_cm,b_buffer_recv_vector_cm)
      deallocate(b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core)
      deallocate(b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core)
    endif
    deallocate(b_request_send_vector_cm,b_request_recv_vector_cm)
    deallocate(b_request_send_scalar_oc,b_request_recv_scalar_oc)
    deallocate(b_request_send_vector_ic,b_request_recv_vector_ic)
  endif

  if (allocated(my_neighbors_crust_mantle)) then
    deallocate(my_neighbors_crust_mantle,nibool_interfaces_crust_mantle)
    deallocate(ibool_interfaces_crust_mantle)
  endif
  if (allocated(my_neighbors_outer_core)) then
    deallocate(my_neighbors_outer_core,nibool_interfaces_outer_core)
    deallocate(ibool_interfaces_outer_core)
  endif
  if (allocated(my_neighbors_inner_core)) then
    deallocate(my_neighbors_inner_core,nibool_interfaces_inner_core)
    deallocate(ibool_interfaces_inner_core)
  endif

  ! inner/outer elements
  if (allocated(phase_ispec_inner_crust_mantle)) deallocate(phase_ispec_inner_crust_mantle)
  if (allocated(phase_ispec_inner_outer_core)) deallocate(phase_ispec_inner_outer_core)
  if (allocated(phase_ispec_inner_inner_core)) deallocate(phase_ispec_inner_inner_core)

  ! coloring
  if (allocated(num_elem_colors_crust_mantle)) deallocate(num_elem_colors_crust_mantle)
  if (allocated(num_elem_colors_outer_core)) deallocate(num_elem_colors_outer_core)
  if (allocated(num_elem_colors_inner_core)) deallocate(num_elem_colors_inner_core)

  ! sources
  deallocate(islice_selected_source, &
             ispec_selected_source, &
             Mxx,Myy,Mzz,Mxy,Mxz,Myz)
  deallocate(xi_source,eta_source,gamma_source)
  deallocate(tshift_src,hdur,hdur_Gaussian)
  deallocate(nu_source)
  ! point force source
  if (USE_FORCE_POINT_SOURCE) then
    deallocate(force_stf)
    deallocate(factor_force_source)
    deallocate(comp_dir_vect_source_E)
    deallocate(comp_dir_vect_source_N)
    deallocate(comp_dir_vect_source_Z_UP)
  endif
  deallocate(sourcearrays)
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    deallocate(iadj_vec)
    if (nadj_rec_local > 0) then
      deallocate(source_adjoint)
      if (IO_ASYNC_COPY .and. NSTEP_SUB_ADJ > 1 ) deallocate(buffer_source_adjoint)
      deallocate(iadjsrc,iadjsrc_len)
    endif
  endif

  ! receivers
  deallocate(islice_selected_rec,ispec_selected_rec, &
             xi_receiver,eta_receiver,gamma_receiver)
  if (myrank == 0 .and. WRITE_SEISMOGRAMS_BY_MAIN) deallocate(islice_num_rec_local)
  deallocate(station_name,network_name, &
             stlat,stlon,stele,stbur)
  deallocate(nu_rec,number_receiver_global)
  if (nrec_local > 0) then
    deallocate(hxir_store, &
               hetar_store, &
               hgammar_store)
    if (SIMULATION_TYPE == 2) then
      deallocate(moment_der,stshift_der)
    endif
  endif
  deallocate(seismograms)
  if (SAVE_SEISMOGRAMS_STRAIN) deallocate(seismograms_eps)

  ! kernels
  if (SIMULATION_TYPE == 3) then
    if (APPROXIMATE_HESS_KL) then
      deallocate(hess_kl_crust_mantle)
      deallocate(hess_rho_kl_crust_mantle)
      deallocate(hess_kappa_kl_crust_mantle)
      deallocate(hess_mu_kl_crust_mantle)
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
    deallocate(mask_3dmovie,muvstore_crust_mantle_3dmovie)
  endif

  ! noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then
    deallocate(noise_sourcearray, &
               normal_x_noise,normal_y_noise,normal_z_noise, &
               mask_noise,noise_surface_movie)
    ! file i/o buffer
    deallocate(noise_buffer)
  endif

  ! oceans
  if (OCEANS_VAL) then
    ! frees memory
    deallocate(ibool_ocean_load,rmass_ocean_load_selected,normal_ocean_load)
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
