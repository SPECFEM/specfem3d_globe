  if(SIMULATION_TYPE > 1) then
    allocate(b_buffer_send_faces(NDIM,npoin2D_max_all_CM_IC,NUMFACES_SHARED), &
             b_buffer_received_faces(NDIM,npoin2D_max_all_CM_IC,NUMFACES_SHARED),stat=ier)
  else
! dummy allocation of unusued arrays
    allocate(b_buffer_send_faces(1,1,1), &
             b_buffer_received_faces(1,1,1),stat=ier)
  endif

  ! mass matrix including central cube
  if(INCLUDE_CENTRAL_CUBE) then

    if(myrank == 0) write(IMAIN,*) 'including central cube'

    if(SIMULATION_TYPE > 1) then
      allocate(b_buffer_all_cube_from_slices(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices,NDIM), &
               b_buffer_slices(npoin2D_cube_from_slices,NDIM))
    else
! dummy allocation of unusued arrays
      allocate(b_buffer_all_cube_from_slices(1,1,1), &
               b_buffer_slices(1,1))
    endif
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating backward cube buffers')

  else

    ! allocate fictitious buffers for cube and slices with a dummy size
    ! just to be able to use them as arguments in subroutine calls
    allocate(b_buffer_all_cube_from_slices(1,1,1), &
            b_buffer_slices(1,1),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating dummy buffers')

  endif

  if (SIMULATION_TYPE == 3) then
    b_div_displ_outer_core(:,:,:,:) = 0._CUSTOM_REAL
  endif
