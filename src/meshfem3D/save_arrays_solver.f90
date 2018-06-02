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

  subroutine save_arrays_solver(idoubling,ibool, &
                                iregion_code,xstore,ystore,zstore, &
                                NSPEC2D_TOP,NSPEC2D_BOTTOM)

  use constants

  use meshfem3D_models_par, only: &
    OCEANS,TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE, &
    ANISOTROPIC_INNER_CORE,ATTENUATION

  use meshfem3D_par, only: &
    NCHUNKS,ABSORBING_CONDITIONS,SAVE_MESH_FILES, &
    ROTATION,EXACT_MASS_MATRIX_FOR_ROTATION, &
    nspec,nglob, &
    OUTPUT_FILES,xstore_glob,ystore_glob,zstore_glob

  use regions_mesh_par2, only: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore, &
    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
    rmassx,rmassy,rmassz,rmass_ocean_load, &
    b_rmassx,b_rmassy, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
    jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
    jacobian2D_bottom,jacobian2D_top, &
    rho_vp,rho_vs, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
    ispec_is_tiso,tau_s,T_c_source,tau_e_store,Qmu_store, &
    prname

  implicit none

  ! doubling mesh flag
  integer,dimension(nspec),intent(in) :: idoubling
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  integer,intent(in) :: iregion_code

  ! arrays with the mesh in double precision
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: xstore,ystore,zstore

  ! boundary parameters locator
  integer,intent(in) :: NSPEC2D_TOP,NSPEC2D_BOTTOM

  ! local parameters
  integer :: i,j,k,ispec,iglob,ier
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: tmp_array

  ! debug file output
  character(len=MAX_STRING_LEN) :: filename
  logical,parameter :: DEBUG = .false.

  ! mesh databases
  open(unit=IOUT,file=prname(1:len_trim(prname))//'solver_data.bin', &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data.bin file')

  ! save nspec and nglob, to be used in combine_paraview_data
  write(IOUT) nspec
  write(IOUT) nglob

  ! mesh topology

  ! mesh arrays used in the solver to locate source and receivers
  ! and for anisotropy and gravity, save in single precision
  ! use tmp_array for temporary storage to perform conversion
  allocate(tmp_array(nglob),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array for mesh topology')

  !--- x coordinate
  tmp_array(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          tmp_array(iglob) = real(xstore(i,j,k,ispec), kind=CUSTOM_REAL)
        enddo
      enddo
    enddo
  enddo
  write(IOUT) tmp_array ! xstore

  !--- y coordinate
  tmp_array(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          tmp_array(iglob) = real(ystore(i,j,k,ispec), kind=CUSTOM_REAL)
        enddo
      enddo
    enddo
  enddo
  write(IOUT) tmp_array ! ystore

  !--- z coordinate
  tmp_array(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          tmp_array(iglob) = real(zstore(i,j,k,ispec), kind=CUSTOM_REAL)
        enddo
      enddo
    enddo
  enddo
  write(IOUT) tmp_array ! zstore

  deallocate(tmp_array)

  ! local to global indexing
  write(IOUT) ibool
  write(IOUT) idoubling
  write(IOUT) ispec_is_tiso

  ! local GLL points
  write(IOUT) xixstore
  write(IOUT) xiystore
  write(IOUT) xizstore
  write(IOUT) etaxstore
  write(IOUT) etaystore
  write(IOUT) etazstore
  write(IOUT) gammaxstore
  write(IOUT) gammaystore
  write(IOUT) gammazstore

  write(IOUT) rhostore
  write(IOUT) kappavstore

  ! other terms needed in the solid regions only
  if (iregion_code /= IREGION_OUTER_CORE) then

    ! note: muvstore needed for Q_mu shear attenuation in inner core
    if (.not. (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) then
      write(IOUT) muvstore
    endif

    !   save anisotropy in the mantle only
    if (TRANSVERSE_ISOTROPY) then
      if (iregion_code == IREGION_CRUST_MANTLE .and. .not. ANISOTROPIC_3D_MANTLE) then
        write(IOUT) kappahstore
        write(IOUT) muhstore
        write(IOUT) eta_anisostore
      endif
    endif

    !   save anisotropy in the inner core only
    if (ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
      write(IOUT) c11store
      write(IOUT) c33store
      write(IOUT) c12store
      write(IOUT) c13store
      write(IOUT) c44store
    endif

    if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
        write(IOUT) c11store
        write(IOUT) c12store
        write(IOUT) c13store
        write(IOUT) c14store
        write(IOUT) c15store
        write(IOUT) c16store
        write(IOUT) c22store
        write(IOUT) c23store
        write(IOUT) c24store
        write(IOUT) c25store
        write(IOUT) c26store
        write(IOUT) c33store
        write(IOUT) c34store
        write(IOUT) c35store
        write(IOUT) c36store
        write(IOUT) c44store
        write(IOUT) c45store
        write(IOUT) c46store
        write(IOUT) c55store
        write(IOUT) c56store
        write(IOUT) c66store
    endif

  endif

  ! Stacey
  if (ABSORBING_CONDITIONS) then
    if (iregion_code == IREGION_CRUST_MANTLE) then
      write(IOUT) rho_vp
      write(IOUT) rho_vs
    else if (iregion_code == IREGION_OUTER_CORE) then
      write(IOUT) rho_vp
    endif
  endif

  ! mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  if (((NCHUNKS /= 6 .and. ABSORBING_CONDITIONS) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) .and. iregion_code == IREGION_INNER_CORE)) then
     write(IOUT) rmassx
     write(IOUT) rmassy
  endif

  write(IOUT) rmassz

  ! mass matrices for backward simulation when ROTATION is .true.
  if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
    if (iregion_code == IREGION_CRUST_MANTLE .or. iregion_code == IREGION_INNER_CORE) then
       write(IOUT) b_rmassx
       write(IOUT) b_rmassy
    endif
  endif

  ! additional ocean load mass matrix if oceans and if we are in the crust
  if (OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) write(IOUT) rmass_ocean_load

  close(IOUT) ! solver_data.bin

  ! absorbing boundary parameters
  open(unit=IOUT,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening boundary.bin file')

  write(IOUT) nspec2D_xmin
  write(IOUT) nspec2D_xmax
  write(IOUT) nspec2D_ymin
  write(IOUT) nspec2D_ymax
  write(IOUT) NSPEC2D_BOTTOM
  write(IOUT) NSPEC2D_TOP

  write(IOUT) ibelm_xmin
  write(IOUT) ibelm_xmax
  write(IOUT) ibelm_ymin
  write(IOUT) ibelm_ymax
  write(IOUT) ibelm_bottom
  write(IOUT) ibelm_top

  write(IOUT) normal_xmin
  write(IOUT) normal_xmax
  write(IOUT) normal_ymin
  write(IOUT) normal_ymax
  write(IOUT) normal_bottom
  write(IOUT) normal_top

  write(IOUT) jacobian2D_xmin
  write(IOUT) jacobian2D_xmax
  write(IOUT) jacobian2D_ymin
  write(IOUT) jacobian2D_ymax
  write(IOUT) jacobian2D_bottom
  write(IOUT) jacobian2D_top

  close(IOUT)

  if (ATTENUATION) then
    open(unit=IOUT, file=prname(1:len_trim(prname))//'attenuation.bin', &
          status='unknown', form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening attenuation.bin file')
    write(IOUT) tau_s
    write(IOUT) tau_e_store
    write(IOUT) Qmu_store
    write(IOUT) T_c_source
    close(IOUT)
  endif

  if (HETEROGEN_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    open(unit=IOUT,file=prname(1:len_trim(prname))//'dvp.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening dvp.bin file')

    write(IOUT) dvpstore
    close(IOUT)
  endif

  ! uncomment for vp & vs model storage
  if (SAVE_MESH_FILES) then
    ! outputs model files in binary format
    call save_arrays_solver_meshfiles()
  endif

  ! debug outputs flags as vtk file
  if (DEBUG) then
    if (iregion_code == IREGION_CRUST_MANTLE) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/ispec_is_tiso',myrank
      call write_VTK_data_elem_l(nspec,nglob,xstore_glob,ystore_glob,zstore_glob,ibool, &
                                 ispec_is_tiso,filename)
    endif
  endif

  end subroutine save_arrays_solver

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_arrays_solver_meshfiles()

! outputs model files in binary format

  use constants

  use meshfem3D_par, only: nspec

  use meshfem3D_models_par, only: &
    TRANSVERSE_ISOTROPY,ATTENUATION,ATTENUATION_3D,ATTENUATION_1D_WITH_3D_STORAGE

  use regions_mesh_par2, only: &
    rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    Qmu_store, &
    prname

  implicit none

  ! local parameters
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store

  ! scaling factors to re-dimensionalize units
  scaleval1 = sngl( sqrt(PI*GRAV*RHOAV)*(R_EARTH/1000.0d0) )
  scaleval2 = sngl( RHOAV/1000.0d0 )

  ! uses temporary array
  allocate(temp_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating temp_store array'
  temp_store(:,:,:,:) = 0._CUSTOM_REAL

  ! isotropic model
  ! vp
  open(unit=IOUT,file=prname(1:len_trim(prname))//'vp.bin', &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening vp.bin file')

  temp_store(:,:,:,:) = sqrt((kappavstore(:,:,:,:) + 4.0_CUSTOM_REAL * muvstore(:,:,:,:)/3.0_CUSTOM_REAL)/rhostore(:,:,:,:)) &
                        * scaleval1
  write(IOUT) temp_store
  close(IOUT)
  ! vs
  open(unit=IOUT,file=prname(1:len_trim(prname))//'vs.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening vs.bin file')

  temp_store(:,:,:,:) = sqrt( muvstore(:,:,:,:)/rhostore(:,:,:,:) )*scaleval1
  write(IOUT) temp_store
  close(IOUT)
  ! rho
  open(unit=IOUT,file=prname(1:len_trim(prname))//'rho.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening rho.bin file')

  temp_store(:,:,:,:) = rhostore(:,:,:,:) * scaleval2
  write(IOUT) temp_store
  close(IOUT)

  ! transverse isotropic model
  if (TRANSVERSE_ISOTROPY) then
    ! vpv
    open(unit=IOUT,file=prname(1:len_trim(prname))//'vpv.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening vpv.bin file')

    temp_store(:,:,:,:) = sqrt((kappavstore(:,:,:,:) + 4.0_CUSTOM_REAL * muvstore(:,:,:,:)/3.0_CUSTOM_REAL)/rhostore(:,:,:,:)) &
                          * scaleval1
    write(IOUT) temp_store
    close(IOUT)

    ! vph
    open(unit=IOUT,file=prname(1:len_trim(prname))//'vph.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening vph.bin file')

    temp_store(:,:,:,:) = sqrt((kappahstore(:,:,:,:) + 4.0_CUSTOM_REAL * muhstore(:,:,:,:)/3.0_CUSTOM_REAL)/rhostore(:,:,:,:)) &
                          * scaleval1
    write(IOUT) temp_store
    close(IOUT)

    ! vsv
    open(unit=IOUT,file=prname(1:len_trim(prname))//'vsv.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening vsv.bin file')

    temp_store(:,:,:,:) = sqrt( muvstore(:,:,:,:)/rhostore(:,:,:,:) )*scaleval1
    write(IOUT) temp_store
    close(IOUT)

    ! vsh
    open(unit=IOUT,file=prname(1:len_trim(prname))//'vsh.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening vsh.bin file')

    temp_store(:,:,:,:) = sqrt( muhstore(:,:,:,:)/rhostore(:,:,:,:) )*scaleval1
    write(IOUT) temp_store
    close(IOUT)

    ! rho
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rho.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening rho.bin file')

    temp_store(:,:,:,:) = rhostore(:,:,:,:) * scaleval2
    write(IOUT) temp_store
    close(IOUT)

    ! eta
    open(unit=IOUT,file=prname(1:len_trim(prname))//'eta.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening eta.bin file')
    write(IOUT) eta_anisostore
    close(IOUT)
  endif ! TRANSVERSE_ISOTROPY

  ! shear attenuation
  if (ATTENUATION) then
    ! saves Qmu_store to full CUSTOM_REAL array
    if (ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      ! attenuation arrays are fully 3D
      temp_store(:,:,:,:) = Qmu_store(:,:,:,:)
    else
      ! attenuation array dimensions: Q_mustore(1,1,1,nspec)
      do ispec = 1,nspec
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              temp_store(i,j,k,ispec) = Qmu_store(1,1,1,ispec)
            enddo
          enddo
        enddo
      enddo
    endif

    ! Qmu
    open(unit=IOUT,file=prname(1:len_trim(prname))//'qmu.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening qmu.bin file')
    write(IOUT) temp_store
    close(IOUT)
  endif ! ATTENUATION

  ! frees temporary memory
  deallocate(temp_store)

  end subroutine save_arrays_solver_meshfiles
!
!-------------------------------------------------------------------------------------------------
!


  subroutine save_arrays_solver_MPI(iregion_code)

  use meshfem3D_par, only: &
    LOCAL_PATH, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    ADIOS_FOR_MPI_ARRAYS

!  use MPI_interfaces_par

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  implicit none

  integer,intent(in):: iregion_code

  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust mantle
    if (ADIOS_FOR_MPI_ARRAYS) then
      call save_MPI_arrays_adios(IREGION_CRUST_MANTLE,LOCAL_PATH, &
                                  num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                                  my_neighbors_crust_mantle,nibool_interfaces_crust_mantle, &
                                  ibool_interfaces_crust_mantle, &
                                  nspec_inner_crust_mantle,nspec_outer_crust_mantle, &
                                  num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
                                  num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
                                  num_elem_colors_crust_mantle)
    else
      call save_MPI_arrays(IREGION_CRUST_MANTLE,LOCAL_PATH, &
                            num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                            my_neighbors_crust_mantle,nibool_interfaces_crust_mantle, &
                            ibool_interfaces_crust_mantle, &
                            nspec_inner_crust_mantle,nspec_outer_crust_mantle, &
                            num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
                            num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
                            num_elem_colors_crust_mantle)
    endif

  case (IREGION_OUTER_CORE)
    ! outer core
    if (ADIOS_FOR_MPI_ARRAYS) then
      call save_MPI_arrays_adios(IREGION_OUTER_CORE,LOCAL_PATH, &
                                  num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                  my_neighbors_outer_core,nibool_interfaces_outer_core, &
                                  ibool_interfaces_outer_core, &
                                  nspec_inner_outer_core,nspec_outer_outer_core, &
                                  num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
                                  num_colors_outer_outer_core,num_colors_inner_outer_core, &
                                  num_elem_colors_outer_core)
    else
      call save_MPI_arrays(IREGION_OUTER_CORE,LOCAL_PATH, &
                            num_interfaces_outer_core,max_nibool_interfaces_oc, &
                            my_neighbors_outer_core,nibool_interfaces_outer_core, &
                            ibool_interfaces_outer_core, &
                            nspec_inner_outer_core,nspec_outer_outer_core, &
                            num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
                            num_colors_outer_outer_core,num_colors_inner_outer_core, &
                            num_elem_colors_outer_core)
    endif

  case (IREGION_INNER_CORE)
    ! inner core
    if (ADIOS_FOR_MPI_ARRAYS) then
      call save_MPI_arrays_adios(IREGION_INNER_CORE,LOCAL_PATH, &
                                  num_interfaces_inner_core,max_nibool_interfaces_ic, &
                                  my_neighbors_inner_core,nibool_interfaces_inner_core, &
                                  ibool_interfaces_inner_core, &
                                  nspec_inner_inner_core,nspec_outer_inner_core, &
                                  num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                                  num_colors_outer_inner_core,num_colors_inner_inner_core, &
                                  num_elem_colors_inner_core)
    else
      call save_MPI_arrays(IREGION_INNER_CORE,LOCAL_PATH, &
                            num_interfaces_inner_core,max_nibool_interfaces_ic, &
                            my_neighbors_inner_core,nibool_interfaces_inner_core, &
                            ibool_interfaces_inner_core, &
                            nspec_inner_inner_core,nspec_outer_inner_core, &
                            num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                            num_colors_outer_inner_core,num_colors_inner_inner_core, &
                            num_elem_colors_inner_core)
    endif
  end select

  end subroutine save_arrays_solver_MPI


!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_MPI_arrays(iregion_code,LOCAL_PATH, &
                                  num_interfaces,max_nibool_interfaces, &
                                  my_neighbors,nibool_interfaces, &
                                  ibool_interfaces, &
                                  nspec_inner,nspec_outer, &
                                  num_phase_ispec,phase_ispec_inner, &
                                  num_colors_outer,num_colors_inner, &
                                  num_elem_colors)
  use constants

  implicit none

  integer :: iregion_code

  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! MPI interfaces
  integer :: num_interfaces,max_nibool_interfaces
  integer, dimension(num_interfaces) :: my_neighbors
  integer, dimension(num_interfaces) :: nibool_interfaces
  integer, dimension(max_nibool_interfaces,num_interfaces) :: &
    ibool_interfaces

  ! inner/outer elements
  integer :: nspec_inner,nspec_outer
  integer :: num_phase_ispec
  integer,dimension(num_phase_ispec,2) :: phase_ispec_inner

  ! mesh coloring
  integer :: num_colors_outer,num_colors_inner
  integer, dimension(num_colors_outer + num_colors_inner) :: &
    num_elem_colors

  ! local parameters
  character(len=MAX_STRING_LEN) :: prname
  integer :: ier

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

  open(unit=IOUT,file=prname(1:len_trim(prname))//'solver_data_mpi.bin', &
       status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data_mpi.bin')

  ! MPI interfaces
  write(IOUT) num_interfaces
  if (num_interfaces > 0) then
    write(IOUT) max_nibool_interfaces
    write(IOUT) my_neighbors
    write(IOUT) nibool_interfaces
    write(IOUT) ibool_interfaces
  endif

  ! inner/outer elements
  write(IOUT) nspec_inner,nspec_outer
  write(IOUT) num_phase_ispec
  if (num_phase_ispec > 0 ) write(IOUT) phase_ispec_inner

  ! mesh coloring
  if (USE_MESH_COLORING_GPU) then
    write(IOUT) num_colors_outer,num_colors_inner
    write(IOUT) num_elem_colors
  endif

  close(IOUT)

  end subroutine save_MPI_arrays


!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_arrays_boundary()

! saves arrays for boundaries such as MOHO, 400 and 670 discontinuities

  use constants, only: IOUT

  use meshfem3d_par, only: &
    myrank

  use meshfem3D_models_par, only: &
    SAVE_BOUNDARY_MESH,HONOR_1D_SPHERICAL_MOHO,SUPPRESS_CRUSTAL_MESH

! boundary kernels
  use regions_mesh_par2, only: &
    NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670, &
    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
    ibelm_670_top,ibelm_670_bot,normal_moho,normal_400,normal_670, &
    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
    ispec2D_670_top,ispec2D_670_bot, &
    prname

  implicit none

  ! local parameters
  integer :: ier

  ! first check the number of surface elements are the same for Moho, 400, 670
  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    if (ispec2D_moho_top /= NSPEC2D_MOHO .or. ispec2D_moho_bot /= NSPEC2D_MOHO) &
           call exit_mpi(myrank, 'Not the same number of Moho surface elements')
  endif
  if (ispec2D_400_top /= NSPEC2D_400 .or. ispec2D_400_bot /= NSPEC2D_400) &
           call exit_mpi(myrank,'Not the same number of 400 surface elements')
  if (ispec2D_670_top /= NSPEC2D_670 .or. ispec2D_670_bot /= NSPEC2D_670) &
           call exit_mpi(myrank,'Not the same number of 670 surface elements')

  ! writing surface topology databases
  open(unit=IOUT,file=prname(1:len_trim(prname))//'boundary_disc.bin', &
       status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening boundary_disc.bin file')

  write(IOUT) NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670

  write(IOUT) ibelm_moho_top
  write(IOUT) ibelm_moho_bot

  write(IOUT) ibelm_400_top
  write(IOUT) ibelm_400_bot

  write(IOUT) ibelm_670_top
  write(IOUT) ibelm_670_bot

  write(IOUT) normal_moho
  write(IOUT) normal_400
  write(IOUT) normal_670

  close(IOUT)

  end subroutine save_arrays_boundary

