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

! sum_kernels_globe
!
! this program can be used for event kernel summation,
! where it sums up transverse isotropic kernel files:
!
!   - proc***_reg1_bulk_c_kernel.bin
!   - proc***_reg1_bulk_betav_kernel.bin
!   - proc***_reg1_bulk_betah_kernel.bin
!   - proc***_reg1_eta_kernel.bin
!
! input file: kernels_list.txt
!   lists all event kernel directories which should be summed together
!
! input directory:  INPUT_KERNELS/
!    contains links to all event kernel directories (listed in "kernels_list.txt")
!
! output directory: OUTPUT_SUM/
!    the resulting kernel files will be stored in this directory


program sum_kernels_globe

  use tomography_par

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  use adios_helpers_mod
  use manager_adios
#endif

  implicit none

  character(len=MAX_STRING_LEN),dimension(:),allocatable :: kernel_list
  character(len=MAX_STRING_LEN) :: kernel_name,sline
  integer :: nker,iker
  integer :: ier,i

  ! kernel names
  ! consider adding more when implemented...
  !
  ! iso (bulk c, bulk beta, rho)
  character(len=MAX_STRING_LEN),dimension(3),parameter :: kernel_names_iso = &
      (/character(len=MAX_STRING_LEN) :: &
       "bulk_c_kernel", &
       "bulk_beta_kernel", &
       "rho_kernel" &
       /)
  ! iso (alpha,beta,rho)
  character(len=MAX_STRING_LEN),dimension(3),parameter :: kernel_names_alpha_beta_rho = &
      (/character(len=MAX_STRING_LEN) :: &
       "alpha_kernel", &
       "beta_kernel", &
       "rho_kernel" &
       /)
  ! tiso
  character(len=MAX_STRING_LEN),dimension(4),parameter :: kernel_names_tiso = &
      (/character(len=MAX_STRING_LEN) :: &
       "bulk_c_kernel", &
       "bulk_betav_kernel", &
       "bulk_betah_kernel", &
       "eta_kernel" &
       /)

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! corresponding kernel names of adios variables
  character(len=MAX_STRING_LEN),dimension(10),parameter :: adios_kl_names = &
      (/character(len=MAX_STRING_LEN) :: &
       "bulk_c_kl_crust_mantle", &
       "bulk_beta_kl_crust_mantle", &
       "rho_kl_crust_mantle", &
       "alpha_kl_crust_mantle", &
       "beta_kl_crust_mantle", &
       "rho_kl_crust_mantle", &
       "bulk_c_kl_crust_mantle", &
       "bulk_betav_kl_crust_mantle", &
       "bulk_betah_kl_crust_mantle", &
       "eta_kl_crust_mantle" &
       /) ! "rho_kl_crust_mantle","hess_kl_crust_mantle"
#endif

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    write(*,*) 'sum_kernels_globe (ADIOS version):'
#else
    write(*,*) 'sum_kernels_globe:'
#endif
    write(*,*)
    write(*,*) 'reading kernel list: ',trim(kernel_file_list)
  endif

  ! allocates array
  allocate(kernel_list(MAX_KERNEL_PATHS),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel_list array'
  kernel_list(:) = ''

  ! reads in event list
  nker = 0
  open(unit=IIN,file=trim(kernel_file_list),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening ',trim(kernel_file_list),myrank
    stop 1
  endif
  do while (1 == 1)
    read(IIN,'(a)',iostat=ier) sline
    if (ier /= 0) exit
    nker = nker+1
    if (nker > MAX_KERNEL_PATHS) stop 'Error number of kernels exceeds MAX_KERNEL_PATHS'
    kernel_list(nker) = sline
  enddo
  close(IIN)
  if (myrank == 0) then
    write(*,*) '  found ',nker,' events'
    write(*,*)
  endif
  ! check
  if (nker == 0) stop 'No event kernel directories found, please check...'

  ! reads mesh parameters
  if (myrank == 0) then
    ! reads mesh_parameters.bin file from topo/
    LOCAL_PATH = INPUT_DATABASES_DIR        ! 'topo/' should hold mesh_parameters.bin file
    call read_mesh_parameters()
  endif
  ! broadcast parameters to all processes
  call bcast_mesh_parameters()

  ! user output
  if (myrank == 0) then
    print *,'mesh parameters (from topo/ directory):'
    print *,'  NSPEC_CRUST_MANTLE = ',NSPEC_CRUST_MANTLE
    print *,'  NPROCTOT           = ',NPROCTOT_VAL
    print *
  endif

  ! checks if number of MPI process as specified
  if (sizeprocs /= NPROCTOT_VAL) then
    if (myrank == 0) then
      print *
      print *,'Error: run xsum_kernels with the same number of MPI processes '
      print *,'       as specified when slices were created'
      print *
      print *,'for example: mpirun -np ',NPROCTOT_VAL,' ./xsum_kernels ...'
      print *
    endif
    call synchronize_all()
    stop 'Error MPI processes not equal to total number of slices'
  endif
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print *,'kernel directories: ',nker
    print *,'summing kernels in following INPUT_KERNELS/ sub-directories:'
    do iker = 1,nker
      print *,'  ',trim(kernel_list(iker))
    enddo
    print *
  endif

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! initializes
  if (myrank == 0) then
    print *, 'initializing ADIOS...'
    print *, ' '
  endif
  call initialize_adios()

  ! opens file for summed kernel result
  call open_sum_file_adios(adios_kl_names)
#endif

  ! synchronizes
  call synchronize_all()

  ! sums up kernels
  if (USE_ISO_KERNELS) then

    !  isotropic kernels
    if (myrank == 0) then
      write(*,*) 'isotropic kernels: bulk_c, bulk_beta, rho'
      print *
    endif
    do i = 1,3
      ! kernel name
      kernel_name = kernel_names_iso(i)

#ifdef USE_ADIOS_INSTEAD_OF_MESH
      call get_adios_kernel_name(kernel_name,adios_kl_names)
#endif
      ! stores sum over all event kernels
      call sum_kernel(kernel_name,kernel_list,nker)
    enddo

  else if (USE_ALPHA_BETA_RHO) then

    ! isotropic kernels
    if (myrank == 0) then
      write(*,*) 'isotropic kernels: alpha, beta, rho'
      print *
    endif
    do i = 1,3
      ! kernel name
      kernel_name = kernel_names_alpha_beta_rho(i)

#ifdef USE_ADIOS_INSTEAD_OF_MESH
      call get_adios_kernel_name(kernel_name,adios_kl_names)
#endif
      ! stores sum over all event kernels
      call sum_kernel(kernel_name,kernel_list,nker)
    enddo

  else

    ! transverse isotropic kernels
    if (myrank == 0) then
      write(*,*) 'transverse isotropic kernels: bulk_c, bulk_betav, bulk_betah, eta'
      print *
    endif
    do i = 1,4
      ! kernel name
      kernel_name = kernel_names_tiso(i)

#ifdef USE_ADIOS_INSTEAD_OF_MESH
      call get_adios_kernel_name(kernel_name,adios_kl_names)
#endif
      ! stores sum over all event kernels
      call sum_kernel(kernel_name,kernel_list,nker)
    enddo

  endif

  ! synchronizes
  call synchronize_all()

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! closes summed kernel file
  call close_sum_file_adios()
  ! finalizes adios
  call finalize_adios()
#endif

  ! synchronizes
  call synchronize_all()

  if (myrank == 0) write(*,*) 'done writing all kernels, see directory OUTPUT_SUM/'

  ! stop all the processes, and exit
  call finalize_mpi()

end program sum_kernels_globe

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_kernel(kernel_name,kernel_list,nker)

  use tomography_par

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  use adios_helpers_mod
  use manager_adios
#endif

  implicit none

  character(len=MAX_STRING_LEN), intent(in) :: kernel_name,kernel_list(MAX_KERNEL_PATHS)
  integer, intent(in) :: nker

  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel,total_kernel
  double precision :: norm,norm_sum
  integer :: iker,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: mask_source

  ! initializes arrays
  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           total_kernel(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'
  kernel(:,:,:,:) = 0.0_CUSTOM_REAL
  total_kernel(:,:,:,:) = 0.0_CUSTOM_REAL

  if (USE_SOURCE_MASK) then
    allocate( mask_source(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) )
    mask_source(:,:,:,:) = 1.0_CUSTOM_REAL
  endif

  ! loops over all event kernels
  do iker = 1, nker
    ! user output
    if (myrank == 0) then
    write(*,*) 'reading in event kernel for: ',trim(kernel_name)
    write(*,*) '    ',iker, ' out of ', nker
    endif

    ! sensitivity kernel / frechet derivative
    kernel(:,:,:,:) = 0.0_CUSTOM_REAL

    ! reads in kernel values
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    ! reads existing event kernel file
    write(file_name,'(a,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)),'/kernels'
    file_name = get_adios_filename(trim(file_name))

    ! debug
    !print *,'debug: ',myrank,'adios file  : ',trim(file_name)
    !print *,'debug: ',myrank,'adios kernel: ',trim(kernel_name)

    ! reads adios file
    call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

    ! gets kernel value
    call read_adios_array(myadios_file,myadios_group,myrank,NSPEC_CRUST_MANTLE,trim(kernel_name),kernel)

    ! closes read file
    call close_file_adios_read_and_finalize_method(myadios_file)
#else
    ! default binary
    write(file_name,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                                  //'/proc',myrank,trim(REG)//trim(kernel_name)//'.bin'
    ! reads binary file
    call read_kernel_binary(file_name,NSPEC_CRUST_MANTLE,kernel)
#endif

    ! outputs norm of kernel
    norm = sum( kernel * kernel )
    call sum_all_dp(norm,norm_sum)
    if (myrank == 0) then
      print *,'  norm kernel: ',sqrt(norm_sum)
      print *
    endif
    call synchronize_all()

    ! source mask
    if (USE_SOURCE_MASK) then
      ! reads in mask
      write(file_name,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                                    //'/proc',myrank,trim(REG)//'mask_source.bin'
      call read_kernel_binary(file_name,NSPEC_CRUST_MANTLE,mask_source)

      ! masks source elements
      kernel = kernel * mask_source
    endif

    ! sums all kernels from each event
    total_kernel = total_kernel + kernel
  enddo

  ! user output
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  if (myrank == 0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name),' into file adios kernels_sum file'
#else
  if (myrank == 0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)
#endif
  call synchronize_all()

  ! stores summed kernels
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! writes summed kernel result
  call write_adios_array_gll(myadios_val_file,myadios_val_group,myrank,sizeprocs_adios,NSPEC_CRUST_MANTLE, &
                             trim(kernel_name),total_kernel(:,:,:,:))
  ! sync adios2 writes before loosing the temporary scope of the arrays
  if (is_adios_version2) call write_adios_perform(myadios_val_file)
#else
  ! default binary
  write(file_name,'(a,i6.6,a)') 'OUTPUT_SUM/proc',myrank,trim(REG)//trim(kernel_name)//'.bin'
  call write_kernel_binary(file_name,NSPEC_CRUST_MANTLE,total_kernel)
#endif

  if (myrank == 0) then
    ! min/max
    write(*,*) '  slice rank 0 has min/max value = ',minval(total_kernel(:,:,:,:)),"/",maxval(total_kernel(:,:,:,:))
    write(*,*)
  endif
  call synchronize_all()

  ! frees memory
  deallocate(kernel,total_kernel)
  if (USE_SOURCE_MASK) deallocate(mask_source)

  end subroutine sum_kernel


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_kernel_binary(filename,nspec,kernel)

  use tomography_par, only: NGLLX,NGLLY,NGLLZ,IIN,CUSTOM_REAL,MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: filename

  integer,intent(in) :: nspec
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(out) :: kernel

  ! local parameters
  integer :: ier

  ! initializes
  kernel(:,:,:,:) = 0._CUSTOM_REAL

  ! reads in binary file
  open(IIN,file=trim(filename),status='old',form='unformatted',action='read',iostat=ier)
  if (ier /= 0) then
   write(*,*) '  file not found: ',trim(filename)
   stop 'Error file not found'
  endif

  read(IIN) kernel(:,:,:,:)

  close(IIN)

  end subroutine read_kernel_binary

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernel_binary(filename,nspec,kernel)

  use tomography_par, only: NGLLX,NGLLY,NGLLZ,IOUT,CUSTOM_REAL,MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: filename

  integer,intent(in) :: nspec
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: kernel

  ! local parameters
  integer :: ier

  ! writes out in binary format
  open(IOUT,file=trim(filename),form='unformatted',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'Error file not written: ',trim(filename)
    stop 'Error file write'
  endif

  write(IOUT) kernel

  close(IOUT)

  end subroutine write_kernel_binary

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_ADIOS_INSTEAD_OF_MESH
! ADIOS only

  subroutine open_sum_file_adios(adios_kl_names)

  use tomography_par

  use adios_helpers_mod
  use manager_adios

  implicit none

  character(len=MAX_STRING_LEN),dimension(10), intent(in) :: adios_kl_names

  ! local parameters
  integer :: is,ie,iker,ier
  integer(kind=8) :: local_dim
  integer(kind=8) :: group_size_inc
  character(len=MAX_STRING_LEN) :: file_name
  character(len=MAX_STRING_LEN) :: writer_group_name
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dummy_real4d

  ! dummy for definitions
  allocate(dummy_real4d(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating dummy array'
  dummy_real4d(:,:,:,:) = 0.0_CUSTOM_REAL

  ! ADIOS
  ! start setting up full file group size
  ! i/o group for reading in event kernels
  call init_adios_group(myadios_group,"KernelReader")

  ! i/o group to write out summed kernel values
  writer_group_name = "KERNELS_GROUP"
  call init_adios_group(myadios_val_group,writer_group_name)

  ! defines variables and group size
  group_size_inc = 0
  call define_adios_scalar(myadios_val_group, group_size_inc, '', "NSPEC", NSPEC_CRUST_MANTLE)
  call define_adios_scalar(myadios_val_group, group_size_inc, '', "reg1/nspec", NSPEC_CRUST_MANTLE)

  ! defines all arrays
  if (USE_ISO_KERNELS) then
    ! for 3 parameters: (bulk_c, bulk_beta, rho)
    is = 1; ie = 3
  else if (USE_ALPHA_BETA_RHO) then
    ! for 3 parameters: (alpha, beta, rho)
    is = 4; ie = 6
  else
    ! for 4 parameters: (bulk_c, bulk_betav, bulk_betah, eta)
    is = 7; ie = 10
  endif
  do iker = is,ie
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, '', trim(adios_kl_names(iker)), dummy_real4d)
  enddo

  deallocate(dummy_real4d)

  ! opens new adios model file
  file_name = get_adios_filename('OUTPUT_SUM/' // 'kernels_sum')

  call open_file_adios_write(myadios_val_file,myadios_val_group,file_name,writer_group_name)

  call set_adios_group_size(myadios_val_file,group_size_inc)

  ! initially writes nspec (for checking and backward compatibility)
  call write_adios_scalar(myadios_val_file,myadios_val_group,"NSPEC",NSPEC_CRUST_MANTLE)
  call write_adios_scalar(myadios_val_file,myadios_val_group,"reg1/nspec",NSPEC_CRUST_MANTLE)

  end subroutine open_sum_file_adios

!-------------------------------------------------------------------------------------------------

  subroutine close_sum_file_adios()

  use manager_adios

  implicit none

  ! closes summed kernel file
  call close_file_adios(myadios_val_file)

  end subroutine close_sum_file_adios

!-------------------------------------------------------------------------------------------------

  subroutine get_adios_kernel_name(kernel_name,adios_kl_names)

  use constants, only: MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN), intent(inout) :: kernel_name
  character(len=MAX_STRING_LEN),dimension(10), intent(in) :: adios_kl_names

  ! choose corresponding ADIOS kernel name
  select case (trim(kernel_name))
  case ('bulk_c_kernel')
    kernel_name = trim(adios_kl_names(1))
  case ('bulk_beta_kernel')
    kernel_name = trim(adios_kl_names(2))
  case ('rho_kernel')
    kernel_name = trim(adios_kl_names(3))
  case ('alpha_kernel')
    kernel_name = trim(adios_kl_names(4))
  case ('beta_kernel')
    kernel_name = trim(adios_kl_names(5))
  !case ('rho_kernel')
  !  kernel_name = trim(adios_kl_names(6))
  !case ('bulk_c_kernel')
  !  kernel_name = trim(adios_kl_names(7))
  case ('bulk_betav_kernel')
    kernel_name = trim(adios_kl_names(8))
  case ('bulk_betah_kernel')
    kernel_name = trim(adios_kl_names(9))
  case ('eta_kernel')
    kernel_name = trim(adios_kl_names(10))
  case default
    print *,'Error kernel name not recognized for ADIOS: ',trim(kernel_name)
    stop 'Kernel name not recognized for ADIOS'
  end select

  end subroutine get_adios_kernel_name

#endif  /* USE_ADIOS_INSTEAD_OF_MESH */
