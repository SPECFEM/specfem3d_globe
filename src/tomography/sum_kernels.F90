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

#ifdef ADIOS_INPUT
  use manager_adios
#endif

  implicit none

  character(len=MAX_STRING_LEN) :: kernel_list(MAX_KERNEL_PATHS), sline, kernel_name
  integer :: nker
  integer :: ier

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
#ifdef ADIOS_INPUT
    write(*,*) 'sum_kernels_globe (ADIOS version):'
#else
    write(*,*) 'sum_kernels_globe:'
#endif
    write(*,*)
    write(*,*) 'reading kernel list: '
  endif

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
    write(*,*) '  ',nker,' events'
    write(*,*)
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
    stop 'Error total number of slices'
  endif
  call synchronize_all()

#ifdef ADIOS_INPUT
  call synchronize_all()
  ! initializes ADIOS
  if (myrank == 0) then
    print *, 'initializing ADIOS...'
    print *, ' '
  endif
  call initialize_adios()
#endif

  ! user output
  if (myrank == 0) then
    print *,'summing kernels in INPUT_KERNELS/ directories:'
    print *,kernel_list(1:nker)
    print *
  endif

  ! synchronizes
  call synchronize_all()

  ! sums up kernels
  if (USE_ISO_KERNELS) then

    !  isotropic kernels
    if (myrank == 0) write(*,*) 'isotropic kernels: bulk_c, bulk_beta, rho'

    kernel_name = 'bulk_c_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'bulk_beta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'rho_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

  else if (USE_ALPHA_BETA_RHO) then

    ! isotropic kernels
    if (myrank == 0) write(*,*) 'isotropic kernels: alpha, beta, rho'

    kernel_name = 'alpha_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'beta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'rho_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

  else

    ! transverse isotropic kernels
    if (myrank == 0) write(*,*) 'transverse isotropic kernels: bulk_c, bulk_betav, bulk_betah, eta'

    kernel_name = 'bulk_c_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'bulk_betav_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'bulk_betah_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'eta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

  endif

  if (myrank == 0) write(*,*) 'done writing all kernels, see directory OUTPUT_SUM/'

#ifdef ADIOS_INPUT
  ! finalizes adios
  call finalize_adios()
#endif

  ! stop all the processes, and exit
  call finalize_mpi()

end program sum_kernels_globe

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_kernel(kernel_name,kernel_list,nker)

  use tomography_par

#ifdef ADIOS_INPUT
  use manager_adios
  use adios_helpers_mod, only: define_adios_scalar,define_adios_global_array1D
#endif

  implicit none

  character(len=MAX_STRING_LEN) :: kernel_name,kernel_list(MAX_KERNEL_PATHS)
  integer :: nker

  ! local parameters
  character(len=MAX_STRING_LEN) :: k_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel,total_kernel
  double precision :: norm,norm_sum
  integer :: iker,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: mask_source

  ! ADIOS
#ifdef ADIOS_INPUT
  integer :: is,ie
  integer :: local_dim
  integer(kind=8) :: group
  integer(kind=8),save :: group_size_inc
  character(len=MAX_STRING_LEN) :: kernel_name_adios
  logical,save :: is_first_call = .true.

  character(len=MAX_STRING_LEN),parameter :: group_name = "KERNELS_GROUP"
  character(len=MAX_STRING_LEN),dimension(10),parameter :: kl_name = &
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

  ! initializes arrays
  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           total_kernel(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'
  total_kernel(:,:,:,:) = 0._CUSTOM_REAL

  if (USE_SOURCE_MASK) then
    allocate( mask_source(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) )
    mask_source(:,:,:,:) = 1.0_CUSTOM_REAL
  endif

  ! ADIOS setup group size
#ifdef ADIOS_INPUT
  ! start setting up full file group size when first kernel is called
  if (is_first_call) then
    ! reset flag
    is_first_call = .false.

    ! sets up adios group
    group_size_inc = 0
    call init_adios_group(group,group_name)

    ! defines group size
    call define_adios_scalar(group, group_size_inc, '', "NSPEC", NSPEC_CRUST_MANTLE)

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
      call define_adios_global_array1D(group, group_size_inc, local_dim, '', trim(kl_name(iker)), total_kernel(:,:,:,:))
    enddo

    ! opens new adios model file
    write(k_file,'(a)') 'OUTPUT_SUM/' // 'kernels_sum.bp'
    call open_file_adios_write(k_file,group_name)
    call set_adios_group_size(group_size_inc)

    ! writes nspec
    call write_adios_scalar_int("NSPEC",NSPEC_CRUST_MANTLE)

    ! closes file
    call close_file_adios()
  endif

  ! choose corresponding ADIOS kernel name
  select case (trim(kernel_name))
  case ('bulk_c_kernel')
    kernel_name_adios = kl_name(1)
  case ('bulk_beta_kernel')
    kernel_name_adios = kl_name(2)
  case ('rho_kernel')
    kernel_name_adios = kl_name(3)
  case ('alpha_kernel')
    kernel_name_adios = kl_name(4)
  case ('beta_kernel')
    kernel_name_adios = kl_name(5)
  !case ('rho_kernel')
  !  kernel_name_adios = kl_name(6)
  !case ('bulk_c_kernel')
  !  kernel_name_adios = kl_name(7)
  case ('bulk_betav_kernel')
    kernel_name_adios = kl_name(8)
  case ('bulk_betah_kernel')
    kernel_name_adios = kl_name(9)
  case ('eta_kernel')
    kernel_name_adios = kl_name(10)
  case default
    print *,'Error kernel name not recognized for ADIOS: ',trim(kernel_name)
    stop 'Kernel name not recognized for ADIOS'
  end select
#endif

  ! loops over all event kernels
  do iker = 1, nker
    ! user output
    if (myrank == 0) then
#ifdef ADIOS_INPUT
    write(*,*) 'reading in event kernel for: ',trim(kernel_name_adios)
#else
    write(*,*) 'reading in event kernel for: ',trim(kernel_name)
#endif
    write(*,*) '    ',iker, ' out of ', nker
    endif

    ! sensitivity kernel / frechet derivative
#ifdef ADIOS_INPUT
    write(k_file,'(a,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)),'/kernels.bp'
    ! debug
    !write(*,*) 'adios kernel name: ',trim(kernel_name_adios)
    !write(*,*) 'adios file: ',trim(k_file)
    ! reads adios file
    call open_file_adios_read(k_file)
    call read_adios_array_gll(myrank,NSPEC_CRUST_MANTLE,trim(kernel_name_adios),kernel(:,:,:,:))
    call close_file_adios_read()
#else
    write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                               //'/proc',myrank,trim(REG)//trim(kernel_name)//'.bin'
    ! reads binary file
    call read_kernel_binary(k_file,NSPEC_CRUST_MANTLE,kernel)
#endif

    ! outputs norm of kernel
    norm = sum( kernel * kernel )
    call sum_all_dp(norm,norm_sum)
    if (myrank == 0) then
      print *,'  norm kernel: ',sqrt(norm_sum)
      print *
    endif

    ! source mask
    if (USE_SOURCE_MASK) then
      ! reads in mask
      write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                                 //'/proc',myrank,trim(REG)//'mask_source.bin'
      call read_kernel_binary(k_file,NSPEC_CRUST_MANTLE,mask_source)

      ! masks source elements
      kernel = kernel * mask_source
    endif

    ! sums all kernels from each event
    total_kernel = total_kernel + kernel
  enddo

  ! user output
#ifdef ADIOS_INPUT
  if (myrank == 0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name_adios),' into file kernels_sum.bp'
#else
  if (myrank == 0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)
#endif

  ! stores summed kernels
#ifdef ADIOS_INPUT
  ! appends to adios file
  write(k_file,'(a)') 'OUTPUT_SUM/' // 'kernels_sum.bp'
  call open_file_adios_write_append(k_file,group_name)
  call set_adios_group_size(group_size_inc)
  ! writes previously defined ADIOS variables
  call write_adios_array_gll(myrank,NSPEC_CRUST_MANTLE,trim(kernel_name_adios),total_kernel(:,:,:,:))
  ! closing performs actual write
  call close_file_adios()
#else
  write(k_file,'(a,i6.6,a)') 'OUTPUT_SUM/proc',myrank,trim(REG)//trim(kernel_name)//'.bin'
  call write_kernel_binary(k_file,NSPEC_CRUST_MANTLE,total_kernel)
#endif

  if (myrank == 0) write(*,*)

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
