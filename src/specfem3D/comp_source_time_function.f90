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

  double precision function comp_source_time_function(t,hdur,it_index)

  use constants, only: EXTERNAL_SOURCE_TIME_FUNCTION

  implicit none

  double precision,intent(in) :: t,hdur
  integer, intent(in) :: it_index

  ! local parameters
  double precision, external :: comp_source_time_function_heavi
  double precision, external :: comp_source_time_function_ext

  if (EXTERNAL_SOURCE_TIME_FUNCTION) then
    ! external stf
    comp_source_time_function = comp_source_time_function_ext(it_index)
  else
    ! quasi Heaviside
    comp_source_time_function = comp_source_time_function_heavi(t,hdur)
  endif

  end function comp_source_time_function

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_heavi(t,hdur)

  implicit none

  double precision,intent(in) :: t,hdur

  double precision, external :: netlib_specfun_erf

  ! quasi Heaviside, small Gaussian moment-rate tensor with hdur
  ! note: netlib_specfun_erf is calculating the (Gauss) error function based on netlib's implementation:
  !         http://www.netlib.org/specfun/erf
  !         https://en.wikipedia.org/wiki/Error_function
  !       the error function erf(x) would be defined between [-1,1], here we scale it to be within [0,1]
  !       to have a normalized source time function
  comp_source_time_function_heavi = 0.5d0*(1.0d0 + netlib_specfun_erf(t/hdur))

  end function comp_source_time_function_heavi


!
!-------------------------------------------------------------------------------------------------
!


  double precision function comp_source_time_function_mono(t,f0)

  use constants, only: PI,TAPER_MONOCHROMATIC_SOURCE

  implicit none

  double precision,intent(in) :: t,f0
  double precision :: tt
  integer :: taper

  tt = 2 * PI * f0 * t
  taper = ceiling(TAPER_MONOCHROMATIC_SOURCE * f0)

  if (t < taper / f0) then
    comp_source_time_function_mono = sin(tt) * (0.5 - 0.5 * cos(tt / taper / 2.0))
  else
    comp_source_time_function_mono = sin(tt)
  endif

  ! monochromatic source time function

  end function comp_source_time_function_mono

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_rickr(t,f0)

  use constants, only: PI

  implicit none

  double precision,intent(in) :: t,f0

  ! local parameters
  double precision :: a

  ! Ricker wavelet
  a = PI*PI * f0*f0

  comp_source_time_function_rickr = (1.d0 - 2.d0 * a * t*t ) * exp( -a * t*t )

  !!! another source time function they have called 'Ricker' in some old papers,
  !!! e.g., 'Finite-Frequency Kernels Based on Adjoint Methods' by Liu & Tromp, BSSA (2006)
  !!! in order to benchmark those simulations, the following formula is needed.
  ! comp_source_time_function_rickr = -2.d0*PI*PI*f0*f0*f0*t * exp(-PI*PI*f0*f0*t*t)

  end function comp_source_time_function_rickr

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_gauss(t,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur

  ! local parameters
  double precision :: hdur_decay,a

  ! note: hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.628
  hdur_decay = hdur

  ! this here uses a stronger Gaussian decay rate (empirical value) to avoid non-zero onset times;
  ! however, it should mimik a triangle source time function...
  !hdur_decay = hdur  / SOURCE_DECAY_STRONG

  ! note: a nonzero time to start the simulation with would lead to more high-frequency noise
  !          due to the (spatial) discretization of the point source on the mesh

  ! Gaussian wavelet
  a = 1.d0 / (hdur_decay**2)

  comp_source_time_function_gauss = exp(-a * t**2) / (sqrt(PI) * hdur_decay)

  end function comp_source_time_function_gauss


!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_gauss_2(t,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur

  ! local parameters
  double precision :: a

  ! source time function as defined in:
  ! M.A. Meschede, C.L. Myhrvold and J. Tromp, 2011.
  ! Antipodal focusing of seismic waves due to large meteorite impacts on Earth,
  ! GJI, 187, p. 529-537
  !
  ! equation (2):
  ! S(t) = sqrt(pi/tau**2) * exp(-pi**2 * t**2 / tau**2)

  ! factor
  a = sqrt(PI / hdur**2)

  comp_source_time_function_gauss_2 = a * exp(-PI**2 * t**2 / hdur**2)

  end function comp_source_time_function_gauss_2

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_ext(it_index)

  use specfem_par, only: SIMULATION_TYPE, NSTEP, user_source_time_function

  implicit none

  integer,intent(in) :: it_index

  ! local parameters
  integer :: it_tmp

  ! On the first iteration, go get the ASCII file.
  if (.not. allocated (user_source_time_function)) then
    call read_external_source_time_function()
  endif

  ! sets index in user stf
  if (SIMULATION_TYPE == 1) then
    ! forward simulation
    it_tmp = it_index
  else
    ! reverse time
    ! iteration step
    ! (see compute_add_sources_backward() in compute_add_sources.f90)
    ! to match a reversed time value: time_t = dble(NSTEP-it_tmp)*DT - t0
    it_tmp = NSTEP - it_index
  endif

  ! checks bounds
  if (it_tmp < 1) it_tmp = 1
  if (it_tmp > NSTEP) it_tmp = NSTEP

  ! gets stored STF
  comp_source_time_function_ext = user_source_time_function(it_tmp)

  end function comp_source_time_function_ext

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_external_source_time_function()

  use constants, only: IIN, MAX_STRING_LEN
  use specfem_par, only: NSTEP, user_source_time_function, USE_LDDRK

  implicit none

  ! local parameters
  integer :: i,ier
  character(len=MAX_STRING_LEN) :: external_source_time_function_filename
  character(len=256) :: line

  ! hardcoded name for now...
  external_source_time_function_filename = 'DATA/stf'

  ! saftey check
  if (USE_LDDRK) then
    print *,'Error: external source time function is not supported yet for LDDRK scheme'
    stop 'Error external source time function is not supported yet for LDDRK scheme'
  endif

  ! Allocate the source time function array to the number of time steps.
  allocate(user_source_time_function(NSTEP),stat=ier)
  if (ier /= 0) stop 'Error allocating external user source time function array'
  user_source_time_function(:) = 0.d0

  ! Read in source time function.
  open(unit=IIN, file=trim(external_source_time_function_filename), &
       status='old', form='formatted', action='read', iostat=ier)
  if (ier /= 0) then
    print *,'Error could not open external source file: ',trim(external_source_time_function_filename)
    stop 'Error opening external source time function file DATA/stf'
  endif

  ! gets number of file entries
  i = 0
  do while (ier == 0)
    read(IIN,"(a256)",iostat=ier) line
    if (ier == 0) then
      ! suppress leading white spaces, if any
      line = adjustl(line)

      ! skip empty/comment lines
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle
      ! stop 'error in format of external_source_time_function_filename, no comments are allowed in it'

      ! increases counter
      i = i + 1
    endif
  enddo
  rewind(IIN)

  ! checks number of lines
  if (i < 1) then
    print *,'Error: External source time function file ',trim(trim(external_source_time_function_filename)),'has no valid data;'
    print *,'       the number of time steps is < 1. Please check the file...'
    stop 'Error: the number of time steps in external_source_time_function_filename is < 1'
  endif

  if (i > NSTEP) then
    print *
    print *,'****************************************************************************************'
    print *,'Warning: ',trim(external_source_time_function_filename),' contains more than NSTEP time steps,'
    print *,'         only the first NSTEP=',NSTEP,' will be read, all the others will be ignored.'
    print *,'****************************************************************************************'
    print *
  endif

  ! checks number of time steps read
  if (i < NSTEP) then
    print *,'Problem when reading external source time file: ', trim(external_source_time_function_filename)
    print *,'  number of time steps in the simulation = ',NSTEP
    print *,'  number of time steps read from the source time function = ',i
    print *,'Please make sure that the number of time steps in the external source file read is greater or &
             &equal to the number of time steps in the simulation'
    stop 'Error invalid number of time steps in external source time file'
  endif

  ! file format: DATA/stf allows for comment lines
  ! # comment
  ! stf_val
  ! stf_val
  ! ...

  ! read the time step used and check that it is the same as DT used for the code
  ier = 0
  i = 0
  do while (ier == 0)
    read(IIN,"(a256)",iostat=ier) line
    if (ier == 0) then
      ! suppress leading white spaces, if any
      line = adjustl(line)

      ! skip empty/comment lines
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle

      ! increases counter
      i = i + 1

      ! gets STF
      if (i <= NSTEP) then
        ! read the source values
        read(line,*,iostat=ier) user_source_time_function(i)
        if (ier /= 0) then
          print *,'Problem when reading external source time file: ', trim(external_source_time_function_filename)
          print *,'  line   : ',trim(line)
          print *,'  counter: ',i
          print *,'Please check, file format should be: '
          print *,'  #(optional)comment '
          print *,'  stf-value'
          print *,'  stf-value'
          print *,'  ..'
          stop 'Error reading external source time file with invalid format'
        endif
      endif

    endif
  enddo

  ! closes external STF file
  close(IIN)

  end subroutine read_external_source_time_function
