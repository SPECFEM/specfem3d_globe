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

  double precision function comp_source_time_function(t,hdur)

  use constants, only: EXTERNAL_SOURCE_TIME_FUNCTION

  implicit none

  double precision,intent(in) :: t,hdur

  ! local parameters
  double precision, external :: comp_source_time_function_heavi
  double precision, external :: comp_source_time_function_ext

  if (EXTERNAL_SOURCE_TIME_FUNCTION) then
    ! external stf
    comp_source_time_function = comp_source_time_function_ext()
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
  comp_source_time_function_heavi = 0.5d0*(1.0d0 + netlib_specfun_erf(t/hdur))

  end function comp_source_time_function_heavi


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

  double precision function comp_source_time_function_ext()

  use specfem_par, only: it, stfArray_external
  implicit none

  ! On the first iteration, go get the ASCII file.
  if (.not. allocated (stfArray_external)) then
    call get_EXTERNAL_SOURCE_TIME_FUNCTION()
  endif

  comp_source_time_function_ext = stfArray_external(it)

  end function comp_source_time_function_ext

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_EXTERNAL_SOURCE_TIME_FUNCTION()

  use specfem_par, only: NSTEP, stfArray_external, IIN
  implicit none

  integer :: iterator,ier
  character(len=256) :: line

  ! Allocate the source time function array to the number of time steps.
  allocate( stfArray_external(NSTEP),stat=ier)
  if (ier /= 0 ) stop 'Error allocating external source time function array'

  print *, NSTEP

  ! Read in source time function.
  open(unit=IIN, file='DATA/stf', status='old', form='formatted',iostat=ier)
  if (ier /= 0 ) stop 'Error opening file DATA/stf'

  read_loop: do iterator=1,NSTEP

    read(IIN, '(A)', iostat = ier) line

    if (ier /= 0) then
      print *, "Error in external source time function."
      stop 'Error reading external stf file'
    endif

    ! Ignore lines with a hash (comments)
    if (index(line, "#") /= 0) cycle read_loop

    read(line, *) stfArray_external(iterator)

  enddo read_loop

  close(IIN)

  end subroutine get_EXTERNAL_SOURCE_TIME_FUNCTION
