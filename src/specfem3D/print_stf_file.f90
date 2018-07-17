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


  subroutine print_stf_file()

! prints source time function and spectrum

  use constants
  use shared_input_parameters

  use specfem_par, only: NSOURCES,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
    factor_force_source, &
    NSTEP,DT,t0,tshift_src,hdur

  implicit none

  ! local parameters
  integer :: isource,it,iom,ier
  double precision :: timeval,time_t
  double precision :: stf
  double precision :: strength,scaleF,om

  character(len=MAX_STRING_LEN) :: plot_file

  double precision, external :: get_cmt_scalar_moment
  double precision, external :: get_stf_viscoelastic
  double precision, external :: comp_source_spectrum

  ! number of points to plot the source time function and spectrum
  integer, parameter :: NSAMP_PLOT_SOURCE = 1000

  ! only master process outputs stf
  if (myrank /= 0) return

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'printing the source-time function'
  call flush_IMAIN()

  ! source time function
  do isource = 1,NSOURCES

    ! print the source-time function
    if (NSOURCES == 1) then
      plot_file = '/plot_source_time_function.txt'
    else
      if (isource < 10) then
        write(plot_file,"('/plot_source_time_function',i1,'.txt')") isource
      else if (isource < 100) then
        write(plot_file,"('/plot_source_time_function',i2,'.txt')") isource
      else
        write(plot_file,"('/plot_source_time_function',i3,'.txt')") isource
      endif
    endif

    ! output file
    open(unit=IOUT,file=trim(OUTPUT_FILES)//plot_file,status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_mpi(0,'Error opening plot_source_time_function file')

    ! source strength
    if (USE_FORCE_POINT_SOURCE) then
      ! scale and non-dimensionalize the factor_force_source
      ! factor_force_source in FORCESOLUTION file is in Newton
      ! 1 Newton is 1 kg * 1 m / (1 second)^2
      scaleF = RHOAV * (R_EARTH**4) * PI*GRAV*RHOAV
      strength = factor_force_source(isource) * scaleF
    else
      ! calculates scalar moment M0
      strength = get_cmt_scalar_moment(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource))
    endif

    ! file header
    write(IOUT,*) "# Source Time function"
    if (USE_FORCE_POINT_SOURCE) then
      write(IOUT,*) "# source type point force"
      write(IOUT,*) "# scaling factor force = ",sngl(strength),"(N)"
      write(IOUT,*) "# format : #time(s)   #stf   #factor(Newton)"
    else
      write(IOUT,*) "# source type CMT"
      write(IOUT,*) "# scaling scalar moment (M0) = ",sngl(strength)
      write(IOUT,*) "# format: #time(s)   #stf   #scalar_moment"
    endif

    ! writes out source time function to file
    do it = 1,NSTEP
      ! time for simulation
      time_t = dble(it-1)*DT - t0

      ! sets current time for this source
      timeval = time_t - tshift_src(isource)

      ! source time function value (in range [-1,1]
      stf = get_stf_viscoelastic(timeval,isource)

      ! format: #time(s)   #stf   #scaling (M0 or force factor)
      write(IOUT,*) time_t, sngl(stf), sngl(strength)

    enddo
    close(IOUT)

  enddo ! NSOURCES

  ! source spectra
  write(IMAIN,*)
  write(IMAIN,*) 'printing the source spectrum'
  call flush_IMAIN()

  do isource = 1,NSOURCES
    ! print the spectrum of the derivative of the source from 0 to 1/8 Hz
    if (NSOURCES == 1) then
      plot_file = '/plot_source_spectrum.txt'
    else
     if (isource < 10) then
        write(plot_file,"('/plot_source_spectrum',i1,'.txt')") isource
      else if (isource < 100) then
        write(plot_file,"('/plot_source_spectrum',i2,'.txt')") isource
      else
        write(plot_file,"('/plot_source_spectrum',i3,'.txt')") isource
      endif
    endif

    open(unit=IOUT,file=trim(OUTPUT_FILES)//plot_file,status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_mpi(0,'Error opening plot_source_spectrum file')

    do iom = 1,NSAMP_PLOT_SOURCE
      om = TWO_PI*(1.0d0/8.0d0)*(iom-1)/dble(NSAMP_PLOT_SOURCE-1)
      write(IOUT,*) sngl(om/TWO_PI),sngl(strength * om * comp_source_spectrum(om,hdur(isource)))
    enddo
    close(IOUT)

  enddo ! NSOURCES

  end subroutine print_stf_file
