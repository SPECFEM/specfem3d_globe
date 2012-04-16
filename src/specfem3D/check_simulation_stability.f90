!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  subroutine check_simulation_stability(it,displ_crust_mantle,displ_inner_core,displ_outer_core, &
                          b_displ_crust_mantle,b_displ_inner_core,b_displ_outer_core, &
                          eps_trace_over_3_crust_mantle, &
                          epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                          epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                          SIMULATION_TYPE,OUTPUT_FILES,time_start,DT,t0,NSTEP, &
                          myrank)

  use constants_solver
  use specfem_par,only: GPU_MODE,Mesh_pointer

  implicit none

  include 'mpif.h'
  include "precision.h"

  ! time step
  integer it,NSTEP,myrank

  ! displacement
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: b_displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: b_displ_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE_ADJOINT) :: b_displ_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY) :: &
    eps_trace_over_3_crust_mantle

!  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT) ::  &
!    epsilondev_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT) ::  &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

  integer SIMULATION_TYPE
  character(len=150) OUTPUT_FILES

  double precision :: time_start,DT,t0

!  logical COMPUTE_AND_STORE_STRAIN

  ! local parameters
  ! maximum of the norm of the displacement and of the potential in the fluid
  real(kind=CUSTOM_REAL) Usolidnorm,Usolidnorm_all,Ufluidnorm,Ufluidnorm_all
  real(kind=CUSTOM_REAL) Strain_norm,Strain_norm_all,Strain2_norm,Strain2_norm_all
  real(kind=CUSTOM_REAL) b_Usolidnorm,b_Usolidnorm_all,b_Ufluidnorm,b_Ufluidnorm_all
  ! names of the data files for all the processors in MPI
  character(len=150) outputname
  ! timer MPI
  double precision :: tCPU,t_remain,t_total
  integer :: ihours,iminutes,iseconds,int_tCPU, &
             ihours_remain,iminutes_remain,iseconds_remain,int_t_remain, &
             ihours_total,iminutes_total,iseconds_total,int_t_total
  ! to determine date and time at which the run will finish
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  character(len=3), dimension(12) :: month_name
  character(len=3), dimension(0:6) :: weekday_name
  data month_name /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
  data weekday_name /'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'/
  integer :: year,mon,day,hr,minutes,timestamp,julian_day_number,day_of_week, &
             timestamp_remote,year_remote,mon_remote,day_remote,hr_remote,minutes_remote,day_of_week_remote
  integer :: ier
  integer, external :: idaywk

  double precision,parameter :: scale_displ = R_EARTH


  ! compute maximum of norm of displacement in each slice
  if( .not. GPU_MODE) then
    ! on CPU
    Usolidnorm = max( &
        maxval(sqrt(displ_crust_mantle(1,:)**2 + &
                    displ_crust_mantle(2,:)**2 + &
                    displ_crust_mantle(3,:)**2)), &
        maxval(sqrt(displ_inner_core(1,:)**2 + &
                    displ_inner_core(2,:)**2 + &
                    displ_inner_core(3,:)**2)))

    Ufluidnorm = maxval(abs(displ_outer_core))

  else
    ! on GPU
    ! way 2: just get maximum of fields from GPU
    call check_norm_elastic_from_device(Usolidnorm,Mesh_pointer,1)
    call check_norm_acoustic_from_device(Ufluidnorm,Mesh_pointer,1)
  endif

  ! check stability of the code, exit if unstable
  ! negative values can occur with some compilers when the unstable value is greater
  ! than the greatest possible floating-point number of the machine
  if(Usolidnorm > STABILITY_THRESHOLD .or. Usolidnorm < 0) &
    call exit_MPI(myrank,'forward simulation became unstable in solid and blew up')
  if(Ufluidnorm > STABILITY_THRESHOLD .or. Ufluidnorm < 0) &
    call exit_MPI(myrank,'forward simulation became unstable in fluid and blew up')

  ! compute the maximum of the maxima for all the slices using an MPI reduction
  call MPI_REDUCE(Usolidnorm,Usolidnorm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
                      MPI_COMM_WORLD,ier)
  call MPI_REDUCE(Ufluidnorm,Ufluidnorm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
                      MPI_COMM_WORLD,ier)

  if (SIMULATION_TYPE == 3) then
    if( .not. GPU_MODE) then
      ! on CPU
      b_Usolidnorm = max( &
             maxval(sqrt(b_displ_crust_mantle(1,:)**2 + &
                          b_displ_crust_mantle(2,:)**2 + b_displ_crust_mantle(3,:)**2)), &
             maxval(sqrt(b_displ_inner_core(1,:)**2  &
                        + b_displ_inner_core(2,:)**2 &
                        + b_displ_inner_core(3,:)**2)))

      b_Ufluidnorm = maxval(abs(b_displ_outer_core))
    else
      ! on GPU
      call check_norm_elastic_from_device(b_Usolidnorm,Mesh_pointer,3)
      call check_norm_acoustic_from_device(b_Ufluidnorm,Mesh_pointer,3)
    endif

    if(b_Usolidnorm > STABILITY_THRESHOLD .or. b_Usolidnorm < 0) &
      call exit_MPI(myrank,'backward simulation became unstable and blew up  in the solid')
    if(b_Ufluidnorm > STABILITY_THRESHOLD .or. b_Ufluidnorm < 0) &
      call exit_MPI(myrank,'backward simulation became unstable and blew up  in the fluid')

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call MPI_REDUCE(b_Usolidnorm,b_Usolidnorm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
             MPI_COMM_WORLD,ier)
    call MPI_REDUCE(b_Ufluidnorm,b_Ufluidnorm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
             MPI_COMM_WORLD,ier)
  endif

  if (COMPUTE_AND_STORE_STRAIN) then
    if( .not. GPU_MODE) then
      ! on CPU
      Strain_norm = maxval(abs(eps_trace_over_3_crust_mantle))
      Strain2_norm= max( maxval(abs(epsilondev_xx_crust_mantle)), &
                         maxval(abs(epsilondev_yy_crust_mantle)), &
                         maxval(abs(epsilondev_xy_crust_mantle)), &
                         maxval(abs(epsilondev_xz_crust_mantle)), &
                         maxval(abs(epsilondev_yz_crust_mantle)) )
    else
      ! on GPU
      call check_norm_strain_from_device(Strain_norm,Strain2_norm,Mesh_pointer)
    endif

    call MPI_REDUCE(Strain_norm,Strain_norm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
             MPI_COMM_WORLD,ier)
    call MPI_REDUCE(Strain2_norm,Strain2_norm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
             MPI_COMM_WORLD,ier)
  endif

  if(myrank == 0) then

    write(IMAIN,*) 'Time step # ',it
    write(IMAIN,*) 'Time: ',sngl(((it-1)*DT-t0)/60.d0),' minutes'

    ! rescale maximum displacement to correct dimensions
    Usolidnorm_all = Usolidnorm_all * sngl(scale_displ)
    write(IMAIN,*) 'Max norm displacement vector U in solid in all slices (m) = ',Usolidnorm_all
    write(IMAIN,*) 'Max non-dimensional potential Ufluid in fluid in all slices = ',Ufluidnorm_all

    if (SIMULATION_TYPE == 3) then
      b_Usolidnorm_all = b_Usolidnorm_all * sngl(scale_displ)
      write(IMAIN,*) 'Max norm displacement vector U in solid in all slices for back prop.(m) = ',b_Usolidnorm_all
      write(IMAIN,*) 'Max non-dimensional potential Ufluid in fluid in all slices for back prop.= ',b_Ufluidnorm_all
    endif

    if(COMPUTE_AND_STORE_STRAIN) then
      write(IMAIN,*) 'Max of strain, eps_trace_over_3_crust_mantle =',Strain_norm_all
      write(IMAIN,*) 'Max of strain, epsilondev_crust_mantle  =',Strain2_norm_all
    endif

    ! elapsed time since beginning of the simulation
    tCPU = MPI_WTIME() - time_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(IMAIN,*) 'Elapsed time in seconds = ',tCPU
    write(IMAIN,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IMAIN,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)

    ! compute estimated remaining simulation time
    t_remain = (NSTEP - it) * (tCPU/dble(it))
    int_t_remain = int(t_remain)
    ihours_remain = int_t_remain / 3600
    iminutes_remain = (int_t_remain - 3600*ihours_remain) / 60
    iseconds_remain = int_t_remain - 3600*ihours_remain - 60*iminutes_remain
    write(IMAIN,*) 'Time steps done = ',it,' out of ',NSTEP
    write(IMAIN,*) 'Time steps remaining = ',NSTEP - it
    write(IMAIN,*) 'Estimated remaining time in seconds = ',t_remain
    write(IMAIN,"(' Estimated remaining time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain

    ! compute estimated total simulation time
    t_total = t_remain + tCPU
    int_t_total = int(t_total)
    ihours_total = int_t_total / 3600
    iminutes_total = (int_t_total - 3600*ihours_total) / 60
    iseconds_total = int_t_total - 3600*ihours_total - 60*iminutes_total
    write(IMAIN,*) 'Estimated total run time in seconds = ',t_total
    write(IMAIN,"(' Estimated total run time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total
    write(IMAIN,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'

    if(it < NSTEP) then

      ! get current date
      call date_and_time(datein,timein,zone,time_values)
      ! time_values(1): year
      ! time_values(2): month of the year
      ! time_values(3): day of the month
      ! time_values(5): hour of the day
      ! time_values(6): minutes of the hour

      ! compute date at which the run should finish; for simplicity only minutes
      ! are considered, seconds are ignored; in any case the prediction is not
      ! accurate down to seconds because of system and network fluctuations
      year = time_values(1)
      mon = time_values(2)
      day = time_values(3)
      hr = time_values(5)
      minutes = time_values(6)

      ! get timestamp in minutes of current date and time
      call convtime(timestamp,year,mon,day,hr,minutes)

      ! add remaining minutes
      timestamp = timestamp + nint(t_remain / 60.d0)

      ! get date and time of that future timestamp in minutes
      call invtime(timestamp,year,mon,day,hr,minutes)

      ! convert to Julian day to get day of the week
      call calndr(day,mon,year,julian_day_number)
      day_of_week = idaywk(julian_day_number)

      write(IMAIN,"(' The run will finish approximately on (in local time): ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
          weekday_name(day_of_week),month_name(mon),day,year,hr,minutes

      ! print date and time estimate of end of run in another country.
      ! For instance: the code runs at Caltech in California but the person
      ! running the code is connected remotely from France, which has 9 hours more
      if(ADD_TIME_ESTIMATE_ELSEWHERE .and. HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE /= 0) then

        ! add time difference with that remote location (can be negative)
        timestamp_remote = timestamp + HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE

        ! get date and time of that future timestamp in minutes
        call invtime(timestamp_remote,year_remote,mon_remote,day_remote,hr_remote,minutes_remote)

        ! convert to Julian day to get day of the week
        call calndr(day_remote,mon_remote,year_remote,julian_day_number)
        day_of_week_remote = idaywk(julian_day_number)

        if(HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE > 0) then
          write(IMAIN,*) 'Adding positive time difference of ',abs(HOURS_TIME_DIFFERENCE),' hours'
        else
          write(IMAIN,*) 'Adding negative time difference of ',abs(HOURS_TIME_DIFFERENCE),' hours'
        endif
        write(IMAIN,*) 'and ',abs(MINUTES_TIME_DIFFERENCE),' minutes to get estimate at a remote location'
        write(IMAIN, &
            "(' The run will finish approximately on: ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
            weekday_name(day_of_week_remote),month_name(mon_remote),day_remote,year_remote,hr_remote,minutes_remote
      endif

      if(it < 100) then
        write(IMAIN,*) '************************************************************'
        write(IMAIN,*) '**** BEWARE: the above time estimates are not reliable'
        write(IMAIN,*) '**** because fewer than 100 iterations have been performed'
        write(IMAIN,*) '************************************************************'
      endif

    endif

    write(IMAIN,*)

    ! write time stamp file to give information about progression of simulation
    write(outputname,"('/timestamp',i6.6)") it

    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write')

    write(IOUT,*) 'Time step # ',it
    write(IOUT,*) 'Time: ',sngl(((it-1)*DT-t0)/60.d0),' minutes'
    write(IOUT,*)
    write(IOUT,*) 'Max norm displacement vector U in solid in all slices (m) = ',Usolidnorm_all
    write(IOUT,*) 'Max non-dimensional potential Ufluid in fluid in all slices = ',Ufluidnorm_all
    write(IOUT,*)

    if (SIMULATION_TYPE == 3) then
      write(IOUT,*) 'Max norm displacement vector U in solid in all slices for back prop.(m) = ',b_Usolidnorm_all
      write(IOUT,*) 'Max non-dimensional potential Ufluid in fluid in all slices for back prop.= ',b_Ufluidnorm_all
      write(IOUT,*)
    endif

    write(IOUT,*) 'Elapsed time in seconds = ',tCPU
    write(IOUT,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    write(IOUT,*)

    write(IOUT,*) 'Time steps done = ',it,' out of ',NSTEP
    write(IOUT,*) 'Time steps remaining = ',NSTEP - it
    write(IOUT,*) 'Estimated remaining time in seconds = ',t_remain
    write(IOUT,"(' Estimated remaining time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain
    write(IOUT,*)

    write(IOUT,*) 'Estimated total run time in seconds = ',t_total
    write(IOUT,"(' Estimated total run time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total
    write(IOUT,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'
    write(IOUT,*)

    if(it < NSTEP) then

      write(IOUT,"(' The run will finish approximately on (in local time): ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
          weekday_name(day_of_week),month_name(mon),day,year,hr,minutes

      ! print date and time estimate of end of run in another country.
      ! For instance: the code runs at Caltech in California but the person
      ! running the code is connected remotely from France, which has 9 hours more
      if(ADD_TIME_ESTIMATE_ELSEWHERE .and. HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE /= 0) then
        if(HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE > 0) then
          write(IOUT,*) 'Adding positive time difference of ',abs(HOURS_TIME_DIFFERENCE),' hours'
        else
          write(IOUT,*) 'Adding negative time difference of ',abs(HOURS_TIME_DIFFERENCE),' hours'
        endif
        write(IOUT,*) 'and ',abs(MINUTES_TIME_DIFFERENCE),' minutes to get estimate at a remote location'
        write(IOUT, &
            "(' The run will finish approximately on (in remote time): ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
            weekday_name(day_of_week_remote),month_name(mon_remote), &
            day_remote,year_remote,hr_remote,minutes_remote
      endif

      if(it < 100) then
        write(IOUT,*)
        write(IOUT,*) '************************************************************'
        write(IOUT,*) '**** BEWARE: the above time estimates are not reliable'
        write(IOUT,*) '**** because fewer than 100 iterations have been performed'
        write(IOUT,*) '************************************************************'
      endif

    endif

    close(IOUT)

    ! check stability of the code, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if(Usolidnorm_all > STABILITY_THRESHOLD .or. Usolidnorm_all < 0) &
      call exit_MPI(myrank,'forward simulation became unstable and blew up in the solid')
    if(Ufluidnorm_all > STABILITY_THRESHOLD .or. Ufluidnorm_all < 0) &
      call exit_MPI(myrank,'forward simulation became unstable and blew up in the fluid')

    if(SIMULATION_TYPE == 3) then
      if(b_Usolidnorm_all > STABILITY_THRESHOLD .or. b_Usolidnorm_all < 0) &
        call exit_MPI(myrank,'backward simulation became unstable and blew up in the solid')
      if(b_Ufluidnorm_all > STABILITY_THRESHOLD .or. b_Ufluidnorm_all < 0) &
        call exit_MPI(myrank,'backward simulation became unstable and blew up in the fluid')
    endif

  endif

  end subroutine check_simulation_stability
