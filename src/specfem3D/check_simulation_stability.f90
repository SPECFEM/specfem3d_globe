!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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
!ZN                          eps_trace_over_3_crust_mantle,epsilondev_crust_mantle, &
                          SIMULATION_TYPE,OUTPUT_FILES,time_start,DT,t0,NSTEP, &
                          it_begin,it_end,NUMBER_OF_THIS_RUN,NUMBER_OF_RUNS,myrank)

  implicit none

  include 'mpif.h'
  include "constants.h"
  include "precision.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  ! time step
  integer it,it_begin,it_end,NUMBER_OF_THIS_RUN,NUMBER_OF_RUNS,NSTEP,myrank

  ! displacement
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: displ_outer_core

!ZN  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY) :: &
!ZN    eps_trace_over_3_crust_mantle
!ZN  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT) ::  &
!ZN    epsilondev_crust_mantle

  integer SIMULATION_TYPE
  character(len=150) OUTPUT_FILES

  double precision :: time_start,DT,t0

  ! local parameters
  ! maximum of the norm of the displacement and of the potential in the fluid
  real(kind=CUSTOM_REAL) Usolidnorm,Usolidnorm_all,Ufluidnorm,Ufluidnorm_all
!ZN  real(kind=CUSTOM_REAL) Strain_norm,Strain_norm_all,strain2_norm,strain2_norm_all
  ! names of the data files for all the processors in MPI
  character(len=150) outputname
  ! timer MPI
  double precision :: tCPU,t_remain,t_total,t_remain_run,t_total_run
  integer :: ihours,iminutes,iseconds,int_tCPU, &
             ihours_remain,iminutes_remain,iseconds_remain,int_t_remain, &
             ihours_total,iminutes_total,iseconds_total,int_t_total
  integer :: it_run,nstep_run, &
             ihours_remain_run,iminutes_remain_run,iseconds_remain_run,int_t_remain_run, &
             ihours_total_run,iminutes_total_run,iseconds_total_run,int_t_total_run
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

  logical :: SHOW_SEPARATE_RUN_INFORMATION

  ! compute maximum of norm of displacement in each slice
  Usolidnorm = max( &
      maxval(sqrt(displ_crust_mantle(1,:)**2 + &
                  displ_crust_mantle(2,:)**2 + displ_crust_mantle(3,:)**2)), &
      maxval(sqrt(displ_inner_core(1,:)**2 + displ_inner_core(2,:)**2 + displ_inner_core(3,:)**2)))

  Ufluidnorm = maxval(abs(displ_outer_core))

  ! compute the maximum of the maxima for all the slices using an MPI reduction
  call MPI_REDUCE(Usolidnorm,Usolidnorm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
                      MPI_COMM_WORLD,ier)
  call MPI_REDUCE(Ufluidnorm,Ufluidnorm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
                      MPI_COMM_WORLD,ier)

!ZN  if (COMPUTE_AND_STORE_STRAIN) then
!ZN    Strain_norm = maxval(abs(eps_trace_over_3_crust_mantle))
!ZN    strain2_norm= maxval(abs(epsilondev_crust_mantle))
!ZN    call MPI_REDUCE(Strain_norm,Strain_norm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
!ZN             MPI_COMM_WORLD,ier)
!ZN    call MPI_REDUCE(Strain2_norm,Strain2_norm_all,1,CUSTOM_MPI_TYPE,MPI_MAX,0, &
!ZN             MPI_COMM_WORLD,ier)
!ZN  endif

  if(myrank == 0) then

    write(IMAIN,*) 'Time step # ',it
    write(IMAIN,*) 'Time: ',sngl(((it-1)*DT-t0)/60.d0),' minutes'

    ! rescale maximum displacement to correct dimensions
    Usolidnorm_all = Usolidnorm_all * sngl(scale_displ)
    if (SIMULATION_TYPE == 1) then
      write(IMAIN,*) 'Max norm displacement vector U in solid in all slices (m) = ',Usolidnorm_all
      write(IMAIN,*) 'Max non-dimensional potential Ufluid in fluid in all slices = ',Ufluidnorm_all
    else
      write(IMAIN,*) 'Max norm displacement vector U in solid in all slices for back prop.(m) = ',Usolidnorm_all
      write(IMAIN,*) 'Max non-dimensional potential Ufluid in fluid in all slices for back prop.= ',Ufluidnorm_all
    endif

!! DK DK UNDO_ATT
!   if(COMPUTE_AND_STORE_STRAIN) then
!ZN    if(SIMULATION_TYPE == 1 .and. COMPUTE_AND_STORE_STRAIN) then
!ZN      write(IMAIN,*) 'Max of strain, eps_trace_over_3_crust_mantle =',Strain_norm_all
!ZN      write(IMAIN,*) 'Max of strain, epsilondev_crust_mantle  =',Strain2_norm_all
!ZN    endif

    ! information about the current run only
    SHOW_SEPARATE_RUN_INFORMATION = NUMBER_OF_RUNS > 1 .and. NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS
    it_run = it - it_begin + 1
    nstep_run = it_end - it_begin + 1

    ! elapsed time since beginning of the simulation
    tCPU = MPI_WTIME() - time_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(IMAIN,*) 'Elapsed time in seconds = ',tCPU
    write(IMAIN,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IMAIN,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it_run)

    ! compute estimated remaining simulation time
    t_remain = (NSTEP - it) * (tCPU/dble(it_run))
    int_t_remain = int(t_remain)
    ihours_remain = int_t_remain / 3600
    iminutes_remain = (int_t_remain - 3600*ihours_remain) / 60
    iseconds_remain = int_t_remain - 3600*ihours_remain - 60*iminutes_remain

    ! compute estimated total simulation time
    t_total = t_remain + tCPU
    int_t_total = int(t_total)
    ihours_total = int_t_total / 3600
    iminutes_total = (int_t_total - 3600*ihours_total) / 60
    iseconds_total = int_t_total - 3600*ihours_total - 60*iminutes_total

    ! calculate times for the *current* run
    if (SHOW_SEPARATE_RUN_INFORMATION) then
      ! compute estimated remaining simulation time
      t_remain_run = (it_end - it) * (tCPU/dble(it_run))
      int_t_remain_run = int(t_remain_run)
      ihours_remain_run = int_t_remain_run / 3600
      iminutes_remain_run = (int_t_remain_run - 3600*ihours_remain_run) / 60
      iseconds_remain_run = int_t_remain_run - 3600*ihours_remain_run - 60*iminutes_remain_run

      ! compute estimated total simulation time
      t_total_run = t_remain_run + tCPU
      int_t_total_run = int(t_total_run)
      ihours_total_run = int_t_total_run / 3600
      iminutes_total_run = (int_t_total_run - 3600*ihours_total_run) / 60
      iseconds_total_run = int_t_total_run - 3600*ihours_total_run - 60*iminutes_total_run
    endif

    ! print time information
    if (SHOW_SEPARATE_RUN_INFORMATION) then
      write(IMAIN,*) 'Time steps done for this run = ',it_run,' out of ',nstep_run
      write(IMAIN,*) 'Time steps done in total = ',it,' out of ',NSTEP
      write(IMAIN,*) 'Time steps remaining for this run = ',it_end - it
      write(IMAIN,*) 'Time steps remaining for all runs = ',NSTEP - it
      write(IMAIN,*) 'Estimated remaining time for this run in seconds = ',t_remain_run
      write(IMAIN,"(' Estimated remaining time for this run in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_remain_run,iminutes_remain_run,iseconds_remain_run
      write(IMAIN,*) 'Estimated remaining time for all runs in seconds = ',t_remain
      write(IMAIN,"(' Estimated remaining time for all runs in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_remain,iminutes_remain,iseconds_remain
      write(IMAIN,*) 'Estimated total run time for this run in seconds = ',t_total_run
      write(IMAIN,"(' Estimated total run time for this run in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_total_run,iminutes_total_run,iseconds_total_run
      write(IMAIN,*) 'We have done ',sngl(100.d0*dble(it_run)/dble(nstep_run)),'% of this run'
      write(IMAIN,*) 'Estimated total run time for all runs in seconds = ',t_total
      write(IMAIN,"(' Estimated total run time for all runs in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_total,iminutes_total,iseconds_total
      write(IMAIN,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of the total'
    else
      write(IMAIN,*) 'Time steps done = ',it,' out of ',NSTEP
      write(IMAIN,*) 'Time steps remaining = ',NSTEP - it
      write(IMAIN,*) 'Estimated total remaining time in seconds = ',t_remain
      write(IMAIN,"(' Estimated total remaining time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_remain,iminutes_remain,iseconds_remain
      write(IMAIN,*) 'Estimated total run time in seconds = ',t_total
      write(IMAIN,"(' Estimated total run time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_total,iminutes_total,iseconds_total
      write(IMAIN,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'
    endif

    if (it < it_end) then

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
      if (SHOW_SEPARATE_RUN_INFORMATION) then
        timestamp = timestamp + nint(t_remain_run / 60.d0)
      else
        timestamp = timestamp + nint(t_remain / 60.d0)
      endif

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

      if (it_run < 100) then
        write(IMAIN,*) '************************************************************'
        write(IMAIN,*) '**** BEWARE: the above time estimates are not reliable'
        write(IMAIN,*) '**** because fewer than 100 iterations have been performed'
        write(IMAIN,*) '************************************************************'
      endif

    endif

    write(IMAIN,*)

    ! write time stamp file to give information about progression of simulation
!! DK DK UNDO_ATT
    if(SIMULATION_TYPE == 1) then
!     write(outputname,"('/timestamp',i6.6)") it
      write(outputname,"('/timestamp_forward',i6.6)") it
    else
      write(outputname,"('/timestamp_backward',i6.6)") it
    endif

    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write')

    write(IOUT,*) 'Time step # ',it
    write(IOUT,*) 'Time: ',sngl(((it-1)*DT-t0)/60.d0),' minutes'
    write(IOUT,*)
    if (SIMULATION_TYPE == 1) then
      write(IOUT,*) 'Max norm displacement vector U in solid in all slices (m) = ',Usolidnorm_all
      write(IOUT,*) 'Max non-dimensional potential Ufluid in fluid in all slices = ',Ufluidnorm_all
    else
      write(IOUT,*) 'Max norm displacement vector U in solid in all slices for back prop.(m) = ',Usolidnorm_all
      write(IOUT,*) 'Max non-dimensional potential Ufluid in fluid in all slices for back prop.= ',Ufluidnorm_all
    endif
    write(IOUT,*)

    write(IOUT,*) 'Elapsed time in seconds = ',tCPU
    write(IOUT,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it_run)
    write(IOUT,*)

    if (SHOW_SEPARATE_RUN_INFORMATION) then

      write(IOUT,*) 'Time steps done for this run = ',it_run,' out of ',nstep_run
      write(IOUT,*) 'Time steps done in total = ',it,' out of ',NSTEP
      write(IOUT,*) 'Time steps remaining for this run = ',it_end - it
      write(IOUT,*) 'Time steps remaining for all runs = ',NSTEP - it
      write(IOUT,*) 'Estimated remaining time for this run in seconds = ',t_remain_run
      write(IOUT,"(' Estimated remaining time for this run in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_remain_run,iminutes_remain_run,iseconds_remain_run
      write(IOUT,*) 'Estimated remaining time for all runs in seconds = ',t_remain
      write(IOUT,"(' Estimated remaining time for all runs in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_remain,iminutes_remain,iseconds_remain
      write(IOUT,*)

      write(IOUT,*) 'Estimated total run time for this run in seconds = ',t_total_run
      write(IOUT,"(' Estimated total run time for this run in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_total_run,iminutes_total_run,iseconds_total_run
      write(IOUT,*) 'We have done ',sngl(100.d0*dble(it_run)/dble(nstep_run)),'% of this run'
      write(IOUT,*) 'Estimated total run time for all runs in seconds = ',t_total
      write(IOUT,"(' Estimated total run time for all runs in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
               ihours_total,iminutes_total,iseconds_total
      write(IOUT,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of the total'
      write(IOUT,*)

    else

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

    endif

    if (it < it_end) then

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

      if (it_run < 100) then
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
    if(Usolidnorm_all > STABILITY_THRESHOLD .or. Usolidnorm_all < 0) then
      if(SIMULATION_TYPE == 1) then
        call exit_MPI(myrank,'forward simulation became unstable and blew up in the solid')
      else
        call exit_MPI(myrank,'backward simulation became unstable and blew up in the solid')
      endif
    endif

    if(Ufluidnorm_all > STABILITY_THRESHOLD .or. Ufluidnorm_all < 0) then
      if(SIMULATION_TYPE == 1) then
        call exit_MPI(myrank,'forward simulation became unstable and blew up in the fluid')
      else
        call exit_MPI(myrank,'backward simulation became unstable and blew up in the fluid')
      endif
    endif

  endif

  end subroutine check_simulation_stability
