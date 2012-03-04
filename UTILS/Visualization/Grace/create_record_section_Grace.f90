
  program create_record_section_Grace

! creates files to plot record sections with "Grace"
! then call with " xmgrace *.sem.recordsection "
! call this program with " xcreate_record_section_Grace < create_record_section_Grace.in "

! Authors: Jeroen Tromp and Dimitri Komatitsch, mid-2000 and September 2007

  implicit none

! total number of time steps in each seismogram
  integer, parameter :: NSTEP = 41030

! total number of seismic stations (receivers)
  integer, parameter :: NSTATIONS = 180

  integer :: it,irec

  real :: epi,epi_start,epi_end,max_amp,t_start,t_end,scale,enhance

  real, dimension(NSTEP) :: time_sem,sem

  character(len=3) component
  character(len=150) station

  print *,'component (LHZ, LHE, or LHN): '
  read *,component

  print *,'starting epicentral distance (degrees):'
  read *,epi_start

  print *,'ending epicentral distance (degrees):'
  read *,epi_end

  print *,'starting time in record section (minutes):'
  read *,t_start

  print *,'ending in record section (minutes):'
  read *,t_end

  print *,'enhancement factor:'
  read *,enhance

! open the station file
  open(unit=1,file='output_list_stations.txt',status='unknown')

  do irec = 1,NSTATIONS

  read(1,*) station,epi

  print *,'processing file ',irec,' out of ',NSTATIONS

  if(epi >= epi_start .and. epi <= epi_end) then

    max_amp = -100000.

! read a given seismogram
    open(unit=10,file=station(1:len_trim(station))//'.'//component(1:3)//'.semd.ascii.convolved',status='old')
    do it = 1,NSTEP
      read(10,*) time_sem(it),sem(it)
! compute maximum of all the traces
      if(time_sem(it)/60.0 >= t_start .and. time_sem(it)/60.0 <= t_end .and. abs(sem(it)) > max_amp) max_amp = abs(sem(it))
    enddo
    close(10)

! compute scaling factor to plot the traces
    scale = enhance/max_amp

! open the record section file
    open(unit=10,file=station(1:len_trim(station))//'.'//component(1:3)//'.sem.recordsection',status='unknown')
    do it = 1,NSTEP
      if(time_sem(it)/60.0 >= t_start .and. time_sem(it)/60.0 <= t_end) write(10,*) time_sem(it)/60.0,epi+sem(it)*scale
    enddo
    close(10)

  endif

  enddo

  close(1)

  end program create_record_section_Grace

