!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! write seismograms to text files

  subroutine write_seismograms(myrank,seismograms,number_receiver_global, &
               station_name,network_name,nrec,nrec_local, &
               DT,NSTEP,hdur,LOCAL_PATH,it_begin,it_end)

  implicit none

  include "constants.h"

  integer nrec,nrec_local,NSTEP,myrank,it_begin,it_end
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms
  double precision hdur,DT
  character(len=150) LOCAL_PATH

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer irec,irec_local,length_station_name,length_network_name
  integer iorientation,isample
  double precision value

  character(len=4) chn
  character(len=150) sisname,clean_LOCAL_PATH,final_LOCAL_PATH

  do irec_local = 1,nrec_local

! get global number of that receiver
    irec = number_receiver_global(irec_local)

    do iorientation = 1,3

      if(iorientation == 1) then
        chn = 'LHN'
      else if(iorientation == 2) then
        chn = 'LHE'
      else if(iorientation == 3) then
        chn = 'LHZ'
      else
        call exit_MPI(myrank,'incorrect channel value')
      endif

! create the name of the seismogram file for each slice
! file name includes the name of the station and the network
      length_station_name = len_trim(station_name(irec))
      length_network_name = len_trim(network_name(irec))

! check that length conforms to standard
      if(length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) &
           call exit_MPI(myrank,'wrong length of station name')

      if(length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) &
           call exit_MPI(myrank,'wrong length of network name')

      write(sisname,"(a,'.',a,'.',a3,'.semd')") station_name(irec)(1:length_station_name), &
                    network_name(irec)(1:length_network_name),chn

! suppress white spaces if any
    clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full final local path
    final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
    open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),status='unknown')

! subtract half duration of the source to make sure travel time is correct
      do isample = it_begin,it_end
        value = dble(seismograms(iorientation,irec_local,isample))
! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',sngl(value)
        else
          write(IOUT,*) dble(isample-1)*DT - hdur,' ',value
        endif
      enddo

      close(IOUT)

      enddo

  enddo

  end subroutine write_seismograms

!=====================================================================

! write adjoint seismograms to text files

  subroutine write_adj_seismograms(myrank,seismograms,number_receiver_global, &
               nrec_local,it,DT,NSTEP,hdur,LOCAL_PATH)

  implicit none

  include "constants.h"

  integer nrec_local,NSTEP,it,myrank
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(6,nrec_local,NSTEP) :: seismograms
  double precision hdur,DT
  character(len=150) LOCAL_PATH

  integer irec,irec_local
  integer iorientation,isample

  character(len=4) chn
  character(len=150) sisname,clean_LOCAL_PATH,final_LOCAL_PATH

  do irec_local = 1,nrec_local

! get global number of that receiver
    irec = number_receiver_global(irec_local)

    do iorientation = 1,6

      if(iorientation == 1) then
        chn = 'SNN'
      else if(iorientation == 2) then
        chn = 'SEE'
      else if(iorientation == 3) then
        chn = 'SZZ'
      else if(iorientation == 4) then
        chn = 'SNE'
      else if(iorientation == 5) then
        chn = 'SNZ'
      else if(iorientation == 6) then
        chn = 'SEZ'
      endif

! create the name of the seismogram file for each slice
! file name includes the name of the station, the network and the component

      write(sisname,"(a,i3.3,'.',a,'.',a3,'.sem')") 'S',irec,&
           'NT',chn

! suppress white spaces if any
    clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full final local path
    final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
      open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),status='unknown')
! make sure we never write more than the maximum number of time steps
! subtract half duration of the source to make sure travel time is correct
      do isample = 1,min(it,NSTEP)
! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',seismograms(iorientation,irec_local,isample)
        else
          write(IOUT,*) dble(isample-1)*DT - hdur,' ',seismograms(iorientation,irec_local,isample)
        endif
      enddo

      close(IOUT)

      enddo

  enddo

  end subroutine write_adj_seismograms
