!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology August 2003
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
  real(kind=CUSTOM_REAL), dimension(3,nrec_local,NSTEP) :: seismograms
  character(len=8), dimension(nrec) :: station_name,network_name
  double precision hdur,DT
  character(len=150) LOCAL_PATH

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
! file name includes the name of the station, the network and the component
      length_station_name = len_trim(station_name(irec))
      length_network_name = len_trim(network_name(irec))

! check that length conforms to standard
      if(length_station_name < 1 .or. length_station_name > 8 .or. &
         length_network_name < 1 .or. length_network_name > 2) &
           call exit_MPI(myrank,'wrong length of station or network name')

    if(length_network_name == 1) then

      if(length_station_name == 1) then
        write(sisname,"(a1,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 2) then
        write(sisname,"(a2,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 3) then
        write(sisname,"(a3,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 4) then
        write(sisname,"(a4,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 5) then
        write(sisname,"(a5,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 6) then
        write(sisname,"(a6,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 7) then
        write(sisname,"(a7,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else
        write(sisname,"(a8,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      endif

    else

      if(length_station_name == 1) then
        write(sisname,"(a1,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 2) then
        write(sisname,"(a2,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 3) then
        write(sisname,"(a3,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 4) then
        write(sisname,"(a4,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 5) then
        write(sisname,"(a5,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 6) then
        write(sisname,"(a6,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 7) then
        write(sisname,"(a7,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else
        write(sisname,"(a8,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      endif

    endif

! suppress white spaces if any
    clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full final local path
    final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
      open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname,status='unknown')

! subtract half duration of the source to make sure travel time is correct
      do isample = it_begin,it_end
        value = dble(seismograms(iorientation,irec_local,isample))
        write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',sngl(value)
      enddo

      close(IOUT)

      enddo

  enddo

  end subroutine write_seismograms

