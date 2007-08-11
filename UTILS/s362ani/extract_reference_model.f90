
! convert reference model file into Fortran statements

! Author: Dimitri Komatitsch, University of Pau, France, August 2007

  program extract_reference_model

  implicit none

  include "../../constants.h"

! number of arrays to define
  integer, parameter :: NARRAYS = 9

  integer :: i,i1,i2,iarray

  double precision, dimension(NR_REF) :: radius_ref,density_ref,vpv_ref,vph_ref,vsv_ref,vsh_ref,eta_ref,Qkappa_ref,Qmu_ref

  if(NR_REF /= 750) stop 'NR_REF should be equal to 750 in the Harvard model'

  open(unit=10,file='./REF',status='old',action='read')

! skip the first three lines
  read(10,*)
  read(10,*)
  read(10,*)

! read the card deck model
  do i=1,NR_REF
    read(10,*) radius_ref(i),density_ref(i),vpv_ref(i),vsv_ref(i),Qkappa_ref(i),Qmu_ref(i),vph_ref(i),vsh_ref(i),eta_ref(i)
  enddo

  close(10)

! write the model in array constructor Fortran90 format
  do iarray = 1,NARRAYS

  i = 0

  do i1 = 1,25

    if(iarray == 1) then
      print *,'Mref_V%radius_ref(',i+1,':',i+30,') = (/ &'
    else if(iarray == 2) then
      print *,'Mref_V%density_ref(',i+1,':',i+30,') = (/ &'
    else if(iarray == 3) then
      print *,'Mref_V%vpv_ref(',i+1,':',i+30,') = (/ &'
    else if(iarray == 4) then
      print *,'Mref_V%vsv_ref(',i+1,':',i+30,') = (/ &'
    else if(iarray == 5) then
      print *,'Mref_V%Qkappa_ref(',i+1,':',i+30,') = (/ &'
    else if(iarray == 6) then
      print *,'Mref_V%Qmu_ref(',i+1,':',i+30,') = (/ &'
    else if(iarray == 7) then
      print *,'Mref_V%vph_ref(',i+1,':',i+30,') = (/ &'
    else if(iarray == 8) then
      print *,'Mref_V%vsh_ref(',i+1,':',i+30,') = (/ &'
    else if(iarray == 9) then
      print *,'Mref_V%eta_ref(',i+1,':',i+30,') = (/ &'
    else
      stop 'incorrect array number'
    endif

! define groups of 30 lines for each constructor
  do i2 = 1,30

    i = i + 1

! continue to next line or close constructor if this is the last line
    if(i2 /= 30) then

    if(iarray == 1) then
      print *,radius_ref(i),', &'
    else if(iarray == 2) then
      print *,density_ref(i),', &'
    else if(iarray == 3) then
      print *,vpv_ref(i),', &'
    else if(iarray == 4) then
      print *,vsv_ref(i),', &'
    else if(iarray == 5) then
      print *,Qkappa_ref(i),', &'
    else if(iarray == 6) then
      print *,Qmu_ref(i),', &'
    else if(iarray == 7) then
      print *,vph_ref(i),', &'
    else if(iarray == 8) then
      print *,vsh_ref(i),', &'
    else if(iarray == 9) then
      print *,eta_ref(i),', &'
    else
      stop 'incorrect array number'
    endif

    else

    if(iarray == 1) then
      print *,radius_ref(i),' /)'
    else if(iarray == 2) then
      print *,density_ref(i),' /)'
    else if(iarray == 3) then
      print *,vpv_ref(i),' /)'
    else if(iarray == 4) then
      print *,vsv_ref(i),' /)'
    else if(iarray == 5) then
      print *,Qkappa_ref(i),' /)'
    else if(iarray == 6) then
      print *,Qmu_ref(i),' /)'
    else if(iarray == 7) then
      print *,vph_ref(i),' /)'
    else if(iarray == 8) then
      print *,vsh_ref(i),' /)'
    else if(iarray == 9) then
      print *,eta_ref(i),' /)'
    else
      stop 'incorrect array number'
    endif

    endif

  enddo

    print *

  enddo

  enddo

  end program extract_reference_model

