
! convert reference model file into Fortran statements

! Author: Dimitri Komatitsch, University of Pau, France, August 2007

  program extract_model_reference

  implicit none

  include "constants.h"

  integer :: i

  double precision :: radius_ref
  double precision :: density_ref
  double precision :: vpv_ref
  double precision :: vph_ref
  double precision :: vsv_ref
  double precision :: vsh_ref
  double precision :: eta_ref
  double precision :: Qkappa_ref
  double precision :: Qmu_ref

  open(unit=10,file='DATA/s362ani/REF',status='old',action='read')
! skip the first three lines
  read(10,*)
  read(10,*)
  read(10,*)
! read the card deck model
  do i=1,NR_REF

    read(10,*) radius_ref,density_ref,vpv_ref, &
               vsv_ref,Qkappa_ref,Qmu_ref, &
               vph_ref,vsh_ref,eta_ref

    print *,'Mref_V%radius_ref(',i,') = ',radius_ref
    print *,'Mref_V%density_ref(',i,') = ',density_ref
    print *,'Mref_V%vpv_ref(',i,') = ',vpv_ref
    print *,'Mref_V%vsv_ref(',i,') = ',vsv_ref
    print *,'Mref_V%Qkappa_ref(',i,') = ',Qkappa_ref
    print *,'Mref_V%Qmu_ref(',i,') = ',Qmu_ref
    print *,'Mref_V%vph_ref(',i,') = ',vph_ref
    print *,'Mref_V%vsh_ref(',i,') = ',vsh_ref
    print *,'Mref_V%eta_ref(',i,') = ',eta_ref
    print *

  enddo
  close(10)

  end program extract_model_reference

