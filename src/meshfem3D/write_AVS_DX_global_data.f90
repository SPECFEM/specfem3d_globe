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

! create AVS or DX 3D data for the slice, to be recombined in postprocessing
  subroutine write_AVS_DX_global_data(prname,nspec,ibool,idoubling, &
                 xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot)

  use constants

  implicit none

  integer nspec
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! logical mask used to output global points only once
  integer npointot
  logical mask_ibool(npointot)

! numbering of global AVS or DX points
  integer num_ibool_AVS_DX(npointot)

  integer ispec
  integer iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer npoin,numpoin

! processor identification
  character(len=MAX_STRING_LEN) prname

! writing points
  open(unit=IOUT,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='unknown')

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

! mark global AVS or DX points
  do ispec = 1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)
    mask_ibool(iglob1) = .true.
    mask_ibool(iglob2) = .true.
    mask_ibool(iglob3) = .true.
    mask_ibool(iglob4) = .true.
    mask_ibool(iglob5) = .true.
    mask_ibool(iglob6) = .true.
    mask_ibool(iglob7) = .true.
    mask_ibool(iglob8) = .true.
  enddo

! count global number of AVS or DX points
  npoin = count(mask_ibool(:))

! number of points in AVS or DX file
  write(IOUT,*) npoin

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

! output global AVS or DX points
  numpoin = 0
  do ispec = 1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)
    if (.not. mask_ibool(iglob1)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob1) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(1,1,1,ispec)), &
              sngl(ystore(1,1,1,ispec)),sngl(zstore(1,1,1,ispec))
    endif
    if (.not. mask_ibool(iglob2)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob2) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(NGLLX,1,1,ispec)), &
              sngl(ystore(NGLLX,1,1,ispec)),sngl(zstore(NGLLX,1,1,ispec))
    endif
    if (.not. mask_ibool(iglob3)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob3) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(NGLLX,NGLLY,1,ispec)), &
              sngl(ystore(NGLLX,NGLLY,1,ispec)),sngl(zstore(NGLLX,NGLLY,1,ispec))
    endif
    if (.not. mask_ibool(iglob4)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob4) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(1,NGLLY,1,ispec)), &
              sngl(ystore(1,NGLLY,1,ispec)),sngl(zstore(1,NGLLY,1,ispec))
    endif
    if (.not. mask_ibool(iglob5)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob5) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(1,1,NGLLZ,ispec)), &
              sngl(ystore(1,1,NGLLZ,ispec)),sngl(zstore(1,1,NGLLZ,ispec))
    endif
    if (.not. mask_ibool(iglob6)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob6) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(NGLLX,1,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,1,NGLLZ,ispec)),sngl(zstore(NGLLX,1,NGLLZ,ispec))
    endif
    if (.not. mask_ibool(iglob7)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob7) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec)),sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
    endif
    if (.not. mask_ibool(iglob8)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob8) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(1,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(1,NGLLY,NGLLZ,ispec)),sngl(zstore(1,NGLLY,NGLLZ,ispec))
    endif
    mask_ibool(iglob1) = .true.
    mask_ibool(iglob2) = .true.
    mask_ibool(iglob3) = .true.
    mask_ibool(iglob4) = .true.
    mask_ibool(iglob5) = .true.
    mask_ibool(iglob6) = .true.
    mask_ibool(iglob7) = .true.
    mask_ibool(iglob8) = .true.
  enddo

! check that number of global points output is okay
  if (numpoin /= npoin) &
    call exit_MPI(myrank,'incorrect number of global points in AVS or DX file creation')

  close(IOUT)

! writing elements
  open(unit=IOUT,file=prname(1:len_trim(prname))//'AVS_DXelements.txt',status='unknown')

! number of elements in AVS or DX file
  write(IOUT,*) nspec

! output global AVS or DX elements
  do ispec = 1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)
    write(IOUT,*) ispec,idoubling(ispec),num_ibool_AVS_DX(iglob1), &
                  num_ibool_AVS_DX(iglob2),num_ibool_AVS_DX(iglob3), &
                  num_ibool_AVS_DX(iglob4),num_ibool_AVS_DX(iglob5), &
                  num_ibool_AVS_DX(iglob6),num_ibool_AVS_DX(iglob7), &
                  num_ibool_AVS_DX(iglob8)
  enddo

  close(IOUT)

  end subroutine write_AVS_DX_global_data

!
!-------------------------------------------------------------------------------------------------
!

!> Hejun
! write material information for GLL points
  subroutine write_AVS_DX_global_data_gll(prname,nspec, &
                 xstore,ystore,zstore,rhostore,kappavstore,muvstore,Qmustore, &
                 ATTENUATION)

  use constants

  implicit none

  integer nspec
  character(len=MAX_STRING_LEN) prname

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  real(kind=CUSTOM_REAL) kappavstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) muvstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) rhostore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision::  Qmustore(NGLLX,NGLLY,NGLLZ,nspec)

  logical :: ATTENUATION

  ! local parameters
  double precision,dimension(8):: vp,vs,rho,Qmu
  double precision:: vp_average,vs_average,rho_average,Qmu_average

  integer flag(NGLLX,NGLLY,NGLLZ,nspec)

  integer ispec,i,j,k
  integer iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer numpoin,nelem


! writing points
  open(unit=IOUT,file=prname(1:len_trim(prname))//'AVS_DXpoints_gll.txt',status='unknown')

! number of points in AVS or DX file
  write(IOUT,*) nspec*NGLLX*NGLLY*NGLLZ


! output global AVS or DX points
  numpoin = 0
  do ispec = 1,nspec
        do k = 1,NGLLZ
        do j = 1,NGLLY
        do i = 1,NGLLX
                numpoin = numpoin + 1
                write(IOUT,*) numpoin,sngl(xstore(i,j,k,ispec)), &
                              sngl(ystore(i,j,k,ispec)),sngl(zstore(i,j,k,ispec))
                flag(i,j,k,ispec) = numpoin
        enddo
        enddo
        enddo
  enddo

  close(IOUT)

! writing elements
  open(unit=IOUT,file=prname(1:len_trim(prname))//'AVS_DXelements_gll.txt',status='unknown')


! number of elements in AVS or DX file
  write(IOUT,*) nspec*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)

  nelem = 0
! output global AVS or DX elements
  do ispec = 1,nspec
    do k = 1,NGLLZ-1
      do j = 1,NGLLY-1
        do i = 1,NGLLX-1
          nelem = nelem + 1
          iglob1=flag(i,j,k,ispec)
          iglob2=flag(i+1,j,k,ispec)
          iglob3=flag(i+1,j+1,k,ispec)
          iglob4=flag(i,j+1,k,ispec)
          iglob5=flag(i,j,k+1,ispec)
          iglob6=flag(i+1,j,k+1,ispec)
          iglob7=flag(i+1,j+1,k+1,ispec)
          iglob8=flag(i,j+1,k+1,ispec)

          write(IOUT,*) nelem,iglob1, &
                        iglob2,iglob3,iglob4, &
                        iglob5,iglob6,iglob7,iglob8
        enddo
      enddo
    enddo
  enddo

  close(IOUT)

! writing elements property
  open(unit=IOUT,file=prname(1:len_trim(prname))//'AVS_DXmaterials_gll.txt',status='unknown')

! number of elements in AVS or DX file
  write(IOUT,*) nspec*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)

  nelem = 0
! output global AVS or DX elements
  do ispec = 1,nspec
        do k = 1,NGLLZ-1
        do j = 1,NGLLY-1
        do i = 1,NGLLX-1
               nelem = nelem + 1
                rho(1)=dble(rhostore(i,j,k,ispec))
                vs(1)=dble(sqrt(muvstore(i,j,k,ispec)/rhostore(i,j,k,ispec)))
                vp(1)=dble(sqrt(kappavstore(i,j,k,ispec)/rhostore(i,j,k,ispec)+4.d0*vs(1)*vs(1)/3.d0))

                rho(2)=dble(rhostore(i+1,j,k,ispec))
                vs(2)=dble(sqrt(muvstore(i+1,j,k,ispec)/rhostore(i+1,j,k,ispec)))
                vp(2)=dble(sqrt(kappavstore(i+1,j,k,ispec)/rhostore(i+1,j,k,ispec)+4.d0*vs(2)*vs(2)/3.d0))

                rho(3)=dble(rhostore(i+1,j+1,k,ispec))
                vs(3)=dble(sqrt(muvstore(i+1,j+1,k,ispec)/rhostore(i+1,j+1,k,ispec)))
                vp(3)=dble(sqrt(kappavstore(i+1,j+1,k,ispec)/rhostore(i+1,j+1,k,ispec)+4.d0*vs(3)*vs(3)/3.d0))

                rho(4)=dble(rhostore(i,j+1,k,ispec))
                vs(4)=dble(sqrt(muvstore(i,j+1,k,ispec)/rhostore(i,j+1,k,ispec)))
                vp(4)=dble(sqrt(kappavstore(i,j+1,k,ispec)/rhostore(i,j+1,k,ispec)+4.d0*vs(4)*vs(4)/3.d0))

                rho(5)=dble(rhostore(i,j,k+1,ispec))
                vs(5)=dble(sqrt(muvstore(i,j,k+1,ispec)/rhostore(i,j,k+1,ispec)))
                vp(5)=dble(sqrt(kappavstore(i,j,k+1,ispec)/rhostore(i,j,k+1,ispec)+4.d0*vs(5)*vs(5)/3.d0))

                rho(6)=dble(rhostore(i+1,j,k+1,ispec))
                vs(6)=dble(sqrt(muvstore(i+1,j,k+1,ispec)/rhostore(i+1,j,k+1,ispec)))
                vp(6)=dble(sqrt(kappavstore(i+1,j,k+1,ispec)/rhostore(i+1,j,k+1,ispec)+4.d0*vs(6)*vs(6)/3.d0))

                rho(7)=dble(rhostore(i+1,j+1,k+1,ispec))
                vs(7)=dble(sqrt(muvstore(i+1,j+1,k+1,ispec)/rhostore(i+1,j+1,k+1,ispec)))
                vp(7)=dble(sqrt(kappavstore(i+1,j+1,k+1,ispec)/rhostore(i+1,j+1,k+1,ispec)+4.d0*vs(7)*vs(7)/3.d0))

                rho(8)=dble(rhostore(i,j+1,k+1,ispec))
                vs(8)=dble(sqrt(muvstore(i,j+1,k+1,ispec)/rhostore(i,j+1,k+1,ispec)))
                vp(8)=dble(sqrt(kappavstore(i,j+1,k+1,ispec)/rhostore(i,j+1,k+1,ispec)+4.d0*vs(8)*vs(8)/3.d0))

                if (ATTENUATION) then
                        Qmu(1)=dble(Qmustore(i,j,k,ispec))
                        Qmu(2)=dble(Qmustore(i+1,j,k,ispec))
                        Qmu(3)=dble(Qmustore(i+1,j+1,k,ispec))
                        Qmu(4)=dble(Qmustore(i,j+1,k,ispec))
                        Qmu(5)=dble(Qmustore(i,j,k+1,ispec))
                        Qmu(6)=dble(Qmustore(i+1,j,k+1,ispec))
                        Qmu(7)=dble(Qmustore(i+1,j+1,k+1,ispec))
                        Qmu(8)=dble(Qmustore(i,j+1,k+1,ispec))
                        Qmu_average=Qmu(1)
                endif
                !rho_average=sum(rho(1:4))/4.d0
                !vp_average=sum(vp(1:4))/4.d0
                !vs_average=sum(vs(1:4))/4.d0
                rho_average=rho(1)
                vp_average=vp(1)
                vs_average=vs(1)

                if (ATTENUATION) then
                        write(IOUT,*) nelem,rho_average,vp_average,vs_average,Qmu_average
                else
                        write(IOUT,*) nelem,rho_average,vp_average,vs_average
                endif

        enddo
        enddo
        enddo
  enddo

  close(IOUT)

  end subroutine write_AVS_DX_global_data_gll


