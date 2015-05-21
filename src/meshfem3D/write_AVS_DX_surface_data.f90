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

! create AVS or DX 2D data for the surface of the model
! to be recombined in postprocessing
  subroutine write_AVS_DX_surface_data(myrank,prname,nspec,iboun, &
     ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot,&
     rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
     ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
     RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
     RMIDDLE_CRUST,ROCEAN,iregion_code)

  use constants

  implicit none

  integer nspec,myrank
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

  logical iboun(6,nspec)
  logical ELLIPTICITY,ISOTROPIC_3D_MANTLE

  double precision RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771, &
       R400,R120,R80,RMOHO,RMIDDLE_CRUST,ROCEAN

  double precision r,rho,vp,vs,Qkappa,Qmu
  double precision vpv,vph,vsv,vsh,eta_aniso
  double precision x,y,z,theta,phi_dummy,cost,p20,ell,factor
  real(kind=CUSTOM_REAL) dvp,dvs

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  real(kind=CUSTOM_REAL) kappavstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) muvstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) rhostore(NGLLX,NGLLY,NGLLZ,nspec)

! logical mask used to output global points only once
  integer npointot
  logical mask_ibool(npointot)

! numbering of global AVS or DX points
  integer num_ibool_AVS_DX(npointot)

  integer ispec
  integer i,j,k,np
  integer, dimension(8) :: iglobval
  integer npoin,numpoin,nspecface,ispecface

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

! processor identification
  character(len=MAX_STRING_LEN) prname

  integer iregion_code

! writing points
  open(unit=IOUT,file=prname(1:len_trim(prname))//'AVS_DXpointssurface.txt',status='unknown')

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  nspecface = 0

! mark global AVS or DX points
  do ispec = 1,nspec
! only if at the surface (top plane)
  if (iboun(6,ispec)) then

    iglobval(5)=ibool(1,1,NGLLZ,ispec)
    iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
    iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

! element is at the surface
    nspecface = nspecface + 1
    mask_ibool(iglobval(5)) = .true.
    mask_ibool(iglobval(6)) = .true.
    mask_ibool(iglobval(7)) = .true.
    mask_ibool(iglobval(8)) = .true.

  endif
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
! only if at the surface
  if (iboun(6,ispec)) then

    iglobval(5)=ibool(1,1,NGLLZ,ispec)
    iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
    iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

! top face
  if (iboun(6,ispec)) then

    if (.not. mask_ibool(iglobval(5))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(5)) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(1,1,NGLLZ,ispec)), &
                    sngl(ystore(1,1,NGLLZ,ispec)),sngl(zstore(1,1,NGLLZ,ispec))
    endif

    if (.not. mask_ibool(iglobval(6))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(6)) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(NGLLX,1,NGLLZ,ispec)), &
                    sngl(ystore(NGLLX,1,NGLLZ,ispec)),sngl(zstore(NGLLX,1,NGLLZ,ispec))
    endif

    if (.not. mask_ibool(iglobval(7))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(7)) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec)), &
                    sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec)),sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
    endif

    if (.not. mask_ibool(iglobval(8))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(8)) = numpoin
      write(IOUT,*) numpoin,sngl(xstore(1,NGLLY,NGLLZ,ispec)), &
                    sngl(ystore(1,NGLLY,NGLLZ,ispec)),sngl(zstore(1,NGLLY,NGLLZ,ispec))
    endif

    mask_ibool(iglobval(5)) = .true.
    mask_ibool(iglobval(6)) = .true.
    mask_ibool(iglobval(7)) = .true.
    mask_ibool(iglobval(8)) = .true.

  endif

  endif
  enddo

! check that number of global points output is okay
  if (numpoin /= npoin) &
    call exit_MPI(myrank,'incorrect number of global points in AVS or DX file creation')

  close(IOUT)

! output global AVS or DX elements

! writing elements
  open(unit=IOUT,file=prname(1:len_trim(prname))//'AVS_DXelementssurface.txt',status='unknown')
  if (ISOTROPIC_3D_MANTLE) &
       open(unit=11,file=prname(1:len_trim(prname))//'AVS_DXelementssurface_dvp_dvs.txt',status='unknown')

! number of elements in AVS or DX file
  write(IOUT,*) nspecface

  ispecface = 0
  do ispec = 1,nspec
! only if at the surface
     if (iboun(6,ispec)) then

        iglobval(5)=ibool(1,1,NGLLZ,ispec)
        iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
        iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
        iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

                if (ISOTROPIC_3D_MANTLE) then
           !   pick a point within the element and get its radius
           r=dsqrt(xstore(2,2,2,ispec)**2+ystore(2,2,2,ispec)**2+zstore(2,2,2,ispec)**2)

           if (r > RCMB/R_EARTH .and. r < R_UNIT_SPHERE) then
              !     average over the element
              dvp = 0.0
              dvs = 0.0
              np  = 0
              do k=2,NGLLZ-1
                 do j=2,NGLLY-1
                    do i=2,NGLLX-1
                       np=np+1
                       x=xstore(i,j,k,ispec)
                       y=ystore(i,j,k,ispec)
                       z=zstore(i,j,k,ispec)
                       r=dsqrt(x*x+y*y+z*z)
                       ! take out ellipticity
                       if (ELLIPTICITY) then
                          call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi_dummy)
                          cost=dcos(theta)
! this is the Legendre polynomial of degree two, P2(cos(theta)), see the discussion above eq (14.4) in Dahlen and Tromp (1998)
                          p20=0.5d0*(3.0d0*cost*cost-1.0d0)
! get ellipticity using spline evaluation
                          call spline_evaluation(rspl,espl,espl2,nspl,r,ell)
! this is eq (14.4) in Dahlen and Tromp (1998)
                          factor=ONE-(TWO/3.0d0)*ell*p20
                          r=r/factor
                       endif


                       ! gets reference model values: rho,vpv,vph,vsv,vsh and eta_aniso
                       call meshfem3D_models_get1D_val(myrank,iregion_code,idoubling(ispec), &
                            r,rho,vpv,vph,vsv,vsh,eta_aniso, &
                            Qkappa,Qmu,RICB,RCMB, &
                            RTOPDDOUBLEPRIME,R80,R120,R220,R400,R600,R670,R771, &
                            RMOHO,RMIDDLE_CRUST,ROCEAN)

                       ! calculates isotropic values
                       vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                            + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
                       vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                            + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)

                       if (abs(rhostore(i,j,k,ispec))< 1.e-20) then
                          print *,' attention: rhostore close to zero',rhostore(i,j,k,ispec),r,i,j,k,ispec
                          dvp = 0.0
                          dvs = 0.0
                       else if (abs(sngl(vp))< 1.e-20) then
                          print *,' attention: vp close to zero',sngl(vp),r,i,j,k,ispec
                          dvp = 0.0
                       else if (abs(sngl(vs))< 1.e-20) then
                          print *,' attention: vs close to zero',sngl(vs),r,i,j,k,ispec
                          dvs = 0.0
                       else
                          dvp = dvp + (sqrt((kappavstore(i,j,k,ispec)+4.*muvstore(i,j,k,ispec)/3.) &
                               /rhostore(i,j,k,ispec)) - sngl(vp))/sngl(vp)
                          dvs = dvs + (sqrt(muvstore(i,j,k,ispec)/rhostore(i,j,k,ispec)) - sngl(vs))/sngl(vs)
                       endif

                    enddo
                 enddo
              enddo
              dvp = dvp / np
              dvs = dvs / np
           else
              dvp = 0.0
              dvs = 0.0
           endif
        endif

        ! top face
        ispecface = ispecface + 1
        write(IOUT,*) ispecface,idoubling(ispec),num_ibool_AVS_DX(iglobval(5)), &
                      num_ibool_AVS_DX(iglobval(6)),num_ibool_AVS_DX(iglobval(7)), &
                      num_ibool_AVS_DX(iglobval(8))
        if (ISOTROPIC_3D_MANTLE) write(11,*) ispecface,dvp,dvs

     endif
  enddo

! check that number of surface elements output is okay
  if (ispecface /= nspecface) &
    call exit_MPI(myrank,'incorrect number of surface elements in AVS or DX file creation')

  close(IOUT)
  if (ISOTROPIC_3D_MANTLE) close(11)

  end subroutine write_AVS_DX_surface_data

