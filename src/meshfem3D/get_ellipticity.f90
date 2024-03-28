!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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


  subroutine get_ellipticity(xelm,yelm,zelm,nspl,rspl,ellipicity_spline,ellipicity_spline2)

  use constants, only: NR_DENSITY,NGNOD,ONE,TWO

  implicit none

  integer,intent(in) :: nspl
  double precision,intent(inout) :: xelm(NGNOD)
  double precision,intent(inout) :: yelm(NGNOD)
  double precision,intent(inout) :: zelm(NGNOD)
  double precision,intent(in) :: rspl(NR_DENSITY),ellipicity_spline(NR_DENSITY),ellipicity_spline2(NR_DENSITY)

  ! local parameters
  integer :: ia
  double precision :: x,y,z

  do ia = 1,NGNOD

    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! adds ellipticity to position x/y/z
    call add_ellipticity(x,y,z,nspl,rspl,ellipicity_spline,ellipicity_spline2)

    xelm(ia) = x
    yelm(ia) = y
    zelm(ia) = z

  enddo

  end subroutine get_ellipticity

!
!-------------------------------------------------------------------------------------------------
!

  !> Hejun
  ! ellipticity at the GLL points
  ! JAN08, 2010
  subroutine get_ellipticity_gll(xstore,ystore,zstore,ispec,nspec,nspl,rspl,ellipicity_spline,ellipicity_spline2)

  use constants, only: NR_DENSITY,NGLLX,NGLLY,NGLLZ,ONE,TWO

  implicit none

  integer,intent(in) :: nspl
  integer,intent(in) :: ispec,nspec
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: xstore,ystore,zstore
  double precision,intent(in) :: rspl(NR_DENSITY),ellipicity_spline(NR_DENSITY),ellipicity_spline2(NR_DENSITY)

  ! local parameters
  integer :: i,j,k
  double precision :: x,y,z

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        x = xstore(i,j,k,ispec)
        y = ystore(i,j,k,ispec)
        z = zstore(i,j,k,ispec)

        ! adds ellipticity to position x/y/z
        call add_ellipticity(x,y,z,nspl,rspl,ellipicity_spline,ellipicity_spline2)

        xstore(i,j,k,ispec) = x
        ystore(i,j,k,ispec) = y
        zstore(i,j,k,ispec) = z

      enddo
    enddo
  enddo

  end subroutine get_ellipticity_gll

