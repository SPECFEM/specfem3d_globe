!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

  subroutine compute_area(myrank,NCHUNKS,iregion_code, &
                                    area_local_bottom,area_local_top,&
                                    volume_local,volume_total, &
                                    RCMB,RICB,R_CENTRAL_CUBE)

  use meshfem3D_models_par

  implicit none

  integer :: myrank,NCHUNKS,iregion_code

  double precision :: area_local_bottom,area_local_top,volume_local
  double precision :: volume_total
  double precision :: RCMB,RICB,R_CENTRAL_CUBE

  ! local parameters
  double precision :: volume_total_region,area_total_bottom,area_total_top

  ! use an MPI reduction to compute the total area and volume
  volume_total_region = ZERO
  area_total_bottom   = ZERO
  area_total_top   = ZERO

  call sum_all_dp(area_local_bottom,area_total_bottom)
  call sum_all_dp(area_local_top,area_total_top)
  call sum_all_dp(volume_local,volume_total_region)

  if(myrank == 0) then
    !   sum volume over all the regions
    volume_total = volume_total + volume_total_region

    !   check volume of chunk, and bottom and top area
    write(IMAIN,*)
    write(IMAIN,*) '   calculated top area: ',area_total_top

    ! compare to exact theoretical value
    if(NCHUNKS == 6 .and. .not. TOPOGRAPHY) then
      select case(iregion_code)
        case(IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*R_UNIT_SPHERE**2
        case(IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2
        case(IREGION_INNER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
    endif

    write(IMAIN,*) 'calculated bottom area: ',area_total_bottom

    ! compare to exact theoretical value
    if(NCHUNKS == 6 .and. .not. TOPOGRAPHY) then
      select case(iregion_code)
        case(IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2
        case(IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2
        case(IREGION_INNER_CORE)
          write(IMAIN,*) '            similar area (central cube): ',dble(NCHUNKS)*(2.*(R_CENTRAL_CUBE / R_EARTH)/sqrt(3.))**2
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
    endif
    call flush_IMAIN()

  endif

  end subroutine compute_area

!=====================================================================

  ! computes total Earth mass

  subroutine compute_total_Earth_mass(myrank,Earth_mass_local,Earth_mass_total)

  use meshfem3D_models_par

  implicit none

  integer :: myrank

  double precision :: Earth_mass_local
  double precision :: Earth_mass_total

  ! local parameters
  double precision :: Earth_mass_total_region

  ! use an MPI reduction to compute the total Earth mass
  Earth_mass_total_region = ZERO
  call sum_all_dp(Earth_mass_local,Earth_mass_total_region)

  !   sum volume over all the regions
  if(myrank == 0) Earth_mass_total = Earth_mass_total + Earth_mass_total_region

  end subroutine compute_total_Earth_mass

