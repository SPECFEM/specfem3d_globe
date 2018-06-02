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

  subroutine finalize_mesher()

  use meshfem3D_par
  use meshfem3D_models_par

  use manager_adios

  implicit none

  ! local parameters
  ! timing
  double precision :: tCPU
  double precision, external :: wtime

  integer :: NT_DUMP_ATTENUATION_optimal
  logical, parameter :: PRINT_INFO_TO_SCREEN = .false.


  !--- print number of points and elements in the mesh for each region
  if (myrank == 0) then

    ! center of mass distance
    distance_to_center_in_km = (sqrt(Earth_center_of_mass_x_total**2 + &
                                     Earth_center_of_mass_y_total**2 + &
                                     Earth_center_of_mass_z_total**2) / Earth_mass_total) / 1000.d0

    numelem_crust_mantle = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    numelem_outer_core = NSPEC_REGIONS(IREGION_OUTER_CORE)
    numelem_inner_core = NSPEC_REGIONS(IREGION_INNER_CORE)

    numelem_total = numelem_crust_mantle + numelem_outer_core + numelem_inner_core

    ! check volume of chunk
    write(IMAIN,*)
    write(IMAIN,*) 'calculated volume: ',volume_total
    if (NCHUNKS == 6 .and. .not. ELLIPTICITY .and. .not. TOPOGRAPHY) &
        write(IMAIN,*) '     exact volume: ',(4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3)

    ! check total Earth mass
    if (NCHUNKS == 6) then
      write(IMAIN,*)
      write(IMAIN,*) 'computed total Earth mass for this density model and mesh: '
      write(IMAIN,*) '   ',sngl(Earth_mass_total),' kg'
      write(IMAIN,*) '   (should be not too far from 5.974E+24 kg)'
      write(IMAIN,*)
      ! take into account the fact that dimensions have been non-dimensionalized by dividing them by R_EARTH
      write(IMAIN,*) 'average density for this density model and mesh: '
      write(IMAIN,*) '   ',sngl(Earth_mass_total / (volume_total * R_EARTH**3)),' kg/m3'
      write(IMAIN,*) '   (should be not too far from 5514 kg/m3)'
      write(IMAIN,*)
      write(IMAIN,*) 'position of the center of mass of the Earth for this density model and mesh: '
      write(IMAIN,*) '   x = ',sngl((Earth_center_of_mass_x_total / Earth_mass_total) / 1000.d0),' km'
      write(IMAIN,*) '   y = ',sngl((Earth_center_of_mass_y_total / Earth_mass_total) / 1000.d0),' km'
      write(IMAIN,*) '   z = ',sngl((Earth_center_of_mass_z_total / Earth_mass_total) / 1000.d0),' km'
      write(IMAIN,*) '   distance to center = ',sngl(distance_to_center_in_km),' km'
    endif

    ! info output
    write(IMAIN,*)
    write(IMAIN,*) 'Repartition of elements in regions:'
    write(IMAIN,*) '----------------------------------'
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in each slice: ',numelem_total
    write(IMAIN,*)
    write(IMAIN,*) ' - crust and mantle: ',sngl(100.d0*dble(numelem_crust_mantle)/dble(numelem_total)),' %'
    write(IMAIN,*) ' - outer core: ',sngl(100.d0*dble(numelem_outer_core)/dble(numelem_total)),' %'
    write(IMAIN,*) ' - inner core: ',sngl(100.d0*dble(numelem_inner_core)/dble(numelem_total)),' %'
    write(IMAIN,*)
    write(IMAIN,*) 'for some mesh statistics, see comments in file OUTPUT_FILES/values_from_mesher.h'
    write(IMAIN,*)

    ! load balancing
    write(IMAIN,*) 'Load balancing = 100 % by definition'
    write(IMAIN,*)

    write(IMAIN,*)
    write(IMAIN,*) 'the time step of the solver will be DT = ',sngl(DT)
    write(IMAIN,*)

    ! write information about precision used for floating-point operations
    if (CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)
    call flush_IMAIN()

    ! create include file for the solver
    call save_header_file(NSPEC_REGIONS,NGLOB_REGIONS,NPROC,NPROCTOT, &
                          static_memory_size, &
                          NSPEC2D_TOP,NSPEC2D_BOTTOM, &
                          NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_XMIN_XMAX, &
                          NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                          NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
                          NSPEC_INNER_CORE_ATTENUATION, &
                          NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
                          NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
                          NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
                          NSPEC_CRUST_MANTLE_ADJOINT, &
                          NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                          NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
                          NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
                          NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
                          NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION,NT_DUMP_ATTENUATION_optimal, &
                          PRINT_INFO_TO_SCREEN)

  endif   ! end of section executed by main process only

  ! finishes with gravity integrals
  if (GRAVITY_INTEGRALS) call finalize_gravity_integrals()

  ! deallocate arrays used for mesh generation
  deallocate(addressing)
  deallocate(ichunk_slice)
  deallocate(iproc_xi_slice)
  deallocate(iproc_eta_slice)

  ! elapsed time since beginning of mesh generation
  if (myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for mesh generation and buffer creation in seconds = ',tCPU
    write(IMAIN,"(' Elapsed time for mesh generation and buffer creation in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
              int(tCPU/3600),int( (tCPU - int(tCPU/3600)*3600)/60 ),int(tCPU - int(tCPU/60) * 60)
    write(IMAIN,*)
    write(IMAIN,*) 'End of mesh generation'
    write(IMAIN,*)
    call flush_IMAIN()

    ! close main output file
    close(IMAIN)
  endif

  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

  if (ADIOS_ENABLED) then
    call finalize_adios()
  endif

  end subroutine finalize_mesher

