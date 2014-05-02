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

  subroutine finalize_mesher()

  use meshfem3D_par
  use meshfem3D_models_par

  implicit none

  ! local parameters
  ! timing
  double precision :: tCPU
  double precision, external :: wtime

  ! for Roland_Sylvain integrals
  double precision :: g_x, g_y, g_z, G_xx, G_yy, G_zz, G_xy, G_xz, G_yz, real_altitude_of_observ_point

  !--- print number of points and elements in the mesh for each region
  if(myrank == 0) then

    ! check volume of chunk
    write(IMAIN,*)
    write(IMAIN,*) 'calculated volume: ',volume_total
    if(NCHUNKS == 6 .and. .not. ELLIPTICITY .and. .not. TOPOGRAPHY) &
        write(IMAIN,*) '     exact volume: ',(4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3)

    ! check total Earth mass
    if(NCHUNKS == 6) then
      write(IMAIN,*)
      write(IMAIN,*) 'computed total Earth mass for this density model and mesh: ',Earth_mass_total,' kg'
      write(IMAIN,*) '   (should be not too far from 5.97E+24 kg)'
      write(IMAIN,*)
      ! take into account the fact that dimensions have been non-dimensionalized by dividing them by R_EARTH
      write(IMAIN,*) 'average density for this density model and mesh: ',Earth_mass_total / (volume_total * R_EARTH**3),' kg/m3'
      write(IMAIN,*) '   (should be not too far from 5514 kg/m3)'
    endif

!! DK DK for Roland_Sylvain
    ! Roland_Sylvain integrals
    if(ROLAND_SYLVAIN) then

! in m.s-2
      g_x  = Roland_Sylvain_integr_total(1)
      g_y  = Roland_Sylvain_integr_total(2)
      g_z  = Roland_Sylvain_integr_total(3)

! in Eotvos = 1.e+9 s-2
      G_xx = Roland_Sylvain_integr_total(4) * SI_UNITS_TO_EOTVOS
      G_yy = Roland_Sylvain_integr_total(5) * SI_UNITS_TO_EOTVOS
      G_zz = Roland_Sylvain_integr_total(6) * SI_UNITS_TO_EOTVOS
      G_xy = Roland_Sylvain_integr_total(7) * SI_UNITS_TO_EOTVOS
      G_xz = Roland_Sylvain_integr_total(8) * SI_UNITS_TO_EOTVOS
      G_yz = Roland_Sylvain_integr_total(9) * SI_UNITS_TO_EOTVOS

      write(IMAIN,*)
      write(IMAIN,*) 'computed total Roland_Sylvain integral g_x  = ',g_x,' m.s-2'
      write(IMAIN,*) 'computed total Roland_Sylvain integral g_y  = ',g_y,' m.s-2'
      write(IMAIN,*) 'computed total Roland_Sylvain integral g_z  = ',g_z,' m.s-2'
      write(IMAIN,*)
      write(IMAIN,*) 'computed norm of g vector = ',sqrt(g_x**2 + g_y**2 + g_z**2),' m.s-2'

      real_altitude_of_observ_point = sqrt(x_observation**2 + y_observation**2 + z_observation**2)
! gravity force vector norm decays approximately as (r / r_prime)^2 above the surface of the Earth
      write(IMAIN,*) '  (should be not too far from ', &
                             sngl(STANDARD_GRAVITY_EARTH * (R_UNIT_SPHERE / real_altitude_of_observ_point)**2),' m.s-2)'

      write(IMAIN,*)
      write(IMAIN,*) 'computed total Roland_Sylvain integral G_xx = ',G_xx,' Eotvos'
      write(IMAIN,*) 'computed total Roland_Sylvain integral G_yy = ',G_yy,' Eotvos'
      write(IMAIN,*) 'computed total Roland_Sylvain integral G_zz = ',G_zz,' Eotvos'
      write(IMAIN,*)
      write(IMAIN,*) 'G tensor should be traceless, G_xx + G_yy + G_zz = 0.'
      write(IMAIN,*) 'Actual sum obtained = ',G_xx + G_yy + G_zz
      if(max(abs(G_xx),abs(G_yy),abs(G_zz)) > TINYVAL) write(IMAIN,*) ' i.e., ', &
             sngl(100.d0*abs(G_xx + G_yy + G_zz) / max(abs(G_xx),abs(G_yy),abs(G_zz))),'% of max(abs(Gxx),abs(Gyy),abs(Gzz))'
      write(IMAIN,*)
      write(IMAIN,*) 'computed total Roland_Sylvain integral G_xy = ',G_xy,' Eotvos'
      write(IMAIN,*) 'computed total Roland_Sylvain integral G_xz = ',G_xz,' Eotvos'
      write(IMAIN,*) 'computed total Roland_Sylvain integral G_yz = ',G_yz,' Eotvos'

    endif

    ! info output
    numelem_crust_mantle = NSPEC(IREGION_CRUST_MANTLE)
    numelem_outer_core = NSPEC(IREGION_OUTER_CORE)
    numelem_inner_core = NSPEC(IREGION_INNER_CORE)

    numelem_total = numelem_crust_mantle + numelem_outer_core + numelem_inner_core

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
    write(IMAIN,*) 'total number of time steps in the solver will be: ',NSTEP
    write(IMAIN,*)

    write(IMAIN,*)
    write(IMAIN,*) 'time-stepping of the solver will be: ',DT
    write(IMAIN,*)

    ! write information about precision used for floating-point operations
    if(CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)
    call flush_IMAIN()

    ! evaluate the amount of static memory needed by the solver
    call memory_eval(doubling_index,this_region_has_a_doubling, &
                     ner,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                     ratio_sampling_array,NPROCTOT, &
                     NSPEC,NGLOB, &
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
                     NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION, &
                     NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                     static_memory_size)

    ! create include file for the solver
    call save_header_file(NSPEC,NGLOB,NPROC,NPROCTOT, &
                          NSOURCES, &
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
                          NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION )

  endif   ! end of section executed by main process only

  ! deallocate arrays used for mesh generation
  deallocate(addressing)
  deallocate(ichunk_slice)
  deallocate(iproc_xi_slice)
  deallocate(iproc_eta_slice)

  ! elapsed time since beginning of mesh generation
  if(myrank == 0) then
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
    call adios_cleanup()
  endif

  end subroutine finalize_mesher

