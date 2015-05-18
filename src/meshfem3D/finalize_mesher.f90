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

  subroutine finalize_mesher()

  use meshfem3D_par
  use meshfem3D_models_par

  implicit none

  ! local parameters
  ! timing
  double precision :: tCPU
  double precision, external :: wtime

  ! for gravity integrals
  ! take into account the fact that the density and the radius of the Earth have previously been non-dimensionalized
  ! for the gravity vector force, a distance is involved in the dimensions
  double precision, parameter :: nondimensionalizing_factor_gi  = RHOAV * R_EARTH
  ! for the second-order gravity tensor, no distance is involved in the dimensions
  double precision, parameter :: nondimensionalizing_factor_Gij = RHOAV

  double precision, parameter :: scaling_factor_gi = GRAV * nondimensionalizing_factor_gi
  double precision, parameter :: scaling_factor_Gij_Eotvos = GRAV * nondimensionalizing_factor_Gij * SI_UNITS_TO_EOTVOS

  double precision :: real_altitude_of_observ_point,distance_to_center_in_km

  integer :: ixval,iyval,ichunkval

  integer :: NT_DUMP_ATTENUATION_optimal
  logical, parameter :: PRINT_INFO_TO_SCREEN = .false.

  if (GRAVITY_INTEGRALS) then

    ! multiply by the gravitational constant in S.I. units i.e. in m3 kg-1 s-2
    ! and also take into account the fact that the density and the radius of the Earth have previously been non-dimensionalized
    ! the final result is in m.s-2 i.e. in S.I. units
    g_x(:,:,:) = g_x(:,:,:) * scaling_factor_gi
    g_y(:,:,:) = g_y(:,:,:) * scaling_factor_gi
    g_z(:,:,:) = g_z(:,:,:) * scaling_factor_gi

    ! the final result is in Eotvos = 1.e+9 s-2
    G_xx(:,:,:) = G_xx(:,:,:) * scaling_factor_Gij_Eotvos
    G_yy(:,:,:) = G_yy(:,:,:) * scaling_factor_Gij_Eotvos
    G_zz(:,:,:) = G_zz(:,:,:) * scaling_factor_Gij_Eotvos
    G_xy(:,:,:) = G_xy(:,:,:) * scaling_factor_Gij_Eotvos
    G_xz(:,:,:) = G_xz(:,:,:) * scaling_factor_Gij_Eotvos
    G_yz(:,:,:) = G_yz(:,:,:) * scaling_factor_Gij_Eotvos

    ! use an MPI reduction to compute the total value of the integral into a temporary array
    ! and then copy it back into the original array
    call sum_all_3Darray_dp(g_x,temporary_array_for_sum,NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX)
    if (myrank == 0) g_x(:,:,:) = temporary_array_for_sum(:,:,:)

    call sum_all_3Darray_dp(g_y,temporary_array_for_sum,NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX)
    if (myrank == 0) g_y(:,:,:) = temporary_array_for_sum(:,:,:)

    call sum_all_3Darray_dp(g_z,temporary_array_for_sum,NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX)
    if (myrank == 0) g_z(:,:,:) = temporary_array_for_sum(:,:,:)

    call sum_all_3Darray_dp(G_xx,temporary_array_for_sum,NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX)
    if (myrank == 0) G_xx(:,:,:) = temporary_array_for_sum(:,:,:)

    call sum_all_3Darray_dp(G_yy,temporary_array_for_sum,NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX)
    if (myrank == 0) G_yy(:,:,:) = temporary_array_for_sum(:,:,:)

    call sum_all_3Darray_dp(G_zz,temporary_array_for_sum,NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX)
    if (myrank == 0) G_zz(:,:,:) = temporary_array_for_sum(:,:,:)

    call sum_all_3Darray_dp(G_xy,temporary_array_for_sum,NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX)
    if (myrank == 0) G_xy(:,:,:) = temporary_array_for_sum(:,:,:)

    call sum_all_3Darray_dp(G_xz,temporary_array_for_sum,NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX)
    if (myrank == 0) G_xz(:,:,:) = temporary_array_for_sum(:,:,:)

    call sum_all_3Darray_dp(G_yz,temporary_array_for_sum,NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX)
    if (myrank == 0) G_yz(:,:,:) = temporary_array_for_sum(:,:,:)

  endif

  !--- print number of points and elements in the mesh for each region
  if (myrank == 0) then

    ! check volume of chunk
    write(IMAIN,*)
    write(IMAIN,*) 'calculated volume: ',volume_total
    if (NCHUNKS == 6 .and. .not. ELLIPTICITY .and. .not. TOPOGRAPHY) &
        write(IMAIN,*) '     exact volume: ',(4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3)

    ! check total Earth mass
    if (NCHUNKS == 6) then
      write(IMAIN,*)
      write(IMAIN,*) 'computed total Earth mass for this density model and mesh: ',Earth_mass_total,' kg'
      write(IMAIN,*) '   (should be not too far from 5.974E+24 kg)'
      write(IMAIN,*)
      ! take into account the fact that dimensions have been non-dimensionalized by dividing them by R_EARTH
      write(IMAIN,*) 'average density for this density model and mesh: ',Earth_mass_total / (volume_total * R_EARTH**3),' kg/m3'
      write(IMAIN,*) '   (should be not too far from 5514 kg/m3)'
      write(IMAIN,*)
      write(IMAIN,*) 'position of the center of mass of the Earth for this density model and mesh: '
      write(IMAIN,*) '   x = ',(Earth_center_of_mass_x_total / Earth_mass_total) / 1000.d0,' km'
      write(IMAIN,*) '   y = ',(Earth_center_of_mass_y_total / Earth_mass_total) / 1000.d0,' km'
      write(IMAIN,*) '   z = ',(Earth_center_of_mass_z_total / Earth_mass_total) / 1000.d0,' km'
      distance_to_center_in_km = (sqrt(Earth_center_of_mass_x_total**2 + Earth_center_of_mass_y_total**2 + &
                                                      Earth_center_of_mass_z_total**2) / Earth_mass_total) / 1000.d0
      write(IMAIN,*) '   distance to center = ',distance_to_center_in_km,' km'
      if (GRAVITY_INTEGRALS .and. .not. ONLY_COMPUTE_CENTER_OF_MASS .and. distance_to_center_in_km > 0.01d0) &
        stop 'Error: center of mass of the model is not located in zero for gravity integrals, aborting...'
      write(IMAIN,*)
    endif

!! DK DK for gravity integrals
    if (GRAVITY_INTEGRALS) then

      temporary_array_for_sum(:,:,:) = sqrt(g_x(:,:,:)**2 + g_y(:,:,:)**2 + g_z(:,:,:)**2)
      write(IMAIN,*)
      write(IMAIN,*) 'minval of norm of g vector on whole observation surface = ',minval(temporary_array_for_sum),' m.s-2'
      write(IMAIN,*) 'maxval of norm of g vector on whole observation surface = ',maxval(temporary_array_for_sum),' m.s-2'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_xx on whole observation surface = ',minval(G_xx),' Eotvos'
      write(IMAIN,*) 'maxval of G_xx on whole observation surface = ',maxval(G_xx),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_yy on whole observation surface = ',minval(G_yy),' Eotvos'
      write(IMAIN,*) 'maxval of G_yy on whole observation surface = ',maxval(G_yy),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_zz on whole observation surface = ',minval(G_zz),' Eotvos'
      write(IMAIN,*) 'maxval of G_zz on whole observation surface = ',maxval(G_zz),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_xy on whole observation surface = ',minval(G_xy),' Eotvos'
      write(IMAIN,*) 'maxval of G_xy on whole observation surface = ',maxval(G_xy),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_xz on whole observation surface = ',minval(G_xz),' Eotvos'
      write(IMAIN,*) 'maxval of G_xz on whole observation surface = ',maxval(G_xz),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_yz on whole observation surface = ',minval(G_yz),' Eotvos'
      write(IMAIN,*) 'maxval of G_yz on whole observation surface = ',maxval(G_yz),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'Minval and maxval of trace of G, which in principle should be zero:'
      write(IMAIN,*)
      temporary_array_for_sum(:,:,:) = abs(G_xx(:,:,:) + G_yy(:,:,:) + G_zz(:,:,:))
      write(IMAIN,*) 'minval of abs(G_xx + G_yy + G_zz) on whole observation surface = ',minval(temporary_array_for_sum),' Eotvos'
      write(IMAIN,*) 'maxval of abs(G_xx + G_yy + G_zz) on whole observation surface = ',maxval(temporary_array_for_sum),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) '-----------------------------'
      write(IMAIN,*)
      write(IMAIN,*) 'displaying the fields computed at:'
      write(IMAIN,*) '    ix_observation = ',ixr,' out of ',NX_OBSERVATION
      write(IMAIN,*) '    iy_observation = ',iyr,' out of ',NY_OBSERVATION
      write(IMAIN,*) '    of mesh chunk ',ichunkr
      write(IMAIN,*)
      write(IMAIN,*) 'computed g_x  = ',g_x(ixr,iyr,ichunkr),' m.s-2'
      write(IMAIN,*) 'computed g_y  = ',g_y(ixr,iyr,ichunkr),' m.s-2'
      write(IMAIN,*) 'computed g_z  = ',g_z(ixr,iyr,ichunkr),' m.s-2'
      write(IMAIN,*)
      write(IMAIN,*) 'computed norm of g vector = ',sqrt(g_x(ixr,iyr,ichunkr)**2 + g_y(ixr,iyr,ichunkr)**2 + &
                                                                 g_z(ixr,iyr,ichunkr)**2),' m.s-2'

      real_altitude_of_observ_point = sqrt(x_observation(ixr,iyr,ichunkr)**2 + y_observation(ixr,iyr,ichunkr)**2 + &
                                                                 z_observation(ixr,iyr,ichunkr)**2)
! gravity force vector norm decays approximately as (r / r_prime)^2 above the surface of the Earth
      write(IMAIN,*) '  (should be not too far from ', &
                             sngl(STANDARD_GRAVITY_EARTH * (R_UNIT_SPHERE / real_altitude_of_observ_point)**2),' m.s-2)'

      write(IMAIN,*)
      write(IMAIN,*) 'computed G_xx = ',G_xx(ixr,iyr,ichunkr),' Eotvos'
      write(IMAIN,*) 'computed G_yy = ',G_yy(ixr,iyr,ichunkr),' Eotvos'
      write(IMAIN,*) 'computed G_zz = ',G_zz(ixr,iyr,ichunkr),' Eotvos'
      write(IMAIN,*)
      write(IMAIN,*) 'G tensor should be traceless, G_xx + G_yy + G_zz = 0.'
      write(IMAIN,*) 'Actual sum obtained = ',G_xx(ixr,iyr,ichunkr) + G_yy(ixr,iyr,ichunkr) + G_zz(ixr,iyr,ichunkr)
      if (max(abs(G_xx(ixr,iyr,ichunkr)),abs(G_yy(ixr,iyr,ichunkr)),abs(G_zz(ixr,iyr,ichunkr))) > TINYVAL) &
           write(IMAIN,*) ' i.e., ',sngl(100.d0*abs(G_xx(ixr,iyr,ichunkr) + G_yy(ixr,iyr,ichunkr) + G_zz(ixr,iyr,ichunkr)) / &
                                     max(abs(G_xx(ixr,iyr,ichunkr)),abs(G_yy(ixr,iyr,ichunkr)),abs(G_zz(ixr,iyr,ichunkr)))), &
                                     '% of max(abs(G_xx),abs(G_yy),abs(G_zz))'
      write(IMAIN,*)
      write(IMAIN,*) 'computed G_xy = ',G_xy(ixr,iyr,ichunkr),' Eotvos'
      write(IMAIN,*) 'computed G_xz = ',G_xz(ixr,iyr,ichunkr),' Eotvos'
      write(IMAIN,*) 'computed G_yz = ',G_yz(ixr,iyr,ichunkr),' Eotvos'

      ! for future GMT display
      ! loop on all the chunks and then on all the observation nodes in each chunk
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_g_x_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval),g_x(ixval,iyval,ichunkval)
          enddo
        enddo
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_g_y_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval),g_y(ixval,iyval,ichunkval)
          enddo
        enddo
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_g_z_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval),g_z(ixval,iyval,ichunkval)
          enddo
        enddo
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_norm_of_g_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval), &
                             sqrt(g_x(ixval,iyval,ichunkval)**2 + g_y(ixval,iyval,ichunkval)**2 + g_z(ixval,iyval,ichunkval)**2)
          enddo
        enddo
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_xx_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval),G_xx(ixval,iyval,ichunkval)
          enddo
        enddo
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_yy_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval),G_yy(ixval,iyval,ichunkval)
          enddo
        enddo
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_zz_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval),G_zz(ixval,iyval,ichunkval)
          enddo
        enddo
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_xy_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval),G_xy(ixval,iyval,ichunkval)
          enddo
        enddo
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_xz_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval),G_xz(ixval,iyval,ichunkval)
          enddo
        enddo
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_yz_for_GMT.txt',status='unknown',action='write')
      do ichunkval = 1,NCHUNKS_MAX
        do iyval = 1,NY_OBSERVATION
          do ixval = 1,NX_OBSERVATION
            write(IOUT,*) lon_observation(ixval,iyval,ichunkval),lat_observation(ixval,iyval,ichunkval),G_yz(ixval,iyval,ichunkval)
          enddo
        enddo
      enddo
      close(unit=IOUT)

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
    call save_header_file(NSPEC,NGLOB,NPROC,NPROCTOT, &
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
    call adios_cleanup()
  endif

  end subroutine finalize_mesher

