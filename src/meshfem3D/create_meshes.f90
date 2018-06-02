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

  subroutine create_meshes()

  use meshfem3D_par

  implicit none

  ! local parameters
  integer :: ipass
  integer :: ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Radial Meshing parameters:'
    write(IMAIN,*) '  CHUNK WIDTH XI/ETA =', ANGULAR_WIDTH_XI_IN_DEGREES,'/',ANGULAR_WIDTH_ETA_IN_DEGREES
    write(IMAIN,*) '  NEX XI/ETA =', NEX_XI,'/',NEX_ETA
    write(IMAIN,*)
    write(IMAIN,*) '  NER_CRUST:               ', NER_CRUST
    write(IMAIN,*) '  NER_80_MOHO:             ', NER_80_MOHO
    write(IMAIN,*) '  NER_220_80:              ', NER_220_80
    write(IMAIN,*) '  NER_400_220:             ', NER_400_220
    write(IMAIN,*) '  NER_600_400:             ', NER_600_400
    write(IMAIN,*) '  NER_670_600:             ', NER_670_600
    write(IMAIN,*) '  NER_771_670:             ', NER_771_670
    write(IMAIN,*) '  NER_TOPDDOUBLEPRIME_771: ', NER_TOPDDOUBLEPRIME_771
    write(IMAIN,*) '  NER_CMB_TOPDDOUBLEPRIME: ', NER_CMB_TOPDDOUBLEPRIME
    write(IMAIN,*) '  NER_OUTER_CORE:          ', NER_OUTER_CORE
    write(IMAIN,*) '  NER_TOP_CENTRAL_CUBE_ICB:', NER_TOP_CENTRAL_CUBE_ICB
    write(IMAIN,*) '  SUPPRESS_CRUSTAL_MESH:   ', SUPPRESS_CRUSTAL_MESH
    write(IMAIN,*)
    write(IMAIN,*) '  R_CENTRAL_CUBE = ', sngl(R_CENTRAL_CUBE/1000.d0),' km'
    write(IMAIN,*)
    write(IMAIN,*) 'Mesh resolution:'
    write(IMAIN,*) '  DT = ',DT
    write(IMAIN,*) '  Minimum period = ', &
                   max(ANGULAR_WIDTH_ETA_IN_DEGREES,ANGULAR_WIDTH_XI_IN_DEGREES)/90.0 * 256.0/min(NEX_ETA,NEX_XI) * 17.0,' (s)'
    write(IMAIN,*)
    write(IMAIN,*) '  MIN_ATTENUATION_PERIOD = ',MIN_ATTENUATION_PERIOD
    write(IMAIN,*) '  MAX_ATTENUATION_PERIOD = ',MAX_ATTENUATION_PERIOD
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get addressing for this process
  ichunk = ichunk_slice(myrank)
  iproc_xi = iproc_xi_slice(myrank)
  iproc_eta = iproc_eta_slice(myrank)

  ! volume of the final mesh, and Earth mass computed in the final mesh
  ! and gravity integrals
  volume_total = ZERO
  Earth_mass_total = ZERO
  Earth_center_of_mass_x_total = ZERO
  Earth_center_of_mass_y_total = ZERO
  Earth_center_of_mass_z_total = ZERO

  ! make sure everybody is synchronized
  call synchronize_all()

  !----
  !----  loop on all the regions of the mesh
  !----

  ! number of regions in full Earth
  do iregion_code = 1,MAX_NUM_REGIONS

    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '*******************************************'
      write(IMAIN,*) 'creating mesh in region ',iregion_code
      select case (iregion_code)
        case (IREGION_CRUST_MANTLE)
          write(IMAIN,*) 'this region is the crust and mantle'
        case (IREGION_OUTER_CORE)
          write(IMAIN,*) 'this region is the outer core'
        case (IREGION_INNER_CORE)
          write(IMAIN,*) 'this region is the inner core'
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
      write(IMAIN,*) '*******************************************'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! number of spectral elements
    nspec = NSPEC_REGIONS(iregion_code)

    ! number of global GLL points
    nglob = NGLOB_REGIONS(iregion_code)

    ! compute maximum number of points
    npointot = nspec * NGLLX * NGLLY * NGLLZ

    ! use dynamic allocation to allocate memory for arrays
    allocate(idoubling(nspec), &
             ibool(NGLLX,NGLLY,NGLLZ,nspec), &
             xstore(NGLLX,NGLLY,NGLLZ,nspec), &
             ystore(NGLLX,NGLLY,NGLLZ,nspec), &
             zstore(NGLLX,NGLLY,NGLLZ,nspec), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating memory for arrays')

    ! this for non blocking MPI
    allocate(is_on_a_slice_edge(nspec),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating is_on_a_slice_edge array')


    ! create all the regions of the mesh
    ! perform two passes in this part to be able to save memory
    do ipass = 1,2
      call create_regions_mesh(iregion_code,npointot, &
                               NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                               NSPEC2DMAX_XMIN_XMAX(iregion_code),NSPEC2DMAX_YMIN_YMAX(iregion_code), &
                               NSPEC2D_BOTTOM(iregion_code),NSPEC2D_TOP(iregion_code), &
                               mod(iproc_xi_slice(myrank),2),mod(iproc_eta_slice(myrank),2), &
                               ipass)

      ! If we're in the request stage of CEM, exit.
      if (CEM_REQUEST) exit

    enddo

    ! deallocate arrays used for that region
    deallocate(idoubling)
    deallocate(ibool)
    deallocate(xstore)
    deallocate(ystore)
    deallocate(zstore)

    ! this for non blocking MPI
    deallocate(is_on_a_slice_edge)

    ! make sure everybody is synchronized
    call synchronize_all()

  ! end of loop on all the regions
  enddo

  end subroutine create_meshes
