!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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

  subroutine create_central_cube(ichunk,ispec_count,ipass, &
                                 nspec,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                                 iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA, &
                                 iMPIcut_xi,iMPIcut_eta, &
                                 idoubling,iregion_code, &
                                 rmin,rmax,R_CENTRAL_CUBE, &
                                 ispec_is_tiso)

! creates the inner core cube of the mesh

  use constants, only: NGNOD,myrank, &
    IFLAG_BOTTOM_CENTRAL_CUBE,IFLAG_TOP_CENTRAL_CUBE,IFLAG_MIDDLE_CENTRAL_CUBE,IFLAG_IN_FICTITIOUS_CUBE, &
    CHUNK_AB,CHUNK_AB_ANTIPODE, &
    CHUNK_AC,CHUNK_AC_ANTIPODE, &
    CHUNK_BC,CHUNK_BC_ANTIPODE

  use shared_parameters, only: ratio_divide_central_cube,R_PLANET

  use meshfem_par, only: &
    xstore,ystore,zstore

  use regions_mesh_par, only: &
    xigll,yigll,zigll,shape3D

  use regions_mesh_par2, only: &
    iboun

  implicit none

  integer, intent(in) :: ichunk,ipass

  integer, intent(inout) :: ispec_count

! correct number of spectral elements in each block depending on chunk type
  integer, intent(in) :: nspec
  integer, intent(in) :: NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer, intent(in) :: iproc_xi,iproc_eta
  integer, intent(in) :: NPROC_XI,NPROC_ETA

  double precision, intent(in) :: R_CENTRAL_CUBE

! MPI cut-planes parameters along xi and along eta
  logical, dimension(2,nspec), intent(inout) :: iMPIcut_xi,iMPIcut_eta

! code for the four regions of the mesh
  integer, intent(in) :: iregion_code

  integer, intent(inout) :: idoubling(nspec)

  ! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision,intent(inout) :: rmin,rmax

  logical, dimension(nspec), intent(inout) :: ispec_is_tiso

  ! local parameters
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  ! to define the central cube in the inner core
  double precision :: radius_cube
  double precision :: xgrid_central_cube,ygrid_central_cube,zgrid_central_cube
  integer :: ix,iy,iz,ia
  integer :: nx_central_cube,ny_central_cube,nz_central_cube
  ! the height at which the central cube is cut
  integer :: nz_inf_limit
  ! topology of the elements
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

  ! create the shape of a regular mesh element in the inner core
  call hex_nodes(iaddx,iaddy,iaddz)

  ! define vertical slice in central cube on current processor
  ! we can assume that NEX_XI = NEX_ETA, otherwise central cube cannot be defined
  nx_central_cube = NEX_PER_PROC_XI / ratio_divide_central_cube
  ny_central_cube = NEX_PER_PROC_ETA / ratio_divide_central_cube
  nz_central_cube = NEX_XI / ratio_divide_central_cube

  ! size of the cube along Cartesian axes before rotation
  radius_cube = (R_CENTRAL_CUBE / R_PLANET) / sqrt(3.d0)

  ! define spectral elements in central cube
  do iz = 0,2*nz_central_cube-2,2
    do iy = 0,2*ny_central_cube-2,2
      do ix = 0,2*nx_central_cube-2,2

        ! radii that define the shell, we know that we are in the central cube
        rmin = 0.d0
        rmax = R_CENTRAL_CUBE / R_PLANET

        ! loop over the NGNOD nodes
        do ia = 1,NGNOD

          ! flat cubed sphere with correct mapping
          call compute_coord_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia), &
                                          xgrid_central_cube,ygrid_central_cube,zgrid_central_cube, &
                                          iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,nx_central_cube, &
                                          ny_central_cube,nz_central_cube,radius_cube)

          if (ichunk == CHUNK_AB) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = + xgrid_central_cube
            zelm(ia) = + zgrid_central_cube

          else if (ichunk == CHUNK_AB_ANTIPODE) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = - xgrid_central_cube
            zelm(ia) = - zgrid_central_cube

          else if (ichunk == CHUNK_AC) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = - zgrid_central_cube
            zelm(ia) = + xgrid_central_cube

          else if (ichunk == CHUNK_AC_ANTIPODE) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = + zgrid_central_cube
            zelm(ia) = - xgrid_central_cube

          else if (ichunk == CHUNK_BC) then
            xelm(ia) = - zgrid_central_cube
            yelm(ia) = + ygrid_central_cube
            zelm(ia) = + xgrid_central_cube

          else if (ichunk == CHUNK_BC_ANTIPODE) then
            xelm(ia) = + zgrid_central_cube
            yelm(ia) = - ygrid_central_cube
            zelm(ia) = + xgrid_central_cube

          else
            call exit_MPI(myrank,'wrong chunk number in flat cubed sphere definition')
          endif

        enddo

        ! add one spectral element to the list
        ispec_count = ispec_count + 1
        if (ispec_count > nspec) call exit_MPI(myrank,'ispec greater than nspec in central cube creation')

        ! new get_flag_boundaries
        ! xmin & xmax
        if (ix == 0) then
          iMPIcut_xi(1,ispec_count) = .true.
          if (iproc_xi == 0) iboun(1,ispec_count)= .true.
        endif
        if (ix == 2*nx_central_cube-2) then
          iMPIcut_xi(2,ispec_count) = .true.
          if (iproc_xi == NPROC_XI-1) iboun(2,ispec_count)= .true.
        endif
        ! ymin & ymax
        if (iy == 0) then
          iMPIcut_eta(1,ispec_count) = .true.
          if (iproc_eta == 0) iboun(3,ispec_count)= .true.
        endif
        if (iy == 2*ny_central_cube-2) then
          iMPIcut_eta(2,ispec_count) = .true.
          if (iproc_eta == NPROC_ETA-1) iboun(4,ispec_count)= .true.
        endif

        ! define the doubling flag of this element
        ! only two active central cubes, the four others are fictitious

        ! determine where we cut the central cube to share it between CHUNK_AB & CHUNK_AB_ANTIPODE
        ! in the case of mod(NPROC_XI,2) /= 0, the cut is asymmetric and the bigger part is for CHUNK_AB
        nz_inf_limit = nz_central_cube
        if (mod(NPROC_XI,2) /= 0 .and. NPROC_XI > 1) then
          if (ichunk == CHUNK_AB) then
            nz_inf_limit = ((nz_central_cube*2)/NPROC_XI)*floor(NPROC_XI/2.d0)
          else if (ichunk == CHUNK_AB_ANTIPODE) then
            nz_inf_limit = ((nz_central_cube*2)/NPROC_XI)*ceiling(NPROC_XI/2.d0)
          endif
        endif

        if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
          if (iz == nz_inf_limit) then
            idoubling(ispec_count) = IFLAG_BOTTOM_CENTRAL_CUBE
          else if (iz == 2*nz_central_cube-2) then
            idoubling(ispec_count) = IFLAG_TOP_CENTRAL_CUBE
          else if (iz > nz_inf_limit .and. iz < 2*nz_central_cube-2) then
            idoubling(ispec_count) = IFLAG_MIDDLE_CENTRAL_CUBE
          else
            idoubling(ispec_count) = IFLAG_IN_FICTITIOUS_CUBE
          endif
        else
          idoubling(ispec_count) = IFLAG_IN_FICTITIOUS_CUBE
        endif

        ! compute several rheological and geometrical properties for this spectral element
        call compute_element_properties(ispec_count,iregion_code,idoubling,ipass, &
                                        xstore,ystore,zstore,nspec, &
                                        xelm,yelm,zelm,shape3D,rmin,rmax, &
                                        xigll,yigll,zigll,ispec_is_tiso)
      enddo
    enddo
  enddo

  end subroutine create_central_cube
