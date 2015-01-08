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

  subroutine create_central_cube(myrank,ichunk,ispec,iaddx,iaddy,iaddz,ipass, &
                        nspec,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,R_CENTRAL_CUBE, &
                        iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,ratio_divide_central_cube, &
                        iMPIcut_xi,iMPIcut_eta,iboun, &
                        idoubling,iregion_code,xstore,ystore,zstore, &
                        shape3D,rmin,rmax,rhostore,dvpstore,&
                        kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                        xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                        gammaxstore,gammaystore,gammazstore,nspec_actually, &
                        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                        c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                        nspec_ani,nspec_stacey, &
                        rho_vp,rho_vs,xigll,yigll,zigll, &
                        ispec_is_tiso)

! creates the inner core cube of the mesh

  use meshfem3D_models_par

  implicit none

  integer :: ratio_divide_central_cube

! correct number of spectral elements in each block depending on chunk type
  integer nspec,nspec_stacey

  integer NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  integer NPROC_XI,NPROC_ETA

  double precision R_CENTRAL_CUBE

! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! topology of the elements
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

! code for the four regions of the mesh
  integer iregion_code

! Gauss-Lobatto-Legendre points and weights of integration
  double precision xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ)

! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)

  integer idoubling(nspec)

! for model density and anisotropy
  integer nspec_ani
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani) :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

! boundary locator
  logical iboun(6,nspec)

! arrays with mesh parameters
  integer nspec_actually
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_actually) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

! proc numbers for MPI
  integer myrank

! MPI cut-planes parameters along xi and along eta
  logical, dimension(2,nspec) :: iMPIcut_xi,iMPIcut_eta

! Stacey, indices for Clayton-Engquist absorbing conditions
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_stacey) :: rho_vp,rho_vs

  integer ispec
  integer iproc_xi,iproc_eta,ichunk,ipass

  logical, dimension(nspec) :: ispec_is_tiso

  !local parameters
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  ! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision :: rmin,rmax
  ! to define the central cube in the inner core
  double precision :: radius_cube
  double precision :: xgrid_central_cube,ygrid_central_cube,zgrid_central_cube
  integer ix,iy,iz,ia
  integer nx_central_cube,ny_central_cube,nz_central_cube
  ! the height at which the central cube is cut
  integer :: nz_inf_limit



  ! create the shape of a regular mesh element in the inner core
  call hex_nodes(iaddx,iaddy,iaddz)

  ! define vertical slice in central cube on current processor
  ! we can assume that NEX_XI = NEX_ETA, otherwise central cube cannot be defined
  nx_central_cube = NEX_PER_PROC_XI / ratio_divide_central_cube
  ny_central_cube = NEX_PER_PROC_ETA / ratio_divide_central_cube
  nz_central_cube = NEX_XI / ratio_divide_central_cube

  ! size of the cube along Cartesian axes before rotation
  radius_cube = (R_CENTRAL_CUBE / R_EARTH) / sqrt(3.d0)

  ! define spectral elements in central cube
  do iz = 0,2*nz_central_cube-2,2
    do iy = 0,2*ny_central_cube-2,2
      do ix = 0,2*nx_central_cube-2,2

        ! radii that define the shell, we know that we are in the central cube
        rmin = 0.d0
        rmax = R_CENTRAL_CUBE / R_EARTH

        ! loop over the NGNOD nodes
        do ia = 1,NGNOD

          ! flat cubed sphere with correct mapping
          call compute_coord_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia), &
                        xgrid_central_cube,ygrid_central_cube,zgrid_central_cube, &
                        iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,nx_central_cube,&
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
        ispec = ispec + 1
        if (ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in central cube creation')

        ! new get_flag_boundaries
        ! xmin & xmax
        if (ix == 0) then
          iMPIcut_xi(1,ispec) = .true.
          if (iproc_xi == 0) iboun(1,ispec)= .true.
        endif
        if (ix == 2*nx_central_cube-2) then
          iMPIcut_xi(2,ispec) = .true.
          if (iproc_xi == NPROC_XI-1) iboun(2,ispec)= .true.
        endif
        ! ymin & ymax
        if (iy == 0) then
          iMPIcut_eta(1,ispec) = .true.
          if (iproc_eta == 0) iboun(3,ispec)= .true.
        endif
        if (iy == 2*ny_central_cube-2) then
          iMPIcut_eta(2,ispec) = .true.
          if (iproc_eta == NPROC_ETA-1) iboun(4,ispec)= .true.
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
            idoubling(ispec) = IFLAG_BOTTOM_CENTRAL_CUBE
          else if (iz == 2*nz_central_cube-2) then
            idoubling(ispec) = IFLAG_TOP_CENTRAL_CUBE
          else if (iz > nz_inf_limit .and. iz < 2*nz_central_cube-2) then
            idoubling(ispec) = IFLAG_MIDDLE_CENTRAL_CUBE
          else
            idoubling(ispec) = IFLAG_IN_FICTITIOUS_CUBE
          endif
        else
          idoubling(ispec) = IFLAG_IN_FICTITIOUS_CUBE
        endif

        ! compute several rheological and geometrical properties for this spectral element
        call compute_element_properties(ispec,iregion_code,idoubling,ipass, &
                         xstore,ystore,zstore,nspec,myrank, &
                         xelm,yelm,zelm,shape3D,rmin,rmax,rhostore,dvpstore, &
                         kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                         xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,nspec_actually, &
                         c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                         c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                         c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                         nspec_ani,nspec_stacey, &
                         rho_vp,rho_vs, &
                         xigll,yigll,zigll,ispec_is_tiso)
      enddo
    enddo
  enddo

  end subroutine create_central_cube
