!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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


  subroutine create_meshes()

  use meshfem3D_par
  implicit none
  
  ! local parameters
  ! parameters needed to store the radii of the grid points
  ! in the spherically symmetric Earth
  integer, dimension(:), allocatable :: idoubling
  integer, dimension(:,:,:,:), allocatable :: ibool
  ! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore
  integer :: ier
    
  ! get addressing for this process
  ichunk = ichunk_slice(myrank)
  iproc_xi = iproc_xi_slice(myrank)
  iproc_eta = iproc_eta_slice(myrank)

  ! volume of the slice
  volume_total = ZERO

  ! make sure everybody is synchronized
  call sync_all()

  !----
  !----  loop on all the regions of the mesh
  !----

  ! number of regions in full Earth
  do iregion_code = 1,MAX_NUM_REGIONS

    if(myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '*******************************************'
      write(IMAIN,*) 'creating mesh in region ',iregion_code
      select case(iregion_code)
        case(IREGION_CRUST_MANTLE)
          write(IMAIN,*) 'this region is the crust and mantle'
        case(IREGION_OUTER_CORE)
          write(IMAIN,*) 'this region is the outer core'
        case(IREGION_INNER_CORE)
          write(IMAIN,*) 'this region is the inner core'
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
      write(IMAIN,*) '*******************************************'
      write(IMAIN,*)
    endif

    ! compute maximum number of points
    npointot = NSPEC(iregion_code) * NGLLX * NGLLY * NGLLZ

    ! use dynamic allocation to allocate memory for arrays
    allocate(idoubling(NSPEC(iregion_code)), &
            ibool(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)), &
            xstore(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)), &
            ystore(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)), &
            zstore(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)), &
            stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating memory for arrays')
    
    ! this for non blocking MPI
    allocate(is_on_a_slice_edge(NSPEC(iregion_code)), &
            stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating is_on_a_slice_edge array')
    

    ! create all the regions of the mesh
    ! perform two passes in this part to be able to save memory
    do ipass = 1,2
      call create_regions_mesh(iregion_code,ibool,idoubling,is_on_a_slice_edge, &
                          xstore,ystore,zstore,rmins,rmaxs, &
                          iproc_xi,iproc_eta,ichunk,NSPEC(iregion_code),nspec_tiso, &
                          volume_local,area_local_bottom,area_local_top, &
                          NGLOB(iregion_code),npointot, &
                          NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                          NSPEC2DMAX_XMIN_XMAX(iregion_code),NSPEC2DMAX_YMIN_YMAX(iregion_code), &
                          NSPEC2D_BOTTOM(iregion_code),NSPEC2D_TOP(iregion_code), &
                          NPROC_XI,NPROC_ETA, &
                          NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE,NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER, &
                          myrank,LOCAL_PATH,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
                          SAVE_MESH_FILES,NCHUNKS,INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS, &
                          R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,&
                          RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
                          ner,ratio_sampling_array,doubling_index,r_bottom,r_top,&
                          this_region_has_a_doubling,ratio_divide_central_cube, &
                          CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                          mod(iproc_xi_slice(myrank),2),mod(iproc_eta_slice(myrank),2), &
                          ipass)
    enddo

    ! checks number of anisotropic elements found in the mantle
    if(iregion_code /= IREGION_CRUST_MANTLE .and. nspec_tiso /= 0 ) &
      call exit_MPI(myrank,'found anisotropic elements outside of the mantle')

    if( TRANSVERSE_ISOTROPY ) then
      if(iregion_code == IREGION_CRUST_MANTLE .and. nspec_tiso == 0) &
        call exit_MPI(myrank,'found no anisotropic elements in the mantle')
    endif

    ! computes total area and volume
    call compute_area(myrank,NCHUNKS,iregion_code, &
                              area_local_bottom,area_local_top,&
                              volume_local,volume_total, &
                              RCMB,RICB,R_CENTRAL_CUBE)

    ! create chunk buffers if more than one chunk
    if(NCHUNKS > 1) then
      call create_chunk_buffers(iregion_code,NSPEC(iregion_code),ibool,idoubling, &
                              xstore,ystore,zstore, &
                              NGLOB(iregion_code), &
                              NSPEC2DMAX_XMIN_XMAX(iregion_code),NSPEC2DMAX_YMIN_YMAX(iregion_code), &
                              NPROC_XI,NPROC_ETA,NPROC,NPROCTOT, &
                              NGLOB1D_RADIAL_CORNER,maxval(NGLOB1D_RADIAL_CORNER(iregion_code,:)), &
                              NGLOB2DMAX_XMIN_XMAX(iregion_code),NGLOB2DMAX_YMIN_YMAX(iregion_code), &
                              myrank,LOCAL_PATH,addressing, &
                              ichunk_slice,iproc_xi_slice,iproc_eta_slice,NCHUNKS)
    else
      if(myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'only one chunk, no need to create chunk buffers'
        write(IMAIN,*)
      endif
    endif

    ! deallocate arrays used for that region
    deallocate(idoubling)
    deallocate(ibool)
    deallocate(xstore)
    deallocate(ystore)
    deallocate(zstore)

    ! this for non blocking MPI
    deallocate(is_on_a_slice_edge)

    ! make sure everybody is synchronized
    call sync_all()

  ! end of loop on all the regions
  enddo

  end subroutine create_meshes
