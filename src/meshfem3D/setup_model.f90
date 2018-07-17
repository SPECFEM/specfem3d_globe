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

  subroutine setup_model()

  use meshfem3D_par

  implicit none

  ! user output
  if (myrank == 0) call sm_output_info()

  ! dynamic allocation of mesh arrays
  allocate(addressing(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1))
  allocate(ichunk_slice(0:NPROCTOT-1))
  allocate(iproc_xi_slice(0:NPROCTOT-1))
  allocate(iproc_eta_slice(0:NPROCTOT-1))

  ! creates global slice addressing for solver
  call create_addressing(NCHUNKS,NPROC,NPROC_ETA,NPROC_XI,NPROCTOT, &
                         addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice, &
                         OUTPUT_FILES)

  ! this for the different counters (which are now different if the superbrick is cut in the outer core)
  call setup_counters(NSPEC1D_RADIAL,NSPEC2D_XI,NSPEC2D_ETA,NGLOB1D_RADIAL, &
                      DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA, &
                      CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                      NPROCTOT,iproc_xi_slice,iproc_eta_slice, &
                      NSPEC1D_RADIAL_CORNER,NSPEC2D_XI_FACE, &
                      NSPEC2D_ETA_FACE,NGLOB1D_RADIAL_CORNER)

  ! distributes 3D models
  call meshfem3D_models_broadcast(NSPEC_REGIONS,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
                                  R80,R220,R670,RCMB,RICB, &
                                  LOCAL_PATH)

  ! infos about additional mesh optimizations
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'additional mesh optimizations'
    write(IMAIN,*)

    ! moho stretching
    write(IMAIN,*) 'moho:'
    if (CRUSTAL .and. CASE_3D) then
      if (REGIONAL_MOHO_MESH) then
        write(IMAIN,*) '  enforcing regional 3-layer crust'
        if (SUPPRESS_MOHO_STRETCHING) then
          write(IMAIN,*) '  no element stretching for 3-D moho surface'
        else
          if (HONOR_DEEP_MOHO) then
            write(IMAIN,*) '  incorporating element stretching for 3-D moho surface with deep moho'
          else
            write(IMAIN,*) '  incorporating element stretching for 3-D moho surface'
          endif
        endif
      else
        write(IMAIN,'(a,i1,a)') '   default ',NER_CRUST,'-layer crust'
        if (SUPPRESS_MOHO_STRETCHING .or. (.not. TOPOGRAPHY)) then
          write(IMAIN,*) '  no element stretching for 3-D moho surface'
        else
          write(IMAIN,*) '  incorporating element stretching for 3-D moho surface'
        endif
      endif
    else
      write(IMAIN,*) '  no element stretching for 3-D moho surface'
    endif
    write(IMAIN,*)

    ! internal topography
    write(IMAIN,*) 'internal topography 410/660:'
    if ((.not. SUPPRESS_INTERNAL_TOPOGRAPHY) .and. &
        (THREE_D_MODEL == THREE_D_MODEL_S362ANI .or. THREE_D_MODEL == THREE_D_MODEL_S362WMANI &
         .or. THREE_D_MODEL == THREE_D_MODEL_S362ANI_PREM .or. THREE_D_MODEL == THREE_D_MODEL_S29EA &
         .or. THREE_D_MODEL == THREE_D_MODEL_MANTLE_SH)) then
      write(IMAIN,*) '  incorporating element stretching for 3-D internal surfaces'
    else
      write(IMAIN,*) '  no element stretching for 3-D internal surfaces'
    endif
    call flush_IMAIN()
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*)
  endif
  call synchronize_all()

  end subroutine setup_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sm_output_info()

  use meshfem3D_models_par
  use meshfem3D_par, only: &
    MODEL,NEX_XI,NEX_ETA, &
    NPROC_XI,NPROC_ETA,NPROC,NCHUNKS,NPROCTOT, &
    R_CENTRAL_CUBE

  implicit none

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'There are ',NPROCTOT,' MPI processes'
  write(IMAIN,*) 'Processes are numbered from 0 to ',NPROCTOT-1
  write(IMAIN,*)
  write(IMAIN,*) 'There are ',NEX_XI,' elements along xi in each chunk'
  write(IMAIN,*) 'There are ',NEX_ETA,' elements along eta in each chunk'
  write(IMAIN,*)
  write(IMAIN,*) 'There are ',NPROC_XI,' slices along xi in each chunk'
  write(IMAIN,*) 'There are ',NPROC_ETA,' slices along eta in each chunk'
  write(IMAIN,*) 'There is a total of ',NPROC,' slices in each chunk'
  write(IMAIN,*) 'There are ',NCHUNKS,' chunks in the global mesh'
  write(IMAIN,*) 'There is a total of ',NPROCTOT,' slices in the global mesh'
  write(IMAIN,*)
  write(IMAIN,*) 'NGLLX = ',NGLLX
  write(IMAIN,*) 'NGLLY = ',NGLLY
  write(IMAIN,*) 'NGLLZ = ',NGLLZ
  write(IMAIN,*)
  write(IMAIN,*) 'Shape functions defined by NGNOD = ',NGNOD,' control nodes'
  write(IMAIN,*) 'Surface shape functions defined by NGNOD2D = ',NGNOD2D,' control nodes'
  write(IMAIN,*)

  ! model user parameters
  write(IMAIN,*) 'model: ',trim(MODEL)
  if (OCEANS) then
    write(IMAIN,*) '  incorporating the oceans using equivalent load'
  else
    write(IMAIN,*) '  no oceans'
  endif
  if (ELLIPTICITY) then
    write(IMAIN,*) '  incorporating ellipticity'
  else
    write(IMAIN,*) '  no ellipticity'
  endif
  if (TOPOGRAPHY) then
    write(IMAIN,*) '  incorporating surface topography'
  else
    write(IMAIN,*) '  no surface topography'
  endif
  if (GRAVITY) then
    write(IMAIN,*) '  incorporating self-gravitation (Cowling approximation)'
  else
    write(IMAIN,*) '  no self-gravitation'
  endif
  if (ROTATION) then
    write(IMAIN,*) '  incorporating rotation'
  else
    write(IMAIN,*) '  no rotation'
  endif
  if (ATTENUATION) then
    write(IMAIN,*) '  incorporating attenuation using ',N_SLS,' standard linear solids'
    if (ATTENUATION_3D) write(IMAIN,*)'  using 3D attenuation model'
  else
    write(IMAIN,*) '  no attenuation'
  endif
  write(IMAIN,*)

  ! model mesh parameters
  if (ISOTROPIC_3D_MANTLE) then
    write(IMAIN,*) '  incorporating 3-D lateral variations'
  else
    write(IMAIN,*) '  no 3-D lateral variations'
  endif
  if (HETEROGEN_3D_MANTLE) then
    write(IMAIN,*) '  incorporating heterogeneities in the mantle'
  else
    write(IMAIN,*) '  no heterogeneities in the mantle'
  endif
  if (CRUSTAL) then
    write(IMAIN,*) '  incorporating crustal variations'
  else
    write(IMAIN,*) '  no crustal variations'
  endif
  if (ONE_CRUST) then
    write(IMAIN,*) '  using one layer only in crust'
  else
    write(IMAIN,*) '  using unmodified 1D crustal model with two layers'
  endif
  if (TRANSVERSE_ISOTROPY) then
    write(IMAIN,*) '  incorporating anisotropy'
  else
    write(IMAIN,*) '  no anisotropy'
  endif
  if (ANISOTROPIC_INNER_CORE) then
    write(IMAIN,*) '  incorporating anisotropic inner core'
  else
    write(IMAIN,*) '  no inner-core anisotropy'
  endif
  if (ANISOTROPIC_3D_MANTLE) then
    write(IMAIN,*) '  incorporating anisotropic mantle'
  else
    write(IMAIN,*) '  no general mantle anisotropy'
  endif
  write(IMAIN,*)

  write(IMAIN,*) 'Reference radius of the Earth used is ',R_EARTH_KM,' km'
  write(IMAIN,*)
  write(IMAIN,*) 'Central cube is at a radius of ',R_CENTRAL_CUBE/1000.d0,' km'

  ! flushes I/O buffer
  call flush_IMAIN()

  end subroutine sm_output_info
