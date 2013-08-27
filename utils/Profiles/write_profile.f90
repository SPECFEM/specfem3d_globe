!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2013
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

!=====================================================================
! write_profile.f90
!
! This is a utility code written to ouput profiles through the model used by SPECFEM,
! as it were if it was sampled on a radial line through the Earth at 2x2 degree
! intervals, with discontinuities honored by the mesh and those of the MOHO and
! crust/ocean indicated by repeating points with different values.

! The code shortcuts through meshfem3D -> create_regions_mesh ->
! create_regular_elements -> compute_element_properties -> compute_element_properties ->
! get_model, cutting and pasting the relevant parts.

! This code IS NOT MAINTAINED and most likely will not work without modification as
! it calls functions that may have changed since it was last updated.
!
! I have tried to indicate where different parts come from to facilitate updates.
!   - vala hjorleifsdottir (vala@geofisica.unam.mx)
!=====================================================================

  program xwrite_profile

  use meshfem3D_models_par
  use shared_parameters

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"


! parameters needed to store the radii of the grid points
! in the spherically symmetric Earth
  integer idoubling

! proc numbers for MPI
  integer myrank,sizeprocs,ier


! for loop on all the slices
  integer iregion_code

! for ellipticity
!  integer nspl
!  double precision rspl(NR),espl(NR),espl2(NR)

! for stretching
 double precision gamma

  integer i,j
  double precision rho,drhodr,vp,vs,Qkappa,Qmu
  double precision vpv,vph,vsv,vsh,eta_aniso
  double precision dvp,dvs,drho
  real(kind=4) xcolat,xlon,xrad,dvpv,dvph,dvsv,dvsh
  double precision r,r_prem,r_moho,theta,phi,theta_deg,phi_deg,theta_degrees,phi_degrees
  double precision lat,lon,elevation
  double precision vpc,vsc,rhoc,moho
  integer NUMBER_OF_MESH_LAYERS

  double precision c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                   c34,c35,c36,c44,c45,c46,c55,c56,c66

  double precision rmin,rmax,rmax_last

! Attenuation values
  double precision, dimension(N_SLS)                     :: tau_s, tau_e
  double precision  T_c_source

  logical found_crust

  integer nit,ilayer,islice,iline,iline_icb,iline_cmb,iline_moho,iline_ocean
  integer ilayers_ocean,nlayers_ocean
  double precision delta,scaleval,r_ocean
  character(len=200) outfile

!--- from compute_element_properties
  ! Parameter used to decide whether this element is in the crust or not
  logical:: elem_in_crust,elem_in_mantle
!---

!--- from
  ! local parameters
  double precision xmesh,ymesh,zmesh
!---

! ************** PROGRAM STARTS HERE **************
  call MPI_INIT(ier)

!!! -- this part is from meshfem3D (read parfile and output info)

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_profile.txt',status='unknown')
  write(IMAIN,*) 'reading parameter file..'

! read the parameter file and compute additional parameters
  call read_compute_parameters()

  ! distributes 3D models
  print *,'Reading 3D models'
  call meshfem3D_models_broadcast(myrank,NSPEC, &
       MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,&
       R80,R220,R670,RCMB,RICB)
  print *,'Done reading 3D models'

!!    call meshfem3D_output_info(myrank,sizeprocs,NEX_XI,NEX_ETA, &
 !                               NPROC_XI,NPROC_ETA,NPROC,NCHUNKS,NPROCTOT,&
 !                               R_CENTRAL_CUBE)

!!! -- end part from meshfem3D

!!! -- this part is from create_regions_mesh

  if(ONE_CRUST) then
     NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS-1
  else
     NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
  endif
!!! -- end of part from create_regions_mesh

! loop over all theta (colatitude) and phis (longitude), every two degrees
! do i=0,89
  do i=15,15
  theta_degrees = 1.0d0 + i*2.0d0
   do j=45,45
!  do j=0,179
   phi_degrees = 1.0d0 + j*2.0d0

! open output file

  write(*,'(a,i04.4,a,i04.4)') &
     'OUTPUT_FILES/CARDS_th',int(theta_degrees),'_ph',int(phi_degrees)
  write(outfile,'(a,i04.4,a,i04.4)') &
     'OUTPUT_FILES/CARDS_th',int(theta_degrees),'_ph',int(phi_degrees)
  open(unit=57,file=outfile,status='unknown')

    rmax_last = 0.0d0
    theta = theta_degrees*TWO_PI/360.0d0
    phi   = phi_degrees  *TWO_PI/360.0d0

!  keep track of line number to be able to write out locations of discontinuities at end
    iline = 0

! read crustal models and topo models as they are needed to modify the depths of the discontinuities
    if(CRUSTAL) then
       call reduce(theta,phi)
       ! convert from rthetaphi to xyz, with r=-7km
       r = (1.0d0 - 7.0d0/R_EARTH_KM)
       call rthetaphi_2_xyz(xmesh,ymesh,zmesh,r,theta,phi)
       print *, 'xmesh,ymesh,zmesh,r,theta,phi',xmesh,ymesh,zmesh,r,theta,phi
       elem_in_crust = .true.
       call meshfem3D_models_get3Dcrust_val(IREGION_CRUST_MANTLE,xmesh,ymesh,zmesh,r, &
            vpv,vph,vsv,vsh,rho,eta_aniso,dvp, &
            c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
            c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
            elem_in_crust,moho)
       print *, 'moho depth [km]:',moho*R_EARTH_KM, (1.0d0-moho)*R_EARTH_KM,  vpv,vph,vsv,vsh,rho,eta_aniso,dvp !, 'moho radius:',1.0d0 - moho, 'in km: ',(1.0d0-moho)*R_EARTH_KM
    endif
    if(TOPOGRAPHY .or. OCEANS) then
       lat=(PI/2.0d0-theta)*180.0d0/PI
       lon=phi*180.0d0/PI
       if(lon>180.0d0) lon=lon-360.0d0
       print *,'get_topo_bathy(lat,lon,elevation,ibathy_topo',lat,lon,elevation
       call get_topo_bathy(lat,lon,elevation,ibathy_topo)
       print *,'get_topo_bathy(lat,lon,elevation,ibathy_topo',lat,lon,elevation
       print *, 'elevation [km]:',elevation/1000.0d0!, 'surface radius:',1.0d0 + elevation /R_EARTH
    endif

    do ilayer = 1,NUMBER_OF_MESH_LAYERS  ! loop over all layers

       if(ilayer == 1) then
         rmin = 0.0d0
         rmax = rmins(NUMBER_OF_MESH_LAYERS-1)
         idoubling = IFLAG_INNER_CORE_NORMAL
       else
         rmin = rmins(NUMBER_OF_MESH_LAYERS-ilayer+1)
         rmax = rmaxs(NUMBER_OF_MESH_LAYERS-ilayer+1)
         idoubling = doubling_index(NUMBER_OF_MESH_LAYERS-ilayer+1)
       endif

!  make sure that the Moho discontinuity is at the real moho
      if(CRUSTAL) then
        if(rmin == RMOHO_FICTITIOUS_IN_MESHER/R_EARTH) then
          rmin = 1.0d0 - moho
!          write(*,*) 'rmin == RMOHO',iline
        endif
        if(rmax == RMOHO_FICTITIOUS_IN_MESHER/R_EARTH) rmax = 1.0d0 - moho
      endif


      if(rmin == rmax_last) then  !!!! this means that we have just jumped between layers

       ! write values every 10 km in the deep earth and every 1 km in the shallow earth
       if(rmin>(R_EARTH_KM-100.0d0)/R_EARTH_KM) then
         delta = 1.0d0/R_EARTH_KM
       else
         delta = 10.0d0/R_EARTH_KM
       endif

       if((TOPOGRAPHY .or. OCEANS) .or. ((.not. CRUSTAL) .and. (ROCEAN < R_EARTH))) then
        if(rmax == 1.0d0) rmax = ROCEAN /R_EARTH
       endif

       rmax_last = rmax
       nit = floor((rmax - rmin)/delta) + 1
       do islice = 1,nit+1
         r = rmin + (islice-1)*delta
         if(rmin == RICB/R_EARTH .and. islice == 1) iline_icb = iline
         if(rmin == RCMB/R_EARTH .and. islice == 1) iline_cmb = iline
         if(CRUSTAL) then
           if(rmin == (1.0d0 - moho) .and. islice == 1) then
              iline_moho = iline
           endif
         else
           if(rmin == RMOHO/R_EARTH .and. islice == 1) iline_moho = iline
         endif

        ! initializes values
        rho = 0.d0
        vpv = 0.d0
        vph = 0.d0
        vsv = 0.d0
        vsh = 0.d0
        eta_aniso = 0.d0
        c11 = 0.d0
        c12 = 0.d0
        c13 = 0.d0
        c14 = 0.d0
        c15 = 0.d0
        c16 = 0.d0
        c22 = 0.d0
        c23 = 0.d0
        c24 = 0.d0
        c25 = 0.d0
        c26 = 0.d0
        c33 = 0.d0
        c34 = 0.d0
        c35 = 0.d0
        c36 = 0.d0
        c44 = 0.d0
        c45 = 0.d0
        c46 = 0.d0
        c55 = 0.d0
        c56 = 0.d0
        c66 = 0.d0
        Qmu = 0.d0
        Qkappa = 0.d0 ! not used, not stored so far...
        tau_e(:) = 0.d0
        dvp = 0.d0

        ! make sure we are within the right shell in PREM to honor discontinuities
        ! use small geometrical tolerance
        r_prem = r
        if(r <= rmin*1.000001d0) r_prem = rmin*1.000001d0
        if(r >= rmax*0.999999d0) r_prem = rmax*0.999999d0

        ! convert from rthetaphi to xyz to use in function calls.

        call rthetaphi_2_xyz(xmesh,ymesh,zmesh,r_prem,theta,phi)

!!!!-------------------------------------
        !!start GET_MODEL

        ! checks r_prem,rmin/rmax and assigned idoubling
        call get_model_check_idoubling(r_prem,xmesh,ymesh,zmesh,rmin,rmax,idoubling, &
                            RICB,RCMB,RTOPDDOUBLEPRIME, &
                            R220,R670,myrank)

        ! gets reference model values: rho,vpv,vph,vsv,vsh and eta_aniso
        call meshfem3D_models_get1D_val(myrank,iregion_code,idoubling, &
                              r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                              Qkappa,Qmu,RICB,RCMB, &
                              RTOPDDOUBLEPRIME,R80,R120,R220,R400,R600,R670,R771, &
                              RMOHO,RMIDDLE_CRUST,ROCEAN)

        ! gets the 3-D model parameters for the mantle
        call meshfem3D_models_get3Dmntl_val(iregion_code,r_prem,rho,dvp,&
                              vpv,vph,vsv,vsh,eta_aniso, &
                              RCMB,R670,RMOHO, &
                              xmesh,ymesh,zmesh,r, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                              c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

        ! gets the 3-D crustal model
        if( CRUSTAL ) then
          if( .not. elem_in_mantle) &
            call meshfem3D_models_get3Dcrust_val(iregion_code,xmesh,ymesh,zmesh,r, &
                              vpv,vph,vsv,vsh,rho,eta_aniso,dvp, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                              c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                              elem_in_crust,moho)
        endif

!!! VH  commented out following two lines from get_model
!        ! overwrites with tomographic model values (from iteration step) here, given at all GLL points
!        call meshfem3D_models_impose_val(vpv,vph,vsv,vsh,rho,dvp,eta_aniso,&
!                                        myrank,iregion_code,ispec,i,j,k)

        ! checks vpv: if close to zero then there is probably an error
        if( vpv < TINYVAL ) then
          print*,'error vpv: ',vpv,vph,vsv,vsh,rho
          print*,'radius:',r*R_EARTH_KM
          call exit_mpi(myrank,'error get_model values')
        endif

        !> Hejun
        ! New Attenuation assignment
        ! Define 3D and 1D Attenuation after moho stretch
        ! and before TOPOGRAPHY/ELLIPCITY
        !
        !note:  only Qmu attenuation considered, Qkappa attenuation not used so far...
        if( ATTENUATION ) then
          call meshfem3D_models_getatten_val(idoubling,xmesh,ymesh,zmesh,r_prem, &
                              tau_e,tau_s,T_c_source, &
                              moho,Qmu,Qkappa,elem_in_crust) ! R80
        endif

!END GET_MODEL
!!!!--------------------------
!   make sure that the first point of the profile is at zero and the last is at the surface
      if(islice == 1) then
         r_prem = rmin
      else if(islice == nit+1) then
         r_prem = rmax
      endif

!---
! add_topography
if((THREE_D_MODEL/=0 .or. TOPOGRAPHY ) .and. &
   (idoubling == IFLAG_CRUST .or. idoubling == IFLAG_80_MOHO .or. idoubling == IFLAG_220_80)) then

   !print *, 'adding topography.  elevation: ',elevation
   gamma = (r_prem - R220/R_EARTH) / (R_UNIT_SPHERE - R220/R_EARTH)
   if(gamma < -0.02 .or. gamma > 1.02) print *, 'incorrect value of gamma for topograpy'
!   print *,'rprem before: ',r_prem*R_EARTH
   r_prem = r_prem*(ONE + gamma * (elevation/R_EARTH) /r_prem)
!   print *,'r_prem after: ',r_prem*R_EARTH
else if((.not. CRUSTAL) .and. (ROCEAN < R_EARTH)) then
   r_prem = ROCEAN/R_EARTH
endif
! end add_topography
!--------

! scale values read from routines back to true values
      scaleval = dsqrt(PI*GRAV*RHOAV)
      rho = rho*RHOAV/1000.0d0
      vpv  = vpv*R_EARTH*scaleval/1000.0d0
      vph  = vph*R_EARTH*scaleval/1000.0d0
      vsv  = vsv*R_EARTH*scaleval/1000.0d0
      vsh  = vsh*R_EARTH*scaleval/1000.0d0

       iline = iline +1
!       write(57,'(i3,11F10.4)') &
!            iline,sngl(rmin*R_EARTH_KM),sngl(rmax*R_EARTH_KM),sngl(r_prem*R_EARTH_KM),sngl(r*R_EARTH_KM), &
!            sngl(vpv),sngl(vph),sngl(vsv),sngl(vsh),sngl(rho),sngl(eta_aniso),sngl(Qmu)
! finally write the values obtained at the given depth to file

       write(57,'(F8.0,7F9.2,F9.5)') &
           sngl(r_prem*R_EARTH),sngl(rho*1000.d0),sngl(vpv*1000.d0),sngl(vsv*1000.d0), &
            sngl(Qkappa),sngl(Qmu),sngl(vph*1000.d0),sngl(vsh*1000.d0),sngl(eta_aniso)
    enddo !islice
   endif !rmin == rmax_last
  enddo !ilayer

!--------
!!  This part adds the ocean to profile where needed

  if((OCEANS .and. elevation < -500.0) .or. ((.not. CRUSTAL) .and. (ROCEAN < R_EARTH))) then
     iline_ocean = iline
     if(OCEANS .and. elevation < -500.0) then
!        iline_ocean = iline
        nlayers_ocean = floor(-elevation/500.0d0)
     else if((.not. CRUSTAL) .and. (ROCEAN < R_EARTH)) then
        nlayers_ocean = floor((R_EARTH - ROCEAN)/500.0d0)
     endif

     do ilayers_ocean=0,nlayers_ocean
        r_ocean = r_prem + ilayers_ocean*0.5d0/R_EARTH_KM
        write(57,'(F8.0,7F9.2,F9.5)') &
            sngl(r_ocean*R_EARTH),1020.0,1450.,0.0,57822.5,0.0,1450.0,0.0,1.0
        iline = iline +1
     enddo
     write(57,'(F8.0,7F9.2,F9.5)') &
       sngl(1.0d0*R_EARTH),1020.0,1450.,0.0,57822.5,0.0,1450.,0.0,1.0
     iline = iline+1
     write(57,'(5i5)') iline,iline_icb,iline_cmb,iline_moho,iline_ocean
  else
     write(57,'(5i5)') iline,iline_icb,iline_cmb,iline_moho
  endif
!------end adding ocean

  enddo !sum over phi
  enddo !sum over theta
!!!!!!!!

  call MPI_FINALIZE(ier)

  end program xwrite_profile

