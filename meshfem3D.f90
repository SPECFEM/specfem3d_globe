!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================
!
! Copyright September 2002, by the California Institute of Technology.
! ALL RIGHTS RESERVED. United States Government Sponsorship Acknowledged.
!
! Any commercial use must be negotiated with the Office of Technology
! Transfer at the California Institute of Technology. This software may be
! subject to U.S. export control laws and regulations. By accepting
! this software, the user agrees to comply with all applicable U.S. export laws
! and regulations, including the International Traffic and Arms Regulations,
! 22 C.F.R. 120-130 and the Export Administration Regulations,
! 15 C.F.R. 730-744. User has the responsibility to obtain export licenses,
! or other export authority as may be required before exporting such
! information to foreign countries or providing access to foreign nationals.
! In no event shall the California Institute of Technology be liable to any
! party for direct, indirect, special, incidental or consequential damages,
! including lost profits, arising out of the use of this software and its
! documentation, even if the California Institute of Technology has been
! advised of the possibility of such damage.
!
! The California Institute of Technology specifically disclaims any
! warranties, including the implied warranties or merchantability and fitness
! for a particular purpose. The software and documentation provided hereunder
! is on an "as is" basis, and the California Institute of Technology has no
! obligations to provide maintenance, support, updates, enhancements or
! modifications.
!

  program meshfem3D

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

!=====================================================================!
!                                                                     !
!  meshfem3D produces a spectral element grid for the Earth.          !
!  This is accomplished based upon a mapping of the face of a cube    !
!  to a portion of the sphere (Ronchi et al., The Cubed Sphere)       !
!  Grid density is decreased by a factor of two                       !
!  below the moho, below the 670 and at the bottom of the outer core  !
!                                                                     !
!=====================================================================!
!
! If you use this code for your own research, please send an email
! to Jeroen Tromp <jtromp@gps.caltech.edu> for information, and cite:
!
! @ARTICLE{KoRiTr02,
! author={D. Komatitsch and J. Ritsema and J. Tromp},
! year=2002,
! title={The Spectral-Element Method, {B}eowulf Computing, and Global Seismology},
! journal={Science},
! volume=298,
! pages={1737-1742}}
!
! @ARTICLE{KoTr02a,
! author={D. Komatitsch and J. Tromp},
! year=2002,
! title={Spectral-Element Simulations of Global Seismic Wave Propagation{-I. V}alidation},
! journal={Geophys. J. Int.},
! volume=149,
! pages={390-412}}
!
! @ARTICLE{KoTr02b,
! author={D. Komatitsch and J. Tromp},
! year=2002,
! title={Spectral-Element Simulations of Global Seismic Wave Propagation{-II. 3-D} Models, Oceans, Rotation, and Self-Gravitation},
! journal={Geophys. J. Int.},
! volume=150,
! pages={303-318}}
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! year=1999,
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! volume=139,
! pages={806-822}}
!
! If you use 3-D model S20RTS, please cite
!
! @ARTICLE{RiVa00,
! author={J. Ritsema and H. J. {Van Heijst}},
! year=2000,
! title={Seismic imaging of structural heterogeneity in {E}arth's mantle: Evidence for large-scale mantle flow},
! journal={Science Progress},
! volume=83,
! pages={243-259}}
!
! Evolution of the code:
! ---------------------
!
! MPI v. 3.3 Dimitri Komatitsch, Caltech, September 2002:
!      flexible mesh doubling in outer core, inlined code, OpenDX support
! MPI v. 3.2 Jeroen Tromp, Caltech, July 2002:
!      multiple sources and flexible PREM reading
! MPI v. 3.1 Dimitri Komatitsch, Caltech, June 2002:
!      vectorized loops in solver and merged central cube
! MPI v. 3.0 Dimitri Komatitsch and Jeroen Tromp, Caltech, May 2002:
!   ported to SGI and Compaq, double precision solver, more general anisotropy
! MPI v. 2.3 Dimitri Komatitsch and Jeroen Tromp, Caltech, August 2001:
!                       gravity, rotation, oceans and 3-D models
! MPI v. 2.2 Dimitri Komatitsch and Jeroen Tromp, Caltech, March 2001:
!                       final MPI package
! MPI v. 2.0 Dimitri Komatitsch, Harvard, January 2000: MPI code for the globe
! MPI v. 1.0 Dimitri Komatitsch, Mexico, June 1999: first MPI code for a chunk
! Jeroen Tromp, Harvard, July 1998: first chunk solver using OpenMP on Sun
! Dimitri Komatitsch, IPG Paris, December 1996: first 3-D solver for the CM5
!

! correct number of spectral elements in each block depending on chunk type
  integer nspec,nspec_aniso,nspec_aniso_mantle,nspec_aniso_mantle_all,npointot

! meshing parameters
  double precision, dimension(:), allocatable :: rns,rmins,rmaxs

! auxiliary variables to generate the mesh
  integer ix,iy,ir

  double precision xin,etan,xi,eta,rn
  double precision x,y,gamma,rgt,rgb,rg_icb
  double precision x_top,y_top,z_top
  double precision x_bot,y_bot,z_bot
  double precision x_icb,y_icb,z_icb
  double precision x_central_cube,y_central_cube,z_central_cube

  double precision, dimension(:,:,:), allocatable :: xgrid,ygrid,zgrid

! parameters needed to store the radii of the grid points
! in the spherically symmetric Earth
  integer, dimension(:), allocatable :: idoubling
  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
  integer myrank,sizeprocs,ier

! check area and volume of the final mesh
  double precision area_local_bottom,area_total_bottom
  double precision area_local_top,area_total_top
  double precision volume_local,volume_total,volume_total_region

  integer iprocnum,npx,npy

! for loop on all the slices
  integer iregion_code
  integer iproc_xi,iproc_eta,ichunk

! for crustal model Crust2.0
  double precision thlr(NKEYS_CRUST,NLAYERS_CRUST),velocp(NKEYS_CRUST,NLAYERS_CRUST)
  double precision velocs(NKEYS_CRUST,NLAYERS_CRUST),dens(NKEYS_CRUST,NLAYERS_CRUST)
  character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST),code(NKEYS_CRUST)

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

! for Montagner model
  integer npar1_montagner
  double precision beta_montagner(14,34,37,73),pro_montagner(47)

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

! mantle coefficients for 3D models
  double precision dvs_a(0:NK,0:NS,0:NS),dvs_b(0:NK,0:NS,0:NS)
  double precision dvp_a(0:NK,0:NS,0:NS),dvp_b(0:NK,0:NS,0:NS)
  double precision spknt(NK+1),qq0(NK+1,NK+1),qq(3,NK+1,NK+1)

! for some statistics for the mesh
  integer numelem_crust_mantle,numelem_outer_core,numelem_inner_core
  integer numelem_total

! timer MPI
  double precision time_start,tCPU

! addressing for all the slices
  integer, dimension(:), allocatable :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
  integer, dimension(:,:,:), allocatable :: addressing

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST,NER_220_MOHO,NER_400_220, &
             NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
             NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
             NEX_ETA,NEX_XI,NER_DOUBLING_OUTER_CORE, &
             NPROC_ETA,NPROC_XI,NSEIS,NSTEP

  double precision DT

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY, &
             GRAVITY,ONE_CRUST,ROTATION, &
             THREE_D,TOPOGRAPHY,ATTENUATION,OCEANS
  integer NSOURCES,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB
  double precision RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC

  character(len=150) LOCAL_PATH

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_AB,NSPEC_AC,NSPEC_BC, &
               NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
               nglob_AB,nglob_AC,nglob_BC

! ************** PROGRAM STARTS HERE **************

! initialize the MPI communicator and start the NPROCTOT MPI processes.
! sizeprocs returns number of processes started (should be equal to NPROCTOT).
! myrank is the rank of each process, between 0 and NPROCTOT-1.
! as usual in MPI, process 0 is in charge of coordinating everything
! and also takes care of the main output
! do not create anything for the inner core here, will be done in solver
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

! open main output file, only written to by process 0
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file='OUTPUT_FILES/output_mesher.txt',status='unknown')

! get MPI starting time
  time_start = MPI_WTIME()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '****************************'
    write(IMAIN,*) '*** Specfem3D MPI Mesher ***'
    write(IMAIN,*) '****************************'
    write(IMAIN,*)
  endif

! read the parameter file
  call read_parameter_file(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST,NER_220_MOHO,NER_400_220, &
        NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
        NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB,NER_DOUBLING_OUTER_CORE, &
        NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSEIS,NSTEP, &
        DT,TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,CRUSTAL,OCEANS,ELLIPTICITY, &
        GRAVITY,ONE_CRUST,ATTENUATION, &
        ROTATION,THREE_D,TOPOGRAPHY,LOCAL_PATH,NSOURCES,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC)

! compute other parameters based upon values read
  call compute_parameters(NER_CRUST,NER_220_MOHO,NER_400_220, &
      NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
      NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
      NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB, &
      NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NSPEC_AB,NSPEC_AC,NSPEC_BC, &
      NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
      NGLOB_AB,NGLOB_AC,NGLOB_BC,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB)

! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROCTOT) call exit_MPI(myrank,'wrong number of MPI processes')

! dynamic allocation of mesh arrays
  allocate(rns(0:2*NER))
  allocate(rmins(0:2*NER))
  allocate(rmaxs(0:2*NER))

  allocate(xgrid(0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA,0:2*NER))
  allocate(ygrid(0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA,0:2*NER))
  allocate(zgrid(0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA,0:2*NER))

  allocate(addressing(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1))
  allocate(ichunk_slice(0:NPROCTOT-1))
  allocate(iproc_xi_slice(0:NPROCTOT-1))
  allocate(iproc_eta_slice(0:NPROCTOT-1))

! clear arrays
  xgrid(:,:,:) = 0.
  ygrid(:,:,:) = 0.
  zgrid(:,:,:) = 0.

  addressing(:,:,:) = 0
  ichunk_slice(:) = 0
  iproc_xi_slice(:) = 0
  iproc_eta_slice(:) = 0

! loop on all the chunks to create global slice addressing for solver
  if(myrank == 0) then
    open(unit=IOUT,file='OUTPUT_FILES/addressing.txt',status='unknown')
    write(IMAIN,*) 'creating global slice addressing'
    write(IMAIN,*)
  endif
  do ichunk = 1,NCHUNKS
    do iproc_eta=0,NPROC_ETA-1
      do iproc_xi=0,NPROC_XI-1
        iprocnum = (ichunk-1)*NPROC + iproc_eta * NPROC_XI + iproc_xi
        addressing(ichunk,iproc_xi,iproc_eta) = iprocnum
        ichunk_slice(iprocnum) = ichunk
        iproc_xi_slice(iprocnum) = iproc_xi
        iproc_eta_slice(iprocnum) = iproc_eta
        if(myrank == 0) write(IOUT,*) iprocnum,ichunk,iproc_xi,iproc_eta
      enddo
    enddo
  enddo
  if(myrank == 0) close(IOUT)

  if(myrank == 0) then
    write(IMAIN,*) 'This is process ',myrank
    write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
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
    write(IMAIN,*)
  endif

! check that reals are either 4 or 8 bytes
  if(CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) call exit_MPI(myrank,'wrong size of CUSTOM_REAL for reals')

! check that the parameter file is correct
  if(NGNOD /= 27) call exit_MPI(myrank,'number of control nodes must be 27')
  if(NGNOD == 27 .and. NGNOD2D /= 9) call exit_MPI(myrank,'elements with 27 points should have NGNOD2D = 9')

! for the number of standard linear solids for attenuation
  if(N_SLS /= 3) call exit_MPI(myrank,'number of SLS must be 3')

! check number of slices in each direction
  if(NCHUNKS < 1) call exit_MPI(myrank,'must have at least one chunk')
  if(NPROC_XI < 1) call exit_MPI(myrank,'NPROC_XI must be at least 1')
  if(NPROC_ETA < 1) call exit_MPI(myrank,'NPROC_ETA must be at least 1')

! correct topology if more than one chunk
  if(NCHUNKS /= 1 .and. NPROC_XI /= NPROC_ETA) &
    call exit_MPI(myrank,'must have NPROC_XI=NPROC_ETA for more than one chunk')

! check number of elements per processor is the same in both directions
  if(NCHUNKS /= 1 .and. NEX_PER_PROC_XI /= NEX_PER_PROC_ETA) &
    call exit_MPI(myrank,'not the same number of elements per processor in both directions')

! check number of chunks
  if(NCHUNKS /= 1 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) &
    call exit_MPI(myrank,'only one, three or six chunks can be meshed')

  if(INCLUDE_CENTRAL_CUBE .and. NCHUNKS /= 6) &
    call exit_MPI(myrank,'need six chunks to include central cube')

! check that size can be coarsened in depth twice (block size multiple of 8)
  if(mod(NEX_XI/8,NPROC_XI) /= 0) &
    call exit_MPI(myrank,'NEX_XI must be a multiple of 8*NPROC_XI')

  if(mod(NEX_ETA/8,NPROC_ETA) /= 0) &
    call exit_MPI(myrank,'NEX_ETA must be a multiple of 8*NPROC_ETA')

  if(mod(NEX_XI,8) /= 0) call exit_MPI(myrank,'NEX_XI must be a multiple of 8')

  if(mod(NEX_ETA,8) /= 0) call exit_MPI(myrank,'NEX_ETA must be a multiple of 8')

! one doubling layer in outer core (block size multiple of 16)
    if(mod(NEX_PER_PROC_XI,16) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_XI must be a multiple of 16 for outer core doubling')
    if(mod(NEX_PER_PROC_ETA,16) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_ETA must be a multiple of 16 for outer core doubling')

! check that mesh doubling can be implemented
  if(NER_220_MOHO/2<3) call exit_MPI(myrank,'NER_220_MOHO should be at least 3')
  if(NER_400_220/2<2) call exit_MPI(myrank,'NER_400_220 should be at least 2')
  if(NER_670_400/2<2) call exit_MPI(myrank,'NER_670_400 should be at least 2')
  if(NER_CMB_670<16) call exit_MPI(myrank,'NER_CMB_670 should be at least 16')
  if(NER_CENTRAL_CUBE_CMB<4) call exit_MPI(myrank,'NER_CENTRAL_CUBE_CMB should be at least 4')

  if(IASPEI .and. TRANSVERSE_ISOTROPY) call exit_MPI(myrank,'IASPEI is currently isotropic')

  if(myrank == 0) then

  write(IMAIN,*)
  if(ELLIPTICITY) then
    write(IMAIN,*) 'incorporating ellipticity'
  else
    write(IMAIN,*) 'no ellipticity'
  endif

  write(IMAIN,*)
  if(TOPOGRAPHY) then
    write(IMAIN,*) 'incorporating surface topography'
  else
    write(IMAIN,*) 'no surface topography'
  endif

  write(IMAIN,*)
  if(THREE_D) then
    write(IMAIN,*) 'incorporating 3-D lateral variations'
  else
    write(IMAIN,*) 'no 3-D lateral variations'
  endif

  write(IMAIN,*)
  if(CRUSTAL) then
    write(IMAIN,*) 'incorporating crustal variations'
  else
    write(IMAIN,*) 'no crustal variations'
  endif

  write(IMAIN,*)
  if(ONE_CRUST) then
    write(IMAIN,*) 'using one layer only in PREM crust'
  else
    write(IMAIN,*) 'using real PREM crust with two layers'
  endif

  write(IMAIN,*)
  if(GRAVITY) then
    write(IMAIN,*) 'incorporating self-gravitation (Cowling approximation)'
  else
    write(IMAIN,*) 'no self-gravitation'
  endif

  write(IMAIN,*)
  if(ROTATION) then
    write(IMAIN,*) 'incorporating rotation'
  else
    write(IMAIN,*) 'no rotation'
  endif

  write(IMAIN,*)
  if(TRANSVERSE_ISOTROPY) then
    write(IMAIN,*) 'incorporating anisotropy'
  else
    write(IMAIN,*) 'no anisotropy'
  endif

  write(IMAIN,*)
  if(ATTENUATION) then
    write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
  else
    write(IMAIN,*) 'no attenuation'
  endif

  write(IMAIN,*)
  if(OCEANS) then
    write(IMAIN,*) 'incorporating the oceans using equivalent load'
  else
    write(IMAIN,*) 'no oceans'
  endif

  write(IMAIN,*)

  endif

  if(ELLIPTICITY) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

  if(THREE_D) call read_3D_model(dvs_a,dvs_b,dvp_a,dvp_b,spknt,qq0,qq)

! read Montagner's model
  if(ANISOTROPIC_MANTLE) call read_modanis(beta_montagner,pro_montagner,npar1_montagner)

  if(CRUSTAL) call read_3D_crustal_model(abbreviation,code,thlr,velocp,velocs,dens)

! read topography and bathymetry file
  if(TOPOGRAPHY .or. OCEANS) call read_topo_bathy_file(ibathy_topo)

! get addressing for this process
  ichunk = ichunk_slice(myrank)
  iproc_xi = iproc_xi_slice(myrank)
  iproc_eta = iproc_eta_slice(myrank)

  if(myrank == 0) then
    write(IMAIN,*) 'Reference radius of the Earth used is ',R_EARTH_KM,' km'
    write(IMAIN,*)
    write(IMAIN,*) 'Central cube is at a radius of ',R_CENTRAL_CUBE*R_EARTH/1000.d0,' km'
  endif

! create the radial mesh
  call mesh_radial(myrank,rns,rmins,rmaxs,NER,NER_TOP_CENTRAL_CUBE_ICB, &
      NER_CMB_TOPDDOUBLEPRIME,NER_TOPDDOUBLEPRIME_771, &
      NER_771_670,NER_670_600,NER_600_400,NER_400_220,NER_220_MOHO,NER_CRUST, &
      NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
      CRUSTAL,TOPOGRAPHY,ONE_CRUST)

! number of elements in each slice in a given chunk
  npx = 2*NEX_PER_PROC_XI
  npy = 2*NEX_PER_PROC_ETA

! fill the region between the central cube and the free surface
  do iy=0,npy
  do ix=0,npx

! full Earth (cubed sphere)

    xin=dble(ix)/dble(npx)
    xi= - PI_OVER_FOUR + (dble(iproc_xi)+xin)*PI_OVER_TWO/dble(NPROC_XI)

    x=dtan(xi)

    etan=dble(iy)/dble(npy)
    eta= - PI_OVER_FOUR + (dble(iproc_eta)+etan)*PI_OVER_TWO/dble(NPROC_ETA)

    y=dtan(eta)
    gamma=ONE/dsqrt(ONE+x*x+y*y)

    rgt=R_UNIT_SPHERE*gamma
    rgb=R_CENTRAL_CUBE*gamma
    rg_icb = (RICB/R_EARTH)*gamma

!   define the mesh points on the top and the bottom
!   in the six regions of the cubed shpere

    if(ichunk == CHUNK_AB) then
      x_top=-y*rgt
      y_top=x*rgt
      z_top=rgt

      x_bot=-y*rgb
      y_bot=x*rgb
      z_bot=rgb

      x_icb = -y*rg_icb
      y_icb = x*rg_icb
      z_icb = rg_icb
      x_central_cube = -y
      y_central_cube = x
      z_central_cube = 1.

    else if(ichunk == CHUNK_AB_ANTIPODE) then
      x_top=-y*rgt
      y_top=-x*rgt
      z_top=-rgt

      x_bot=-y*rgb
      y_bot=-x*rgb
      z_bot=-rgb

      x_icb = -y*rg_icb
      y_icb = -x*rg_icb
      z_icb = -rg_icb
      x_central_cube = -y
      y_central_cube = -x
      z_central_cube = -1.

    else if(ichunk == CHUNK_AC) then
      x_top=-y*rgt
      y_top=-rgt
      z_top=x*rgt

      x_bot=-y*rgb
      y_bot=-rgb
      z_bot=x*rgb

      x_icb = -y*rg_icb
      y_icb = -rg_icb
      z_icb = x*rg_icb
      x_central_cube = -y
      y_central_cube = -1.
      z_central_cube = x

    else if(ichunk == CHUNK_AC_ANTIPODE) then
      x_top=-y*rgt
      y_top=rgt
      z_top=-x*rgt

      x_bot=-y*rgb
      y_bot=rgb
      z_bot=-x*rgb

      x_icb = -y*rg_icb
      y_icb = rg_icb
      z_icb = -x*rg_icb
      x_central_cube = -y
      y_central_cube = 1.
      z_central_cube = -x

    else if(ichunk == CHUNK_BC) then
      x_top=-rgt
      y_top=y*rgt
      z_top=x*rgt

      x_bot=-rgb
      y_bot=y*rgb
      z_bot=x*rgb

      x_icb = -rg_icb
      y_icb = y*rg_icb
      z_icb = x*rg_icb
      x_central_cube = -1.
      y_central_cube = y
      z_central_cube = x

    else if(ichunk == CHUNK_BC_ANTIPODE) then
      x_top=rgt
      y_top=-y*rgt
      z_top=x*rgt

      x_bot=rgb
      y_bot=-y*rgb
      z_bot=x*rgb

      x_icb = rg_icb
      y_icb = -y*rg_icb
      z_icb = x*rg_icb
      x_central_cube = 1.
      y_central_cube = -y
      z_central_cube = x

    else
      call exit_MPI(myrank,'incorrect chunk numbering in meshfem3D')
    endif

! rescale central cube to match cubed sphere
    x_central_cube = x_central_cube * R_CENTRAL_CUBE / dsqrt(3.d0)
    y_central_cube = y_central_cube * R_CENTRAL_CUBE / dsqrt(3.d0)
    z_central_cube = z_central_cube * R_CENTRAL_CUBE / dsqrt(3.d0)

!   fill the volume
    do ir=0,2*NER
      rn=rns(ir)
      xgrid(ix,iy,ir)=x_top*rn+x_bot*(ONE-rn)
      ygrid(ix,iy,ir)=y_top*rn+y_bot*(ONE-rn)
      zgrid(ix,iy,ir)=z_top*rn+z_bot*(ONE-rn)
    enddo

! modify in the inner core to match the central cube instead of a sphere
! this mesh works because even if ellipticity and/or topography are turned on
! at this stage the reference Earth is still purely spherical
! therefore it will always perfectly match the sphere defined above
    do ir=0,2*NER_TOP_CENTRAL_CUBE_ICB
      rn = dble(ir) / dble(2*NER_TOP_CENTRAL_CUBE_ICB)
      xgrid(ix,iy,ir) = x_icb * rn + x_central_cube * (ONE-rn)
      ygrid(ix,iy,ir) = y_icb * rn + y_central_cube * (ONE-rn)
      zgrid(ix,iy,ir) = z_icb * rn + z_central_cube * (ONE-rn)
    enddo

  enddo
  enddo

! volume of the slice
  volume_total = ZERO

! number of anisotropic elements in the mantle
  nspec_aniso_mantle = 0

! make sure everybody is synchronized
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

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

! volume of bottom and top area of the slice
  volume_local = ZERO
  area_local_bottom = ZERO
  area_local_top = ZERO

! check chunk number and assign theoretical number of elements
  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    nspec = NSPEC_AB(iregion_code)
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    nspec = NSPEC_AC(iregion_code)
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    nspec = NSPEC_BC(iregion_code)
  else
    call exit_MPI(myrank,'incorrect chunk number')
  endif

! compute maximum number of points
  npointot = nspec * NGLLCUBE

! use dynamic allocation to allocate memory for arrays
  allocate(idoubling(nspec))
  allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(ystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(zstore(NGLLX,NGLLY,NGLLZ,nspec))

! create all the regions of the mesh
  call create_regions_mesh(iregion_code,xgrid,ygrid,zgrid,ibool,idoubling, &
         xstore,ystore,zstore,npx,npy,rmins,rmaxs, &
         iproc_xi,iproc_eta,ichunk,nspec,nspec_aniso, &
         volume_local,area_local_bottom,area_local_top, &
         beta_montagner,pro_montagner,npar1_montagner, &
         dvs_a,dvs_b,dvp_a,dvp_b,spknt,qq0,qq,abbreviation,code,thlr,velocp,velocs,dens, &
         nspl,rspl,espl,espl2, &
         nglob_AB(iregion_code),nglob_AC(iregion_code),nglob_BC(iregion_code),npointot, &
         NER,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
         NER_TOP_CENTRAL_CUBE_ICB,NER_CENTRAL_CUBE_CMB,NER_CMB_670,NER_670_400,NER_400_220, &
         NER_220_MOHO,NER_CRUST,NER_DOUBLING_OUTER_CORE, &
         NSPEC2DMAX_XMIN_XMAX(iregion_code), &
         NSPEC2DMAX_YMIN_YMAX(iregion_code),NSPEC2D_BOTTOM(iregion_code),NSPEC2D_TOP(iregion_code), &
         ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE, &
         ANISOTROPIC_INNER_CORE,THREE_D,CRUSTAL,ONE_CRUST, &
         NPROC_XI,NPROC_ETA,NSPEC2D_A_XI(iregion_code),NSPEC2D_B_XI(iregion_code),NSPEC2D_C_XI(iregion_code), &
         NSPEC2D_A_ETA(iregion_code),NSPEC2D_B_ETA(iregion_code), &
         NSPEC2D_C_ETA(iregion_code),NSPEC1D_RADIAL(iregion_code),NPOIN1D_RADIAL(iregion_code), &
         myrank,LOCAL_PATH,OCEANS,ibathy_topo,NER_ICB_BOTTOMDBL)

! store number of anisotropic elements found in the mantle
  if(nspec_aniso /= 0 .and. iregion_code /= IREGION_CRUST_MANTLE) &
    call exit_MPI(myrank,'found anisotropic elements outside of the mantle')

  if(iregion_code == IREGION_CRUST_MANTLE .and. nspec_aniso == 0) &
    call exit_MPI(myrank,'found no anisotropic elements in the mantle')

  if(iregion_code == IREGION_CRUST_MANTLE) nspec_aniso_mantle = nspec_aniso

! use MPI reduction to compute total area and volume
  volume_total_region = ZERO
  area_total_bottom   = ZERO
  area_total_top   = ZERO
  call MPI_REDUCE(area_local_bottom,area_total_bottom,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                          MPI_COMM_WORLD,ier)
  call MPI_REDUCE(area_local_top,area_total_top,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                          MPI_COMM_WORLD,ier)
  call MPI_REDUCE(volume_local,volume_total_region,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                          MPI_COMM_WORLD,ier)

  if(myrank == 0) then
!   sum volume over all the regions
    volume_total = volume_total + volume_total_region

! if only one chunk, create dummy values for array sizes
    if(NCHUNKS == 1) then
      nglob_AC(iregion_code) = nglob_AB(iregion_code)
      nglob_BC(iregion_code) = nglob_AB(iregion_code)
    endif

!   check volume of chunk, and bottom and top area

      write(IMAIN,*)
      write(IMAIN,*) '   calculated top area: ',area_total_top

! compare to exact theoretical value
    if(.not. TOPOGRAPHY) then
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

      write(IMAIN,*)
      write(IMAIN,*) 'calculated bottom area: ',area_total_bottom

! compare to exact theoretical value
  if(.not. TOPOGRAPHY) then

      select case(iregion_code)

        case(IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2

        case(IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2

        case(IREGION_INNER_CORE)
          write(IMAIN,*) '            exact area (central cube): ',dble(NCHUNKS)*(2.*R_CENTRAL_CUBE/sqrt(3.))**2

        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
  endif

  endif

! create chunk buffers if more than one chunk

  if(NCHUNKS > 1) then
    call create_chunk_buffers(iregion_code,nspec,ibool,idoubling,xstore,ystore,zstore, &
      nglob_AB(iregion_code),nglob_AC(iregion_code),nglob_BC(iregion_code), &
      NSPEC2DMAX_XMIN_XMAX(iregion_code),NSPEC2DMAX_YMIN_YMAX(iregion_code), &
      NPROC_XI,NPROC_ETA,NPROC,NPROCTOT,NPOIN1D_RADIAL(iregion_code), &
      NPOIN2DMAX_XMIN_XMAX(iregion_code),NPOIN2DMAX_YMIN_YMAX(iregion_code), &
      NSPEC_AB(iregion_code),NSPEC_AC(iregion_code),NSPEC_BC(iregion_code), &
      myrank,LOCAL_PATH,addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice)
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

! make sure everybody is synchronized
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

! end of loop on all the regions
  enddo

! compute the maximum number of anisotropic elements found in all the slices
  call MPI_ALLREDUCE(nspec_aniso_mantle,nspec_aniso_mantle_all,1,MPI_INTEGER, &
              MPI_MAX,MPI_COMM_WORLD,ier)

  if(myrank == 0) then
! check volume of chunk
      write(IMAIN,*)
      write(IMAIN,*) 'calculated volume: ',volume_total
    if(.not. TOPOGRAPHY) then
! take the central cube into account
! it is counted 6 times because of the fictitious elements
      if(INCLUDE_CENTRAL_CUBE) then
        write(IMAIN,*) '     exact volume: ', &
          dble(NCHUNKS)*((4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3) + 5.*(2.*R_CENTRAL_CUBE/sqrt(3.))**3) / 6.d0
      else
        write(IMAIN,*) '     exact volume: ', &
          dble(NCHUNKS)*((4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3) - (2.*R_CENTRAL_CUBE/sqrt(3.))**3) / 6.d0
      endif
    endif
  endif

!--- print number of points and elements in the mesh for each region
!--- in chunks of type AB

  if(myrank == 0) then

    numelem_crust_mantle = NSPEC_AB(IREGION_CRUST_MANTLE)
    numelem_outer_core = NSPEC_AB(IREGION_OUTER_CORE)
    numelem_inner_core = NSPEC_AB(IREGION_INNER_CORE)

    numelem_total = numelem_crust_mantle + numelem_outer_core + numelem_inner_core

  write(IMAIN,*)
  write(IMAIN,*) 'Repartition of elements in regions:'
  write(IMAIN,*) '----------------------------------'
  write(IMAIN,*)
  write(IMAIN,*) 'total number of elements in each slice in type AB chunk: ',numelem_total
  write(IMAIN,*)
  write(IMAIN,*) ' - crust and mantle: ',sngl(100.d0*dble(numelem_crust_mantle)/dble(numelem_total)),' %'
  write(IMAIN,*) ' - outer core: ',sngl(100.d0*dble(numelem_outer_core)/dble(numelem_total)),' %'
  write(IMAIN,*) ' - inner core: ',sngl(100.d0*dble(numelem_inner_core)/dble(numelem_total)),' %'
  write(IMAIN,*)
  write(IMAIN,*) 'for some mesh statistics, see comments in file OUTPUT_FILES/values_from_mesher.h'
  write(IMAIN,*)

! load balancing with respect to chunk of type BC (largest number of elements)
  if(NCHUNKS > 1) then
    write(IMAIN,*) 'Load balancing AB blocks = ',100.d0*dble(sum(nspec_AB))/dble(sum(nspec_BC)),' %'
    write(IMAIN,*) 'Load balancing AC blocks = ',100.d0*dble(sum(nspec_AC))/dble(sum(nspec_BC)),' %'
    write(IMAIN,*) 'Load balancing BC blocks = 100 % by definition'
  else
    write(IMAIN,*) 'Load balancing = 100 % by definition for one chunk'
  endif
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

! create include file for the solver
  call save_header_file(NSPEC_AB,NSPEC_AC,NSPEC_BC,nglob_AB,nglob_AC,nglob_BC, &
        NEX_XI,NEX_ETA,nspec_aniso_mantle_all,NPROC,NPROCTOT, &
        TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE, &
        ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION)

  endif   ! end of section executed by main process only

! deallocate arrays used for mesh generation
  deallocate(rns)
  deallocate(rmins)
  deallocate(rmaxs)

  deallocate(xgrid)
  deallocate(ygrid)
  deallocate(zgrid)

  deallocate(addressing)
  deallocate(ichunk_slice)
  deallocate(iproc_xi_slice)
  deallocate(iproc_eta_slice)

! elapsed time since beginning of mesh generation
  if(myrank == 0) then
    tCPU = MPI_WTIME() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for mesh generation and buffer creation in seconds = ',tCPU
    write(IMAIN,*) 'End of mesh generation'
    write(IMAIN,*)
  endif

! close main output file
  if(myrank == 0) then
    write(IMAIN,*) 'done'
    write(IMAIN,*)
    close(IMAIN)
  endif

! synchronize all the processes to make sure everybody has finished
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

  end program meshfem3D

