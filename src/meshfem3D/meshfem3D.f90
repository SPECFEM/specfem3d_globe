!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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
!
! United States and French Government Sponsorship Acknowledged.

  program xmeshfem3D

  use meshfem3D_models_par
  use mpi

  implicit none

  include "precision.h"

!=====================================================================!
!                                                                     !
!  meshfem3D produces a spectral element grid for the Earth.          !
!  This is accomplished based upon a mapping of the face of a cube    !
!  to a portion of the sphere (Ronchi et al., The Cubed Sphere).      !
!  Grid density is decreased by a factor of two                       !
!  three times in the radial direction.                               !
!                                                                     !
!=====================================================================!
!
! If you use this code for your own research, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! @ARTICLE{PeKoLuMaLeCaLeMaLiBlNiBaTr11,
! author = {Daniel Peter and Dimitri Komatitsch and Yang Luo and Roland Martin
!     and Nicolas {Le Goff} and Emanuele Casarotti and Pieyre {Le Loher}
!     and Federica Magnoni and Qinya Liu and C\'eline Blitz and Tarje Nissen-Meyer
!     and Piero Basini and Jeroen Tromp},
! title = {Forward and adjoint simulations of seismic wave propagation on fully
!     unstructured hexahedral meshes},
! journal={Geophys. J. Int.},
! year = {2011},
! volume = {186},
! pages = {721-739},
! number = {2},
! doi = {10.1111/j.1365-246X.2011.05044.x}}
!
! or
!
! @INCOLLECTION{ChKoViCaVaFe07,
! author = {Emmanuel Chaljub and Dimitri Komatitsch and Jean-Pierre Vilotte and
! Yann Capdeville and Bernard Valette and Gaetano Festa},
! title = {Spectral Element Analysis in Seismology},
! booktitle = {Advances in Wave Propagation in Heterogeneous Media},
! publisher = {Elsevier - Academic Press},
! year = {2007},
! editor = {Ru-Shan Wu and Val\'erie Maupin},
! volume = {48},
! series = {Advances in Geophysics},
! pages = {365-419}}
!
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! year=1999,
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoRiTr02,
! author={D. Komatitsch and J. Ritsema and J. Tromp},
! year=2002,
! title={The Spectral-Element Method, {B}eowulf Computing, and Global Seismology},
! journal={Science},
! volume=298,
! number=5599,
! pages={1737-1742},
! doi={10.1126/science.1076024}}
!
! @ARTICLE{KoTr02a,
! author={D. Komatitsch and J. Tromp},
! year=2002,
! title={Spectral-Element Simulations of Global Seismic Wave Propagation{-I. V}alidation},
! journal={Geophys. J. Int.},
! volume=149,
! number=2,
! pages={390-412},
! doi={10.1046/j.1365-246X.2002.01653.x}}
!
! @ARTICLE{KoTr02b,
! author={D. Komatitsch and J. Tromp},
! year=2002,
! title={Spectral-Element Simulations of Global Seismic Wave Propagation{-II. 3-D} Models, Oceans, Rotation, and Self-Gravitation},
! journal={Geophys. J. Int.},
! volume=150,
! pages={303-318},
! number=1,
! doi={10.1046/j.1365-246X.2002.01716.x}}
!
! and/or another article from http://web.univ-pau.fr/~dkomati1/publications.html
!
!
! If you use the kernel capabilities of the code, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! @ARTICLE{PeKoLuMaLeCaLeMaLiBlNiBaTr11,
! author = {Daniel Peter and Dimitri Komatitsch and Yang Luo and Roland Martin
!     and Nicolas {Le Goff} and Emanuele Casarotti and Pieyre {Le Loher}
!     and Federica Magnoni and Qinya Liu and C\'eline Blitz and Tarje Nissen-Meyer
!     and Piero Basini and Jeroen Tromp},
! title = {Forward and adjoint simulations of seismic wave propagation on fully
!     unstructured hexahedral meshes},
! journal={Geophys. J. Int.},
! year = {2011},
! volume = {186},
! pages = {721-739},
! number = {2},
! doi = {10.1111/j.1365-246X.2011.05044.x}}
!
! @ARTICLE{LiTr06,
! author={Qinya Liu and Jeroen Tromp},
! title={Finite-frequency kernels based on adjoint methods},
! journal={Bull. Seismol. Soc. Am.},
! year=2006,
! volume=96,
! number=6,
! pages={2383-2397},
! doi={10.1785/0120060041}}
!
! If you use 3-D model S20RTS, please cite:
!
! @ARTICLE{RiVa00,
! author={J. Ritsema and H. J. {Van Heijst}},
! year=2000,
! title={Seismic imaging of structural heterogeneity in {E}arth's mantle: Evidence for large-scale mantle flow},
! journal={Science Progress},
! volume=83,
! pages={243-259}}
!
! Reference frame - convention:
! ----------------------------
!
! The code uses the following convention for the reference frame:
!
!  - X axis is East
!  - Y axis is North
!  - Z axis is up
!
! Note that this convention is different from both the Aki-Richards convention
! and the Harvard CMT convention.
!
! Let us recall that the Aki-Richards convention is:
!
!  - X axis is North
!  - Y axis is East
!  - Z axis is down
!
! and that the Harvard CMT convention is:
!
!  - X axis is South
!  - Y axis is East
!  - Z axis is up
!
! To report bugs or suggest improvements to the code, please send an email
! to Jeroen Tromp <jtromp AT princeton.edu> and/or use our online
! bug tracking system at http://www.geodynamics.org/roundup .
!
! Evolution of the code:
! ---------------------
!
! v. 5.1, Dimitri Komatitsch, University of Toulouse, France and Ebru Bozdag, Princeton University, USA, February 2011:
!     non blocking MPI for much better scaling on large clusters;
!     new convention for the name of seismograms, to conform to the IRIS standard;
!     new directory structure
!
! v. 5.0, many developers, February 2010:
!     new moho mesh stretching honoring crust2.0 moho depths,
!     new attenuation assignment, new SAC headers, new general crustal models,
!     faster performance due to Deville routines and enhanced loop unrolling,
!     slight changes in code structure
!
! v. 4.0 David Michea and Dimitri Komatitsch, University of Pau, France, February 2008:
!      new doubling brick in the mesh, new perfectly load-balanced mesh,
!      more flexible routines for mesh design, new inflated central cube
!      with optimized shape, far fewer mesh files saved by the mesher,
!      global arrays sorted to speed up the simulation, seismos can be
!      written by the master, one more doubling level at the bottom
!      of the outer core if needed (off by default)
!
! v. 3.6 Many people, many affiliations, September 2006:
!      adjoint and kernel calculations, fixed IASP91 model,
!      added AK135 and 1066a, fixed topography/bathymetry routine,
!      new attenuation routines, faster and better I/Os on very large
!      systems, many small improvements and bug fixes, new "configure"
!      script, new user's manual etc.
!
! v. 3.5 Dimitri Komatitsch, Brian Savage and Jeroen Tromp, Caltech, July 2004:
!      any size of chunk, 3D attenuation, case of two chunks,
!      more precise topography/bathymetry model, new Par_file structure
!
! v. 3.4 Dimitri Komatitsch and Jeroen Tromp, Caltech, August 2003:
!      merged global and regional codes, no iterations in fluid, better movies
!
! v. 3.3 Dimitri Komatitsch, Caltech, September 2002:
!      flexible mesh doubling in outer core, inlined code, OpenDX support
!
! v. 3.2 Jeroen Tromp, Caltech, July 2002:
!      multiple sources and flexible PREM reading
!
! v. 3.1 Dimitri Komatitsch, Caltech, June 2002:
!      vectorized loops in solver and merged central cube
!
! v. 3.0 Dimitri Komatitsch and Jeroen Tromp, Caltech, May 2002:
!   ported to SGI and Compaq, double precision solver, more general anisotropy
!
! v. 2.3 Dimitri Komatitsch and Jeroen Tromp, Caltech, August 2001:
!                       gravity, rotation, oceans and 3-D models
!
! v. 2.2 Dimitri Komatitsch and Jeroen Tromp, Caltech, USA, March 2001:
!                       final MPI package
!
! v. 2.0 Dimitri Komatitsch, Harvard, USA, January 2000: MPI code for the globe
!
! v. 1.0 Dimitri Komatitsch, UNAM, Mexico, June 1999: first MPI code for a chunk
!
! Jeroen Tromp and Dimitri Komatitsch, Harvard, USA, July 1998: first chunk solver using OpenMP on a Sun machine
!
! Dimitri Komatitsch, IPG Paris, France, December 1996: first 3-D solver for the CM-5 Connection Machine,
!    parallelized on 128 processors using Connection Machine Fortran
!
! From Dahlen and Tromp (1998):
! ----------------------------
!
! Gravity is approximated by solving eq (3.259) without the Phi_E' term
! The ellipsoidal reference model is that of section 14.1
! The transversely isotropic expression for PREM is that of eq (8.190)
!
! Formulation in the fluid (acoustic) outer core:
! -----------------------------------------------
!
! In case of an acoustic medium, a displacement potential Chi is used
! as in Chaljub and Valette, Geophysical Journal International, vol. 158,
! p. 131-141 (2004) and *NOT* a velocity potential as in Komatitsch and Tromp,
! Geophysical Journal International, vol. 150, p. 303-318 (2002).
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement if we ignore gravity is then: u = grad(Chi)
! (In the context of the Cowling approximation displacement is
! u = grad(rho * Chi) / rho, *not* u = grad(Chi).)
! Velocity is then: v = grad(Chi_dot)       (Chi_dot being the time derivative of Chi)
! and pressure is: p = - rho * Chi_dot_dot  (Chi_dot_dot being the time second derivative of Chi).
! The source in an acoustic element is a pressure source.
! The potential in the outer core is called displ_outer_core for simplicity.
! Its first time derivative is called veloc_outer_core.
! Its second time derivative is called accel_outer_core.


! correct number of spectral elements in each block depending on chunk type
  integer nspec_tiso,npointot

! parameters needed to store the radii of the grid points
! in the spherically symmetric Earth
  integer, dimension(:), allocatable :: idoubling
  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
  integer myrank,sizeprocs,ier

! check area and volume of the final mesh
  double precision area_local_bottom
  double precision area_local_top
  double precision volume_local,volume_total

  !integer iprocnum

! for loop on all the slices
  integer iregion_code
  integer iproc_xi,iproc_eta,ichunk

! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! for some statistics for the mesh
  integer numelem_crust_mantle,numelem_outer_core,numelem_inner_core
  integer numelem_total

! timer MPI
  double precision time_start,tCPU

! addressing for all the slices
  integer, dimension(:), allocatable :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
  integer, dimension(:,:,:), allocatable :: addressing

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,NOISE_TOMOGRAPHY, &
          ATT1,ATT2,ATT3,ATT4,ATT5

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
          RMOHO_FICTITIOUS_IN_MESHER, &
          RATIO_BY_WHICH_TO_INCREASE_IT

  logical MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,SAVE_REGULAR_KL, &
          USE_LDDRK,INCREASE_CFL_FOR_LDDRK,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY,APPROXIMATE_HESS_KL, &
          USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK,GPU_MODE,ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
          ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_AVS_DX, &
          EXACT_MASS_MATRIX_FOR_ROTATION

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
               NGLOB

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

! memory size of all the static arrays
  double precision :: static_memory_size

  integer :: ipass

  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
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
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION

! this for the different corners of the slice (which are different if the superbrick is cut)
! 1 : xi_min, eta_min
! 2 : xi_max, eta_min
! 3 : xi_max, eta_max
! 4 : xi_min, eta_max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER

! 1 -> min, 2 -> max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE

  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(MAX_NUM_REGIONS) :: NGLOB1D_RADIAL_TEMP

! this for non blocking MPI
  logical, dimension(:), allocatable :: is_on_a_slice_edge

  logical :: PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION
  integer :: NT_DUMP_ATTENUATION

! ************** PROGRAM STARTS HERE **************

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
! trivia about the programming style adopted here
!
! note 1: in general, we do not use modules in the Fortran codes. this seems to
!             be mainly a performance reason. changing the codes to adopt modules
!             will have to prove that it performs as fast as it does without now.
!
!             another reason why modules are avoided, is to make the code thread safe.
!             having different threads access the same data structure and modifying it at the same time
!             would lead to problems. passing arguments is a way to avoid such complications.
!
!             however, the mesher makes one exception here: it uses the
!             module "meshfem3D_models_par" defined in the 'meshfem3D_models.f90' file.
!             the exception is based on the fact, that when one wants to incorporate
!             a new 3D/1D velocity model, it became tedious to change so many routines hardly
!             related to any model specific need.
!
! note 2: adding a new velocity model should become easier. the module tries to help with
!             that task. basically, you would follow the comments "ADD YOUR MODEL HERE"
!             to have an idea where you will have to put some new code:
!
!                 - meshfem3D_models.f90: main file for models
!                     put your model structure into the module "meshfem3D_models_par"
!                     and add your specific routine calls to get 1D/3D/attenuation values.
!
!                 - get_model_parameters.f90:
!                     set your specific model flags and radii
!
!                 - read_compute_parameters.f90:
!                     some models need to explicitly set smaller time steps which
!                     can be done in routine rcp_set_timestep_and_layers()
!
!                 - add your model implementation into a new file named model_***.f90:
!                     in general, this file should have as first routine the model_***_broadcast() routine
!                     implemented which deals with passing the model structure to all processes.
!                     this involves reading in model specific data which is normally put in directory DATA/
!                     then follows a routine that returns the velocity values
!                     (as perturbation to the associated 1D reference model) for a given point location.
!
!             finally, in order to compile the new mesher with your new file(s),
!             you will add it to the list in the 'Makefile.in' file and run
!             `configure` to recreate a new Makefile.
!
!
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

! sizeprocs returns number of processes started (should be equal to NPROCTOT).
! myrank is the rank of each process, between 0 and NPROCTOT-1.
! as usual in MPI, process 0 is in charge of coordinating everything
! and also takes care of the main output
! do not create anything for the inner core here, will be done in solver
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! open main output file, only written to by process 0
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_mesher.txt',status='unknown')

! get MPI starting time
  time_start = MPI_WTIME()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '****************************'
    write(IMAIN,*) '*** Specfem3D MPI Mesher ***'
    write(IMAIN,*) '****************************'
    write(IMAIN,*)
  endif

  if (myrank==0) then
    ! reads the parameter file and computes additional parameters
    call read_compute_parameters(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,DT, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE,MOVIE_VOLUME_TYPE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,MOVIE_START,MOVIE_STOP, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
          ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD, &
          NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
          NSPEC,NSPEC2D_XI,NSPEC2D_ETA,NSPEC2DMAX_XMIN_XMAX, &
          NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
          NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
          ratio_sampling_array, ner, doubling_index,r_bottom,r_top,&
          this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube, &
          HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
          DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
          WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE, &
          USE_BINARY_FOR_LARGE_FILE,.false.,NOISE_TOMOGRAPHY,SAVE_REGULAR_KL, &
          PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION,NT_DUMP_ATTENUATION, &
          USE_LDDRK,INCREASE_CFL_FOR_LDDRK,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY,APPROXIMATE_HESS_KL, &
          USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK,GPU_MODE,ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
          ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_AVS_DX,RATIO_BY_WHICH_TO_INCREASE_IT, &
          ATT1,ATT2,ATT3,ATT4,ATT5,EXACT_MASS_MATRIX_FOR_ROTATION)
  endif

  ! distributes parameters from master to all processes
  call broadcast_computed_parameters(myrank,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
                NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
                NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
                NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
                NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
                NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
                NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
                MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
                DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
                CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
                RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
                R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
                MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
                RMOHO_FICTITIOUS_IN_MESHER, &
                MOVIE_SURFACE,MOVIE_VOLUME,RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
                SAVE_MESH_FILES,ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
                SAVE_ALL_SEISMOS_IN_ONE_FILE,MOVIE_COARSE,OUTPUT_SEISMOS_ASCII_TEXT, &
                OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
                ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,USE_BINARY_FOR_LARGE_FILE, &
                LOCAL_PATH,MODEL, &
                NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
                NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
                ratio_sampling_array, ner, doubling_index,r_bottom,r_top, &
                this_region_has_a_doubling,rmins,rmaxs, &
                ratio_divide_central_cube,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA, &
                REFERENCE_1D_MODEL,THREE_D_MODEL,ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &
                HONOR_1D_SPHERICAL_MOHO,CRUSTAL,ONE_CRUST,CASE_3D,TRANSVERSE_ISOTROPY, &
                ISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
                ATTENUATION,ATTENUATION_3D,ANISOTROPIC_INNER_CORE,NOISE_TOMOGRAPHY, &
                SAVE_REGULAR_KL,PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION,NT_DUMP_ATTENUATION, &
          USE_LDDRK,INCREASE_CFL_FOR_LDDRK,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY,APPROXIMATE_HESS_KL, &
          USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK,GPU_MODE,ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
          ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_AVS_DX,RATIO_BY_WHICH_TO_INCREASE_IT, &
          ATT1,ATT2,ATT3,ATT4,ATT5,EXACT_MASS_MATRIX_FOR_ROTATION)

  ! check that the code is running with the requested number of processes
  if(sizeprocs /= NPROCTOT) call exit_MPI(myrank,'wrong number of MPI processes')

  ! compute rotation matrix from Euler angles
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * PI / 180.d0
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * PI / 180.d0
  if(NCHUNKS /= 6) call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  ! dynamic allocation of mesh arrays
  allocate(addressing(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1))
  allocate(ichunk_slice(0:NPROCTOT-1))
  allocate(iproc_xi_slice(0:NPROCTOT-1))
  allocate(iproc_eta_slice(0:NPROCTOT-1))

  ! creates global slice addressing for solver
  call meshfem3D_create_addressing(myrank,NCHUNKS,NPROC,NPROC_ETA,NPROC_XI,NPROCTOT, &
                        addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice, &
                        OUTPUT_FILES)


  ! this for the different counters (which are now different if the superbrick is cut in the outer core)
  call meshfem3D_setup_counters(myrank, &
                        NSPEC1D_RADIAL,NSPEC2D_XI,NSPEC2D_ETA,NGLOB1D_RADIAL, &
                        DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA, &
                        CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                        NPROCTOT,iproc_xi_slice,iproc_eta_slice, &
                        NSPEC1D_RADIAL_CORNER,NSPEC2D_XI_FACE, &
                        NSPEC2D_ETA_FACE,NGLOB1D_RADIAL_CORNER)

  ! user output
  if(myrank == 0) call meshfem3D_output_info(myrank,sizeprocs,NEX_XI,NEX_ETA, &
                                           NPROC_XI,NPROC_ETA,NPROC,NCHUNKS,NPROCTOT, &
                                           R_CENTRAL_CUBE)

  ! distributes 3D models
  call meshfem3D_models_broadcast(myrank,NSPEC, &
                                MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,&
                                R80,R220,R670,RCMB,RICB)


  if(myrank == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) 'model setup successfully read in'
    write(IMAIN,*)
  endif

  ! get addressing for this process
  ichunk = ichunk_slice(myrank)
  iproc_xi = iproc_xi_slice(myrank)
  iproc_eta = iproc_eta_slice(myrank)

  ! volume of the slice
  volume_total = ZERO

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

    ! compute maximum number of points
    npointot = NSPEC(iregion_code) * NGLLX * NGLLY * NGLLZ

    ! use dynamic allocation to allocate memory for arrays
    allocate(idoubling(NSPEC(iregion_code)))
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)))
    allocate(xstore(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)))
    allocate(ystore(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)))
    allocate(zstore(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)))

! this for non blocking MPI
    allocate(is_on_a_slice_edge(NSPEC(iregion_code)))

    ! create all the regions of the mesh
    ! perform two passes in this part to be able to save memory
    do ipass = 1,2
      call create_regions_mesh(iregion_code,ibool,idoubling,is_on_a_slice_edge, &
                          xstore,ystore,zstore,rmins,rmaxs, &
                          iproc_xi,iproc_eta,ichunk,NSPEC(iregion_code),nspec_tiso, &
                          volume_local,area_local_bottom,area_local_top, &
                          nglob(iregion_code),npointot, &
                          NSTEP,DT, &
                          NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                          NSPEC2DMAX_XMIN_XMAX(iregion_code),NSPEC2DMAX_YMIN_YMAX(iregion_code), &
                          NSPEC2D_BOTTOM(iregion_code),NSPEC2D_TOP(iregion_code), &
                          NPROC_XI,NPROC_ETA,NSPEC2D_XI_FACE, &
                          NSPEC2D_ETA_FACE,NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER, &
                          myrank,LOCAL_PATH,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
                          SAVE_MESH_FILES,NCHUNKS,INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS, &
                          R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,&
                          RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
                          ner,ratio_sampling_array,doubling_index,r_bottom,r_top,&
                          this_region_has_a_doubling,ipass,ratio_divide_central_cube, &
                          CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                          mod(iproc_xi_slice(myrank),2),mod(iproc_eta_slice(myrank),2),USE_FULL_TISO_MANTLE, &
                          ATT1,ATT2,ATT3,SIMULATION_TYPE,USE_LDDRK,EXACT_MASS_MATRIX_FOR_ROTATION)
    enddo

    ! checks number of anisotropic elements found in the mantle
    if(iregion_code /= IREGION_CRUST_MANTLE .and. nspec_tiso /= 0 ) &
      call exit_MPI(myrank,'found anisotropic elements outside of the mantle')

    if( TRANSVERSE_ISOTROPY ) then
      if(iregion_code == IREGION_CRUST_MANTLE .and. nspec_tiso == 0) &
        call exit_MPI(myrank,'found no anisotropic elements in the mantle')
    endif

    ! computes total area and volume
    call meshfem3D_compute_area(myrank,NCHUNKS,iregion_code, &
                              area_local_bottom,area_local_top,&
                              volume_local,volume_total, &
                              RCMB,RICB,R_CENTRAL_CUBE,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES)

    ! create chunk buffers if more than one chunk
    if(NCHUNKS > 1) then
      call create_chunk_buffers(iregion_code,NSPEC(iregion_code),ibool,idoubling, &
                              xstore,ystore,zstore, &
                              nglob(iregion_code), &
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
    call MPI_BARRIER(MPI_COMM_WORLD,ier)

! end of loop on all the regions
  enddo

  if(myrank == 0) then
    ! check volume of chunk
    write(IMAIN,*)
    write(IMAIN,*) 'calculated volume: ',volume_total
    if((NCHUNKS == 6 .or. (abs(ANGULAR_WIDTH_XI_IN_DEGREES - 90.d0) < TINYVAL .and. &
                           abs(ANGULAR_WIDTH_ETA_IN_DEGREES - 90.d0) < TINYVAL)) .and. .not. TOPOGRAPHY) then
      ! take the central cube into account
      ! it is counted 6 times because of the fictitious elements
      if(INCLUDE_CENTRAL_CUBE) then
        write(IMAIN,*) '  similar volume: ', &
          dble(NCHUNKS)*((4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3)+5.*(2.*(R_CENTRAL_CUBE/R_EARTH)/sqrt(3.))**3)/6.d0
      else
        write(IMAIN,*) '  similar volume: ', &
          dble(NCHUNKS)*((4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3)-(2.*(R_CENTRAL_CUBE/R_EARTH)/sqrt(3.))**3)/6.d0
      endif
      write(IMAIN,*) '  (but not exact because the central cube is purposely inflated)'
    endif
  endif


!--- print number of points and elements in the mesh for each region

  if(myrank == 0) then

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

    ! evaluate the amount of static memory needed by the solver
    call memory_eval(OCEANS,ABSORBING_CONDITIONS,ATTENUATION,ANISOTROPIC_3D_MANTLE,&
                   TRANSVERSE_ISOTROPY,ANISOTROPIC_INNER_CORE,ROTATION,&
                   ONE_CRUST,doubling_index,this_region_has_a_doubling,&
                   ner,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_sampling_array,&
                   NSPEC,nglob,SIMULATION_TYPE,MOVIE_VOLUME,SAVE_FORWARD, &
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
                   NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION,ATT1,ATT2,ATT3, &
                   APPROXIMATE_HESS_KL,ANISOTROPIC_KL,NOISE_TOMOGRAPHY, &
                   NCHUNKS,USE_LDDRK,EXACT_MASS_MATRIX_FOR_ROTATION,static_memory_size)

    NGLOB1D_RADIAL_TEMP(:) = &
      (/maxval(NGLOB1D_RADIAL_CORNER(1,:)),maxval(NGLOB1D_RADIAL_CORNER(2,:)),maxval(NGLOB1D_RADIAL_CORNER(3,:))/)

    ! create include file for the solver
    call save_header_file(NSPEC,nglob,NEX_XI,NEX_ETA,NPROC,NPROCTOT, &
                    TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
                    ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS,ATTENUATION,ATTENUATION_3D, &
                    ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,NCHUNKS, &
                    INCLUDE_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES,&
                    CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,NSOURCES,NSTEP, &
                    static_memory_size,NGLOB1D_RADIAL_TEMP, &
                    NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NSPEC2D_TOP,NSPEC2D_BOTTOM, &
                    NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_XMIN_XMAX, &
                    NPROC_XI,NPROC_ETA, &
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
                    SIMULATION_TYPE,SAVE_FORWARD,MOVIE_VOLUME,SAVE_REGULAR_KL,NOISE_TOMOGRAPHY, &
                    ATT1,ATT2,ATT3,ATT4,ATT5,APPROXIMATE_HESS_KL,ANISOTROPIC_KL,PARTIAL_PHYS_DISPERSION_ONLY, &
                    SAVE_SOURCE_MASK,ABSORBING_CONDITIONS,USE_LDDRK,EXACT_MASS_MATRIX_FOR_ROTATION)

  endif   ! end of section executed by main process only

  ! deallocate arrays used for mesh generation
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
    ! close main output file
    close(IMAIN)
  endif

  ! synchronize all the processes to make sure everybody has finished
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

  end program xmeshfem3D

!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_create_addressing(myrank,NCHUNKS,NPROC,NPROC_ETA,NPROC_XI,NPROCTOT, &
                        addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice, &
                        OUTPUT_FILES)

  implicit none

  include "constants.h"

  integer :: myrank,NCHUNKS,NPROC,NPROC_ETA,NPROC_XI,NPROCTOT

  integer, dimension(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1) :: addressing
  integer, dimension(0:NPROCTOT-1) :: ichunk_slice,iproc_xi_slice,iproc_eta_slice

  character(len=150) OUTPUT_FILES

  ! local parameters
  integer ichunk,iproc_eta,iproc_xi,iprocnum,ier

  ! initializes
  addressing(:,:,:) = 0
  ichunk_slice(:) = 0
  iproc_xi_slice(:) = 0
  iproc_eta_slice(:) = 0

  ! loop on all the chunks to create global slice addressing for solver
  if(myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/addressing.txt',status='unknown',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening addressing.txt')
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

  end subroutine meshfem3D_create_addressing


!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_setup_counters(myrank, &
                        NSPEC1D_RADIAL,NSPEC2D_XI,NSPEC2D_ETA,NGLOB1D_RADIAL, &
                        DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA, &
                        CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                        NPROCTOT,iproc_xi_slice,iproc_eta_slice, &
                        NSPEC1D_RADIAL_CORNER,NSPEC2D_XI_FACE, &
                        NSPEC2D_ETA_FACE,NGLOB1D_RADIAL_CORNER)

! returns: NSPEC1D_RADIAL_CORNER,NSPEC2D_XI_FACE,
!              NSPEC2D_ETA_FACE,NGLOB1D_RADIAL_CORNER

  implicit none

  include "constants.h"

  integer myrank

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2D_XI,NSPEC2D_ETA, &
                                         NSPEC1D_RADIAL,NGLOB1D_RADIAL

  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  ! addressing for all the slices
  integer :: NPROCTOT
  integer, dimension(0:NPROCTOT-1) :: iproc_xi_slice,iproc_eta_slice

  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA

! this for the different corners of the slice (which are different if the superbrick is cut)
! 1 : xi_min, eta_min
! 2 : xi_max, eta_min
! 3 : xi_max, eta_max
! 4 : xi_min, eta_max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: &
    NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER
! 1 -> min, 2 -> max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE


  ! local parameters
  integer :: iregion

  do iregion=1,MAX_NUM_REGIONS
    NSPEC1D_RADIAL_CORNER(iregion,:) = NSPEC1D_RADIAL(iregion)
    NSPEC2D_XI_FACE(iregion,:) = NSPEC2D_XI(iregion)
    NSPEC2D_ETA_FACE(iregion,:) = NSPEC2D_ETA(iregion)
    NGLOB1D_RADIAL_CORNER(iregion,:) = NGLOB1D_RADIAL(iregion)
  enddo

  if (CUT_SUPERBRICK_XI) then
    if (CUT_SUPERBRICK_ETA) then
      if (mod(iproc_xi_slice(myrank),2) == 0) then
        if (mod(iproc_eta_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,1)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,1)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,1)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,2)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,2)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,2)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      else
        if (mod(iproc_eta_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,3)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,3)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,3)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,3)*(NGLLZ-1))
        else
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,4)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,4)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,4)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,4)*(NGLLZ-1))
        endif
      endif
    else
      if (mod(iproc_xi_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,1)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,1)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,1)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
      else
        NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,2)
        NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,2)
        NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,2)
        NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                      + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
      endif
    endif
  else
    if (CUT_SUPERBRICK_ETA) then
      if (mod(iproc_eta_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,1)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,1)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,1)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
      else
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,2)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,2)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,2)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
      endif
    endif
  endif

  end subroutine meshfem3D_setup_counters


!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_output_info(myrank,sizeprocs,NEX_XI,NEX_ETA, &
                                NPROC_XI,NPROC_ETA,NPROC,NCHUNKS,NPROCTOT,&
                                R_CENTRAL_CUBE)

  use meshfem3D_models_par

  implicit none

  integer :: myrank,sizeprocs,NEX_XI,NEX_ETA, &
           NPROC_XI,NPROC_ETA,NPROC,NCHUNKS,NPROCTOT
  double precision :: R_CENTRAL_CUBE

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
  if(ISOTROPIC_3D_MANTLE) then
    write(IMAIN,*) 'incorporating 3-D lateral variations'
  else
    write(IMAIN,*) 'no 3-D lateral variations'
  endif
  write(IMAIN,*)
  if(HETEROGEN_3D_MANTLE) then
    write(IMAIN,*) 'incorporating heterogeneities in the mantle'
  else
    write(IMAIN,*) 'no heterogeneities in the mantle'
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
    write(IMAIN,*) 'using unmodified 1D crustal model with two layers'
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
    if(ATTENUATION_3D) write(IMAIN,*)'using 3D attenuation'
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
  if(ANISOTROPIC_INNER_CORE) then
    write(IMAIN,*) 'incorporating anisotropic inner core'
  else
    write(IMAIN,*) 'no inner-core anisotropy'
  endif
  write(IMAIN,*)
  if(ANISOTROPIC_3D_MANTLE) then
    write(IMAIN,*) 'incorporating anisotropic mantle'
  else
    write(IMAIN,*) 'no general mantle anisotropy'
  endif
  write(IMAIN,*)
  write(IMAIN,*) 'Reference radius of the Earth used is ',R_EARTH_KM,' km'
  write(IMAIN,*)
  write(IMAIN,*) 'Central cube is at a radius of ',R_CENTRAL_CUBE/1000.d0,' km'

  end subroutine meshfem3D_output_info

!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_compute_area(myrank,NCHUNKS,iregion_code, &
                                    area_local_bottom,area_local_top,&
                                    volume_local,volume_total, &
                                    RCMB,RICB,R_CENTRAL_CUBE,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES)

  use meshfem3D_models_par
  use mpi

  implicit none

  integer :: myrank,NCHUNKS,iregion_code

  double precision :: area_local_bottom,area_local_top,volume_local,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES
  double precision :: volume_total
  double precision :: RCMB,RICB,R_CENTRAL_CUBE

  ! local parameters
  double precision :: volume_total_region,area_total_bottom,area_total_top
  integer :: ier

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

    !   check volume of chunk, and bottom and top area
    write(IMAIN,*)
    write(IMAIN,*) '   calculated top area: ',area_total_top

    ! compare to exact theoretical value
    if((NCHUNKS == 6 .or. (abs(ANGULAR_WIDTH_XI_IN_DEGREES - 90.d0) < TINYVAL .and. &
                           abs(ANGULAR_WIDTH_ETA_IN_DEGREES - 90.d0) < TINYVAL)) .and. .not. TOPOGRAPHY) then
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
    if((NCHUNKS == 6 .or. (abs(ANGULAR_WIDTH_XI_IN_DEGREES - 90.d0) < TINYVAL .and. &
                           abs(ANGULAR_WIDTH_ETA_IN_DEGREES - 90.d0) < TINYVAL)) .and. .not. TOPOGRAPHY) then
      select case(iregion_code)
        case(IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2
        case(IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2
        case(IREGION_INNER_CORE)
          write(IMAIN,*) '            similar (but not exact) area (central cube): ', &
                                                    dble(NCHUNKS)*(2.*(R_CENTRAL_CUBE / R_EARTH)/sqrt(3.))**2
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
    endif

  endif

  end subroutine meshfem3D_compute_area

