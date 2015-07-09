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
!
! United States and French Government Sponsorship Acknowledged.

  program xmeshfem3D

  implicit none

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
! number = 1,
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
! v. 7.0, many developers, January 2015:
!     simultaneous MPI runs, ADIOS file I/O support, ASDF seismograms, new seismogram names, tomography tools,
!     CUDA and OpenCL GPU support, CEM model support, updates AK135 model, binary topography files,
!     fixes geocentric/geographic conversions, updates ellipticity and gravity factors, git versioning system.
!
! v. 6.0, Daniel Peter (ETH Z\"urich, Switzerland), Dimitri Komatitsch and Zhinan Xie (CNRS / University of Marseille, France),
!     Elliott Sales de Andrade (University of Toronto, Canada), and many others, in particular from Princeton University, USA,
!     April 2014:
!     more flexible MPI implementation, GPU support, exact undoing of attenuation, LDDRK4-6 higher-order time scheme, etc...
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
!      first port to GPUs using CUDA, new doubling brick in the mesh, new perfectly load-balanced mesh,
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
!                     can be done in routine get_timestep_and_layers()
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
  call init_mpi()

  ! initializes parameters
  call initialize_mesher()

  ! setup addressing and models
  call setup_model()

  ! creates meshes for regions crust/mantle, outer core and inner core
  call create_meshes()

  ! outputs mesh info and saves new header file
  call finalize_mesher()

  ! stop all the MPI processes, and exit
  call finalize_mpi()

  end program xmeshfem3D

