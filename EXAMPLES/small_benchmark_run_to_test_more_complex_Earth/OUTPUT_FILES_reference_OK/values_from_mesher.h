
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 !
 ! number of chunks =  6
 !
 ! these statistics include the central cube
 !
 ! number of processors =  24
 !
 ! maximum number of points per region =  921753
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =  2765259
 !
 ! total elements per slice =  15275
 ! total points per slice =  1017411
 !
 ! the time step of the solver will be DT =  0.118749999
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !  361600.
 ! approximate total number of points in entire mesh =
 !  24073259.
 ! approximate total number of degrees of freedom in entire mesh =
 !  68665425.
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =  320
 ! GLL points along a great circle =  1280
 ! average distance between points in degrees =  0.28125
 ! average distance between points in km =  31.2735729
 ! average size of a spectral element in km =  125.094292
 !

 ! approximate static memory needed by the solver:
 ! ----------------------------------------------
 !
 ! (lower bound, usually the real amount used is 5% to 10% higher)
 !
 ! (you can get a more precise estimate of the size used per MPI process
 !  by typing "size -d bin/xspecfem3D"
 !  after compiling the code with the DATA/Par_file you plan to use)
 !
 ! size of static arrays per slice =  418.43724799999995  MB
 !                                 =  399.0528564453125  MiB
 !                                 =  0.41843724800000004  GB
 !                                 =  0.38970005512237549  GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =  10.042493952000001  GB
 !                                      =  9.3528013229370117  GiB
 !                                      =  1.00424939519999992E-2  TB
 !                                      =  9.13359504193067551E-3  TiB
 !

 integer, parameter :: NEX_XI_VAL =  80
 integer, parameter :: NEX_ETA_VAL =  80

 integer, parameter :: NSPEC_CRUST_MANTLE =  13900
 integer, parameter :: NSPEC_OUTER_CORE =  1075
 integer, parameter :: NSPEC_INNER_CORE =  300

 integer, parameter :: NGLOB_CRUST_MANTLE =  921753
 integer, parameter :: NGLOB_OUTER_CORE =  74049
 integer, parameter :: NGLOB_INNER_CORE =  21609

 integer, parameter :: NSPECMAX_ANISO_IC =  1

 integer, parameter :: NSPECMAX_ISO_MANTLE =  13900
 integer, parameter :: NSPECMAX_TISO_MANTLE =  13900
 integer, parameter :: NSPECMAX_ANISO_MANTLE =  1

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =  13900
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =  300

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =  13900
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =  300

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT =  1
 integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT =  1

 integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY =  1
 integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY =  1

 integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT =  1
 integer, parameter :: NSPEC_OUTER_CORE_ADJOINT =  1
 integer, parameter :: NSPEC_INNER_CORE_ADJOINT =  1
 integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT =  1
 integer, parameter :: NGLOB_OUTER_CORE_ADJOINT =  1
 integer, parameter :: NGLOB_INNER_CORE_ADJOINT =  1
 integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT =  1

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =  1
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =  1

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =  921753

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .true.

 logical, parameter :: ATTENUATION_3D_VAL = .false.

 logical, parameter :: ELLIPTICITY_VAL = .true.

 logical, parameter :: GRAVITY_VAL = .true.

 logical, parameter :: OCEANS_VAL = .true.

 integer, parameter :: NX_BATHY_VAL =  5400
 integer, parameter :: NY_BATHY_VAL =  2700

 logical, parameter :: ROTATION_VAL = .true.
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =  1075

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.

 integer, parameter :: NPROC_XI_VAL =  2
 integer, parameter :: NPROC_ETA_VAL =  2
 integer, parameter :: NCHUNKS_VAL =  6
 integer, parameter :: NPROCTOT_VAL =  24

 integer, parameter :: ATT1_VAL =  5
 integer, parameter :: ATT2_VAL =  5
 integer, parameter :: ATT3_VAL =  5
 integer, parameter :: ATT4_VAL =  13900
 integer, parameter :: ATT5_VAL =  300

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =  570
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =  570
 integer, parameter :: NSPEC2D_BOTTOM_CM =  100
 integer, parameter :: NSPEC2D_TOP_CM =  1600
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =  60
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =  60
 integer, parameter :: NSPEC2D_BOTTOM_IC =  25
 integer, parameter :: NSPEC2D_TOP_IC =  25
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =  130
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =  135
 integer, parameter :: NSPEC2D_BOTTOM_OC =  25
 integer, parameter :: NSPEC2D_TOP_OC =  100
 integer, parameter :: NSPEC2D_MOHO =  1
 integer, parameter :: NSPEC2D_400 =  1
 integer, parameter :: NSPEC2D_670 =  1
 integer, parameter :: NSPEC2D_CMB =  1
 integer, parameter :: NSPEC2D_ICB =  1

 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 1
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 1

 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 1
 integer, parameter :: NM_KL_REG_PTS_VAL = 1

 integer, parameter :: NGLOB_XY_CM =  1
 integer, parameter :: NGLOB_XY_IC =  1

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 integer, parameter :: NT_DUMP_ATTENUATION =  100000000

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    90.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    90.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =     0.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =     0.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

