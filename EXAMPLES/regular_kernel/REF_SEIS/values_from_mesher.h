
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 !
 ! number of chunks =            6
 !
 ! these statistics include the central cube
 !
 ! number of processors =            6
 !
 ! maximum number of points per region =       576013
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      1728039
 !
 ! total elements per slice =         9424
 ! total points per slice =       631223
 !
 ! the time step of the solver will be DT =   0.118749999      (s)
 ! the (approximate) minimum period resolved will be =    136.000000      (s)
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !    56224.000000000000
 ! approximate total number of points in entire mesh =
 !    3762773.0000000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    10712499.000000000
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          128
 ! GLL points along a great circle =          512
 ! average distance between points in degrees =   0.703125000
 ! average distance between points in km =    78.1839294
 ! average size of a spectral element in km =    312.735718
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
 ! size of static arrays per slice =    419.10752000000002       MB
 !                                 =    399.69207763671875       MiB
 !                                 =   0.41910752000000001       GB
 !                                 =   0.39032429456710815       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    2514.6451200000001       MB
 !                                      =    2398.1524658203125       MiB
 !                                      =    2.5146451200000000       GB
 !                                      =    2.3419457674026489       GiB
 !                                      =    2.5146451200000002E-003  TB
 !                                      =    2.2870564134791493E-003  TiB
 !

 integer, parameter :: NEX_XI_VAL =           32
 integer, parameter :: NEX_ETA_VAL =           32

 integer, parameter :: NSPEC_CRUST_MANTLE =         8640
 integer, parameter :: NSPEC_OUTER_CORE =          688
 integer, parameter :: NSPEC_INNER_CORE =           96

 integer, parameter :: NGLOB_CRUST_MANTLE =       576013
 integer, parameter :: NGLOB_OUTER_CORE =        47985
 integer, parameter :: NGLOB_INNER_CORE =         7225

 integer, parameter :: NSPECMAX_ANISO_IC =            1

 integer, parameter :: NSPECMAX_ISO_MANTLE =         8640
 integer, parameter :: NSPECMAX_TISO_MANTLE =         8640
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =            1
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =         8640
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =           96

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT =            1
 integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY =         8640
 integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY =           96

 integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT =         8640
 integer, parameter :: NSPEC_OUTER_CORE_ADJOINT =          688
 integer, parameter :: NSPEC_INNER_CORE_ADJOINT =           96
 integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT =       576013
 integer, parameter :: NGLOB_OUTER_CORE_ADJOINT =        47985
 integer, parameter :: NGLOB_INNER_CORE_ADJOINT =         7225
 integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT =          688

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =            1
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =            1

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =            1

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .false.

 logical, parameter :: ATTENUATION_3D_VAL = .false.

 logical, parameter :: ELLIPTICITY_VAL = .false.

 logical, parameter :: GRAVITY_VAL = .true.

 logical, parameter :: OCEANS_VAL = .false.

 integer, parameter :: NX_BATHY_VAL = 1
 integer, parameter :: NY_BATHY_VAL = 1

 logical, parameter :: ROTATION_VAL = .true.
 logical, parameter :: EXACT_MASS_MATRIX_FOR_ROTATION_VAL = .true.

 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =          688

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.

 integer, parameter :: NPROC_XI_VAL =            1
 integer, parameter :: NPROC_ETA_VAL =            1
 integer, parameter :: NCHUNKS_VAL =            6
 integer, parameter :: NPROCTOT_VAL =            6

 integer, parameter :: ATT1_VAL =            1
 integer, parameter :: ATT2_VAL =            1
 integer, parameter :: ATT3_VAL =            1
 integer, parameter :: ATT4_VAL =            1
 integer, parameter :: ATT5_VAL =            1

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          440
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          440
 integer, parameter :: NSPEC2D_BOTTOM_CM =           64
 integer, parameter :: NSPEC2D_TOP_CM =         1024
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           24
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           24
 integer, parameter :: NSPEC2D_BOTTOM_IC =           16
 integer, parameter :: NSPEC2D_TOP_IC =           16
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          100
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          100
 integer, parameter :: NSPEC2D_BOTTOM_OC =           16
 integer, parameter :: NSPEC2D_TOP_OC =           64
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1

 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 1
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 1

 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 1
 integer, parameter :: NGLOB_XY_CM =       576013
 integer, parameter :: NGLOB_XY_IC =         7225

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 logical, parameter :: UNDO_ATTENUATION_VAL = .false.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =          403

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    90.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    90.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =     0.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =     0.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

