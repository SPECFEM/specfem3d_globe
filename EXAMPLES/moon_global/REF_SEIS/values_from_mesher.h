
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
 ! maximum number of points per region =      1068221
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      3204663
 !
 ! total elements per slice =        16884
 ! total points per slice =      1121671
 !
 ! the time step of the solver will be DT =   0.150000006      (s)
 ! the (approximate) minimum period resolved will be =    31.5813694      (s)
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !    100224.00000000000
 ! approximate total number of points in entire mesh =
 !    6651901.0000000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    19561803.000000000
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          192
 ! GLL points along a great circle =          768
 ! average distance between points in degrees =   0.468750000
 ! average distance between points in km =    14.2116165
 ! average size of a spectral element in km =    56.8464661
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
 ! size of static arrays per slice =    474.24883999999997       MB
 !                                 =    452.27893829345703       MiB
 !                                 =   0.47424884000000000       GB
 !                                 =   0.44167865067720413       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    2845.4930399999998       MB
 !                                      =    2713.6736297607422       MiB
 !                                      =    2.8454930400000000       GB
 !                                      =    2.6500719040632248       GiB
 !                                      =    2.8454930400000000E-003  TB
 !                                      =    2.5879608438117430E-003  TiB
 !

 integer, parameter :: NEX_XI_VAL =           48
 integer, parameter :: NEX_ETA_VAL =           48

 integer, parameter :: NSPEC_CRUST_MANTLE =        16128
 integer, parameter :: NSPEC_OUTER_CORE =          468
 integer, parameter :: NSPEC_INNER_CORE =          288

 integer, parameter :: NGLOB_CRUST_MANTLE =      1068221
 integer, parameter :: NGLOB_OUTER_CORE =        32825
 integer, parameter :: NGLOB_INNER_CORE =        20625

 integer, parameter :: NSPECMAX_ANISO_IC =            0

 integer, parameter :: NSPECMAX_ISO_MANTLE =        16128
 integer, parameter :: NSPECMAX_TISO_MANTLE =            0
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =        16128
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =          288

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =        16128
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =          288

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT =            0
 integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY =            0
 integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT =            0
 integer, parameter :: NSPEC_OUTER_CORE_ADJOINT =            0
 integer, parameter :: NSPEC_INNER_CORE_ADJOINT =            0
 integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT =            0
 integer, parameter :: NGLOB_OUTER_CORE_ADJOINT =            0
 integer, parameter :: NGLOB_INNER_CORE_ADJOINT =            0
 integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =            0
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =            0

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =            0

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .false.

 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .true.

 logical, parameter :: ATTENUATION_3D_VAL = .false.

 logical, parameter :: ELLIPTICITY_VAL = .true.

 logical, parameter :: GRAVITY_VAL = .true.

 logical, parameter :: FULL_GRAVITY_VAL = .false.

 logical, parameter :: OCEANS_VAL = .false.

 integer, parameter :: NX_BATHY_VAL =         5400
 integer, parameter :: NY_BATHY_VAL =         2700

 logical, parameter :: ROTATION_VAL = .true.
 logical, parameter :: EXACT_MASS_MATRIX_FOR_ROTATION_VAL = .false.

 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =          468

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.

 integer, parameter :: NPROC_XI_VAL =            1
 integer, parameter :: NPROC_ETA_VAL =            1
 integer, parameter :: NCHUNKS_VAL =            6
 integer, parameter :: NPROCTOT_VAL =            6

 integer, parameter :: ATT1_VAL =            5
 integer, parameter :: ATT2_VAL =            5
 integer, parameter :: ATT3_VAL =            5
 integer, parameter :: ATT4_VAL =        16128
 integer, parameter :: ATT5_VAL =          288

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          504
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          504
 integer, parameter :: NSPEC2D_BOTTOM_CM =          144
 integer, parameter :: NSPEC2D_TOP_CM =         2304
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           48
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           48
 integer, parameter :: NSPEC2D_BOTTOM_IC =           36
 integer, parameter :: NSPEC2D_TOP_IC =           36
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =           42
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =           42
 integer, parameter :: NSPEC2D_BOTTOM_OC =           36
 integer, parameter :: NSPEC2D_TOP_OC =          144
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1

 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 0
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 0

 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 0
 integer, parameter :: NGLOB_XY_CM =            0
 integer, parameter :: NGLOB_XY_IC =            0

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 logical, parameter :: UNDO_ATTENUATION_VAL = .false.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =          219

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    90.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    90.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =     0.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =     0.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

