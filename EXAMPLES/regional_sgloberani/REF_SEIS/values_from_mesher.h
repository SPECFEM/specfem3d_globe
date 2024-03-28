
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 !
 ! number of chunks =            1
 !
 ! these statistics do not include the central cube
 !
 ! number of processors =            4
 !
 ! maximum number of points per region =       147645
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =       442935
 !
 ! total elements per slice =         2340
 ! total points per slice =       161231
 !
 ! the time step of the solver will be DT =   0.115000002      (s)
 ! the (approximate) minimum period resolved will be =    138.993652      (s)
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !    9360.0000000000000
 ! approximate total number of points in entire mesh =
 !    644924.00000000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    1831916.0000000000
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    90.0000000
 ! angular size in second direction in degrees =    90.0000000
 !
 ! longitude of center in degrees =    140.000000
 ! latitude of center in degrees =    45.0000000
 !
 ! angle of rotation of the first chunk =   -130.000000
 !
 ! corner            1
 ! longitude in degrees =    343.09168009901288
 ! latitude in degrees =    79.546211208738200
 !
 ! corner            2
 ! longitude in degrees =    79.413162387413678
 ! latitude in degrees =    20.973056977842038
 !
 ! corner            3
 ! longitude in degrees =    206.24829303692437
 ! latitude in degrees =    27.294721195949393
 !
 ! corner            4
 ! longitude in degrees =    144.13890799393090
 ! latitude in degrees =   -9.6079169325811549
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
 ! size of static arrays per slice =    57.531911999999998       MB
 !                                 =    54.866706848144531       MiB
 !                                 =    5.7531911999999998E-002  GB
 !                                 =    5.3580768406391144E-002  GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    230.12764799999999       MB
 !                                      =    219.46682739257812       MiB
 !                                      =   0.23012764799999999       GB
 !                                      =   0.21432307362556458       GiB
 !                                      =    2.3012764799999999E-004  TB
 !                                      =    2.0929987658746541E-004  TiB
 !

 integer, parameter :: NEX_XI_VAL =           32
 integer, parameter :: NEX_ETA_VAL =           32

 integer, parameter :: NSPEC_CRUST_MANTLE =         2160
 integer, parameter :: NSPEC_OUTER_CORE =          172
 integer, parameter :: NSPEC_INNER_CORE =            8

 integer, parameter :: NGLOB_CRUST_MANTLE =       147645
 integer, parameter :: NGLOB_OUTER_CORE =        12857
 integer, parameter :: NGLOB_INNER_CORE =          729

 integer, parameter :: NSPECMAX_ANISO_IC =            0

 integer, parameter :: NSPECMAX_ISO_MANTLE =         2160
 integer, parameter :: NSPECMAX_TISO_MANTLE =         2160
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =            0
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =            0
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =            0

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

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =         2160
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =          172

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =            0

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .false.

 logical, parameter :: ATTENUATION_3D_VAL = .false.

 logical, parameter :: ELLIPTICITY_VAL = .false.

 logical, parameter :: GRAVITY_VAL = .false.

 logical, parameter :: FULL_GRAVITY_VAL = .false.

 logical, parameter :: OCEANS_VAL = .false.

 integer, parameter :: NX_BATHY_VAL = 0
 integer, parameter :: NY_BATHY_VAL = 0

 logical, parameter :: ROTATION_VAL = .false.
 logical, parameter :: EXACT_MASS_MATRIX_FOR_ROTATION_VAL = .false.

 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =            0

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.

 integer, parameter :: NPROC_XI_VAL =            2
 integer, parameter :: NPROC_ETA_VAL =            2
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =            4

 integer, parameter :: ATT1_VAL =            1
 integer, parameter :: ATT2_VAL =            1
 integer, parameter :: ATT3_VAL =            1
 integer, parameter :: ATT4_VAL =            1
 integer, parameter :: ATT5_VAL =            1

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          220
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          220
 integer, parameter :: NSPEC2D_BOTTOM_CM =           16
 integer, parameter :: NSPEC2D_TOP_CM =          256
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =            4
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =            4
 integer, parameter :: NSPEC2D_BOTTOM_IC =            4
 integer, parameter :: NSPEC2D_TOP_IC =            4
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =           50
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =           50
 integer, parameter :: NSPEC2D_BOTTOM_OC =            4
 integer, parameter :: NSPEC2D_TOP_OC =           16
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1

 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 0
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 0

 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 0
 integer, parameter :: NGLOB_XY_CM =       147645
 integer, parameter :: NGLOB_XY_IC =            0

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 logical, parameter :: UNDO_ATTENUATION_VAL = .false.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =         1884

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    90.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    90.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    45.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =   140.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =  -130.000000

