
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
 ! number of processors =            1
 !
 ! maximum number of points per region =      1408505
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      4225515
 !
 ! total elements per slice =        23364
 ! total points per slice =      1548787
 !
 ! the time step of the solver will be DT =   0.165000007      (s)
 ! the (approximate) minimum period resolved will be =    20.5916538      (s)
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !    23364.000000000000
 ! approximate total number of points in entire mesh =
 !    1548787.0000000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    4382047.0000000000
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    20.0000000
 ! angular size in second direction in degrees =    20.0000000
 !
 ! longitude of center in degrees =    25.0000000
 ! latitude of center in degrees =    40.0000000
 !
 ! angle of rotation of the first chunk =    0.00000000
 !
 ! corner            1
 ! longitude in degrees =    13.661882747982684
 ! latitude in degrees =    29.677869656290330
 !
 ! corner            2
 ! longitude in degrees =    36.338117252017327
 ! latitude in degrees =    29.677869656290330
 !
 ! corner            3
 ! longitude in degrees =    9.8824891302835294
 ! latitude in degrees =    49.192981807843722
 !
 ! corner            4
 ! longitude in degrees =    40.117510869716469
 ! latitude in degrees =    49.192981807843722
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          864
 ! GLL points along a great circle =         3456
 ! average distance between points in degrees =   0.104166664
 ! average distance between points in km =    11.5828047
 ! average size of a spectral element in km =    46.3312187
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
 ! size of static arrays per slice =    952.55532400000004       MB
 !                                 =    908.42754745483398       MiB
 !                                 =   0.95255532399999998       GB
 !                                 =   0.88713627681136131       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    952.55532400000004       MB
 !                                      =    908.42754745483398       MiB
 !                                      =   0.95255532399999998       GB
 !                                      =   0.88713627681136131       GiB
 !                                      =    9.5255532400000001E-004  TB
 !                                      =    8.6634402032359503E-004  TiB
 !

 integer, parameter :: NEX_XI_VAL =           48
 integer, parameter :: NEX_ETA_VAL =           48

 integer, parameter :: NSPEC_CRUST_MANTLE =        21312
 integer, parameter :: NSPEC_OUTER_CORE =         1944
 integer, parameter :: NSPEC_INNER_CORE =          108

 integer, parameter :: NGLOB_CRUST_MANTLE =      1408505
 integer, parameter :: NGLOB_OUTER_CORE =       132157
 integer, parameter :: NGLOB_INNER_CORE =         8125

 integer, parameter :: NSPECMAX_ANISO_IC =            1

 integer, parameter :: NSPECMAX_ISO_MANTLE =        21312
 integer, parameter :: NSPECMAX_TISO_MANTLE =        21312
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =        21312
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =          108

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =        21312
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =          108

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT =            1
 integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY =            1
 integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT =            1
 integer, parameter :: NSPEC_OUTER_CORE_ADJOINT =            1
 integer, parameter :: NSPEC_INNER_CORE_ADJOINT =            1
 integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT =            1
 integer, parameter :: NGLOB_OUTER_CORE_ADJOINT =            1
 integer, parameter :: NGLOB_INNER_CORE_ADJOINT =            1
 integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =        21312
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =         1944

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =      1408505

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .true.

 logical, parameter :: ATTENUATION_3D_VAL = .false.

 logical, parameter :: ELLIPTICITY_VAL = .true.

 logical, parameter :: GRAVITY_VAL = .true.

 logical, parameter :: OCEANS_VAL = .true.

 integer, parameter :: NX_BATHY_VAL =         5400
 integer, parameter :: NY_BATHY_VAL =         2700

 logical, parameter :: ROTATION_VAL = .true.
 logical, parameter :: EXACT_MASS_MATRIX_FOR_ROTATION_VAL = .false.

 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =         1944

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.

 integer, parameter :: NPROC_XI_VAL =            1
 integer, parameter :: NPROC_ETA_VAL =            1
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =            1

 integer, parameter :: ATT1_VAL =            5
 integer, parameter :: ATT2_VAL =            5
 integer, parameter :: ATT3_VAL =            5
 integer, parameter :: ATT4_VAL =        21312
 integer, parameter :: ATT5_VAL =          108

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          768
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          768
 integer, parameter :: NSPEC2D_BOTTOM_CM =          144
 integer, parameter :: NSPEC2D_TOP_CM =         2304
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           18
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           18
 integer, parameter :: NSPEC2D_BOTTOM_IC =           36
 integer, parameter :: NSPEC2D_TOP_IC =           36
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          192
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          192
 integer, parameter :: NSPEC2D_BOTTOM_OC =           36
 integer, parameter :: NSPEC2D_TOP_OC =          144
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1

 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 1
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 1

 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 1
 integer, parameter :: NGLOB_XY_CM =      1408505
 integer, parameter :: NGLOB_XY_IC =            1

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 logical, parameter :: UNDO_ATTENUATION_VAL = .true.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =          114

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    40.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =    25.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

