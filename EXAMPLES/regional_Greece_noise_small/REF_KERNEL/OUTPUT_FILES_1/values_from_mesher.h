
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
 ! maximum number of points per region =       358433
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      1075299
 !
 ! total elements per slice =         5841
 ! total points per slice =       395363
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
 !    1581452.0000000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    4466492.0000000000
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
 ! longitude in degrees =    13.682753768668109
 ! latitude in degrees =    29.492489378095556
 !
 ! corner            2
 ! longitude in degrees =    36.317246231331893
 ! latitude in degrees =    29.492489378095556
 !
 ! corner            3
 ! longitude in degrees =    9.9387399017687663
 ! latitude in degrees =    49.011233731524619
 !
 ! corner            4
 ! longitude in degrees =    40.061260098231237
 ! latitude in degrees =    49.011233731524619
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
 ! size of static arrays per slice =    306.35236400000002       MB
 !                                 =    292.16038131713867       MiB
 !                                 =   0.30635236399999999       GB
 !                                 =   0.28531287238001823       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    1225.4094560000001       MB
 !                                      =    1168.6415252685547       MiB
 !                                      =    1.2254094560000000       GB
 !                                      =    1.1412514895200729       GiB
 !                                      =    1.2254094559999999E-003  TB
 !                                      =    1.1145034077344462E-003  TiB
 !

 integer, parameter :: NEX_XI_VAL =           48
 integer, parameter :: NEX_ETA_VAL =           48

 integer, parameter :: NSPEC_CRUST_MANTLE =         5328
 integer, parameter :: NSPEC_OUTER_CORE =          486
 integer, parameter :: NSPEC_INNER_CORE =           27
 integer, parameter :: NSPEC_TRINFINITE =            0
 integer, parameter :: NSPEC_INFINITE =            0

 integer, parameter :: NGLOB_CRUST_MANTLE =       358433
 integer, parameter :: NGLOB_OUTER_CORE =        34733
 integer, parameter :: NGLOB_INNER_CORE =         2197
 integer, parameter :: NGLOB_TRINFINITE =            0
 integer, parameter :: NGLOB_INFINITE =            0

 integer, parameter :: NSPECMAX_ANISO_IC =            0

 integer, parameter :: NSPECMAX_ISO_MANTLE =         5328
 integer, parameter :: NSPECMAX_TISO_MANTLE =            0
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =         5328
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =           27

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =         5328
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =           27

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT =         5328
 integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT =           27

 integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY =         5328
 integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY =           27

 integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT =         5328
 integer, parameter :: NSPEC_OUTER_CORE_ADJOINT =          486
 integer, parameter :: NSPEC_INNER_CORE_ADJOINT =           27
 integer, parameter :: NSPEC_TRINFINITE_ADJOINT =            0
 integer, parameter :: NSPEC_INFINITE_ADJOINT =            0

 integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT =       358433
 integer, parameter :: NGLOB_OUTER_CORE_ADJOINT =        34733
 integer, parameter :: NGLOB_INNER_CORE_ADJOINT =         2197
 integer, parameter :: NGLOB_TRINFINITE_ADJOINT =            0
 integer, parameter :: NGLOB_INFINITE_ADJOINT =            0

 integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT =          486

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =         5328
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =          486

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =       358433

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .false.

 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .true.

 logical, parameter :: ATTENUATION_3D_VAL = .false.

 logical, parameter :: ELLIPTICITY_VAL = .true.

 logical, parameter :: GRAVITY_VAL = .true.

 logical, parameter :: FULL_GRAVITY_VAL = .false.

 logical, parameter :: OCEANS_VAL = .true.

 integer, parameter :: NX_BATHY_VAL =         5400
 integer, parameter :: NY_BATHY_VAL =         2700

 logical, parameter :: ROTATION_VAL = .true.
 logical, parameter :: EXACT_MASS_MATRIX_FOR_ROTATION_VAL = .false.

 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =          486

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .true.

 integer, parameter :: NPROC_XI_VAL =            2
 integer, parameter :: NPROC_ETA_VAL =            2
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =            4

 integer, parameter :: ATT1_VAL =            5
 integer, parameter :: ATT2_VAL =            5
 integer, parameter :: ATT3_VAL =            5
 integer, parameter :: ATT4_VAL =         5328
 integer, parameter :: ATT5_VAL =           27

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          384
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          384
 integer, parameter :: NSPEC2D_BOTTOM_CM =           36
 integer, parameter :: NSPEC2D_TOP_CM =          576

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =            9
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =            9
 integer, parameter :: NSPEC2D_BOTTOM_IC =            9
 integer, parameter :: NSPEC2D_TOP_IC =            9

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =           99
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          102
 integer, parameter :: NSPEC2D_BOTTOM_OC =            9
 integer, parameter :: NSPEC2D_TOP_OC =           36

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_TRINF =            0
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_TRINF =            0
 integer, parameter :: NSPEC2D_BOTTOM_TRINF =            0
 integer, parameter :: NSPEC2D_TOP_TRINF =            0

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_INF =            0
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_INF =            0
 integer, parameter :: NSPEC2D_BOTTOM_INF =            0
 integer, parameter :: NSPEC2D_TOP_INF =            0

 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1

 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 0
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 0

 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 0
 integer, parameter :: NGLOB_XY_CM =       358433
 integer, parameter :: NGLOB_XY_IC =            0

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 logical, parameter :: UNDO_ATTENUATION_VAL = .false.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =          723

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    40.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =    25.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

