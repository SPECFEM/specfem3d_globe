
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
 ! maximum number of points per region =        72961
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =       218883
 !
 ! total elements per slice =         1056
 ! total points per slice =        72961
 !
 ! the time step of the solver will be DT =    5.00000007E-02  (s)
 ! the (approximate) minimum period resolved will be =    38.6093483      (s)
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !    4224.0000000000000
 ! approximate total number of points in entire mesh =
 !    291844.00000000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    875532.00000000000
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    30.0000000
 ! angular size in second direction in degrees =    25.0000000
 !
 ! longitude of center in degrees =   -150.000000
 ! latitude of center in degrees =    64.0000000
 !
 ! angle of rotation of the first chunk =    20.0000000
 !
 ! corner            1
 ! longitude in degrees =   -163.95780150968133
 ! latitude in degrees =    46.636809740502031
 !
 ! corner            2
 ! longitude in degrees =   -118.90706170568004
 ! latitude in degrees =    53.371981420299001
 !
 ! corner            3
 ! longitude in degrees =    165.51003142241825
 ! latitude in degrees =    63.948175763900551
 !
 ! corner            4
 ! longitude in degrees =   -103.80267130025246
 ! latitude in degrees =    76.771785045788704
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          384
 ! GLL points along a great circle =         1536
 ! average distance between points in degrees =   0.234375000
 ! average distance between points in km =    26.0613117
 ! average size of a spectral element in km =    104.245247
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
 ! size of static arrays per slice =    91.972099999999998       MB
 !                                 =    87.711429595947266       MiB
 !                                 =    9.1972100000000001E-002  GB
 !                                 =    8.5655692964792252E-002  GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    367.88839999999999       MB
 !                                      =    350.84571838378906       MiB
 !                                      =   0.36788840000000000       GB
 !                                      =   0.34262277185916901       GiB
 !                                      =    3.6788840000000000E-004  TB
 !                                      =    3.3459255064371973E-004  TiB
 !

 integer, parameter :: NEX_XI_VAL =           32
 integer, parameter :: NEX_ETA_VAL =           32

 integer, parameter :: NSPEC_CRUST_MANTLE =         1056
 integer, parameter :: NSPEC_OUTER_CORE =            0
 integer, parameter :: NSPEC_INNER_CORE =            0

 integer, parameter :: NGLOB_CRUST_MANTLE =        72961
 integer, parameter :: NGLOB_OUTER_CORE =            0
 integer, parameter :: NGLOB_INNER_CORE =            0

 integer, parameter :: NSPECMAX_ANISO_IC =            0

 integer, parameter :: NSPECMAX_ISO_MANTLE =         1056
 integer, parameter :: NSPECMAX_TISO_MANTLE =            0
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =         1056
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =         1056
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

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =         1056
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =            0

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =        72961

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .false.

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

 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =            0

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.

 integer, parameter :: NPROC_XI_VAL =            2
 integer, parameter :: NPROC_ETA_VAL =            2
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =            4

 integer, parameter :: ATT1_VAL =            5
 integer, parameter :: ATT2_VAL =            5
 integer, parameter :: ATT3_VAL =            5
 integer, parameter :: ATT4_VAL =         1056
 integer, parameter :: ATT5_VAL =            0

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =           96
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =           96
 integer, parameter :: NSPEC2D_BOTTOM_CM =           16
 integer, parameter :: NSPEC2D_TOP_CM =          256
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =            0
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =            0
 integer, parameter :: NSPEC2D_BOTTOM_IC =            0
 integer, parameter :: NSPEC2D_TOP_IC =            0
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =            0
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =            0
 integer, parameter :: NSPEC2D_BOTTOM_OC =            0
 integer, parameter :: NSPEC2D_TOP_OC =            0
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1

 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 0
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 0

 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 0
 integer, parameter :: NGLOB_XY_CM =        72961
 integer, parameter :: NGLOB_XY_IC =            0

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .false.

 logical, parameter :: UNDO_ATTENUATION_VAL = .false.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =         4043

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    25.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    30.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    64.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =  -150.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =    20.000000

