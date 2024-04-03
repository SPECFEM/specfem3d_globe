
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
 ! maximum number of points per region =       153157
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =       459471
 !
 ! total elements per slice =         2420
 ! total points per slice =       166743
 !
 ! the time step of the solver will be DT =   0.284999996      (s)
 ! the (approximate) minimum period resolved will be =    30.8874798      (s)
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !    9680.0000000000000
 ! approximate total number of points in entire mesh =
 !    666972.00000000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    1898060.0000000000
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
 ! longitude in degrees =    13.661882747982682
 ! latitude in degrees =    29.513451805567765
 !
 ! corner            2
 ! longitude in degrees =    36.338117252017319
 ! latitude in degrees =    29.513451805567765
 !
 ! corner            3
 ! longitude in degrees =    9.8824891302835312
 ! latitude in degrees =    49.003504872368651
 !
 ! corner            4
 ! longitude in degrees =    40.117510869716469
 ! latitude in degrees =    49.003504872368651
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          576
 ! GLL points along a great circle =         2304
 ! average distance between points in degrees =   0.156250000
 ! average distance between points in km =    17.3742065
 ! average size of a spectral element in km =    69.4968262
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
 ! size of static arrays per slice =    59.599336000000001       MB
 !                                 =    56.838356018066406       MiB
 !                                 =    5.9599336000000003E-002  GB
 !                                 =    5.5506207048892975E-002  GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    238.39734400000000       MB
 !                                      =    227.35342407226562       MiB
 !                                      =   0.23839734400000001       GB
 !                                      =   0.22202482819557190       GiB
 !                                      =    2.3839734400000001E-004  TB
 !                                      =    2.1682112128473818E-004  TiB
 !

 integer, parameter :: NEX_XI_VAL =           32
 integer, parameter :: NEX_ETA_VAL =           32

 integer, parameter :: NSPEC_CRUST_MANTLE =         2240
 integer, parameter :: NSPEC_OUTER_CORE =          172
 integer, parameter :: NSPEC_INNER_CORE =            8

 integer, parameter :: NGLOB_CRUST_MANTLE =       153157
 integer, parameter :: NGLOB_OUTER_CORE =        12857
 integer, parameter :: NGLOB_INNER_CORE =          729

 integer, parameter :: NSPECMAX_ANISO_IC =            0

 integer, parameter :: NSPECMAX_ISO_MANTLE =         2240
 integer, parameter :: NSPECMAX_TISO_MANTLE =         2240
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

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =         2240
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

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          232
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          232
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
 integer, parameter :: NGLOB_XY_CM =       153157
 integer, parameter :: NGLOB_XY_IC =            0

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 logical, parameter :: UNDO_ATTENUATION_VAL = .true.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =         1806

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    40.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =    25.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

