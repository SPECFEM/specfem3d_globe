
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
 ! maximum number of points per region =       215769
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =       647307
 !
 ! total elements per slice =         3280
 ! total points per slice =       223251
 !
 ! the time step of the solver will be DT =   0.189999998      (s)
 ! the (approximate) minimum period resolved will be =    6.93359327      (s)
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !    13120.000000000000
 ! approximate total number of points in entire mesh =
 !    893004.00000000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    2627580.0000000000
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    15.0000000
 ! angular size in second direction in degrees =    15.0000000
 !
 ! longitude of center in degrees =    136.000000
 ! latitude of center in degrees =    4.50000000
 !
 ! angle of rotation of the first chunk =    0.00000000
 !
 ! corner            1
 ! longitude in degrees =    128.55299363952108
 ! latitude in degrees =   -3.0628360932895293
 !
 ! corner            2
 ! longitude in degrees =    143.44700636047895
 ! latitude in degrees =   -3.0628360932895293
 !
 ! corner            3
 ! longitude in degrees =    128.40069001291837
 ! latitude in degrees =    11.982129369326749
 !
 ! corner            4
 ! longitude in degrees =    143.59930998708165
 ! latitude in degrees =    11.982129369326749
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          768
 ! GLL points along a great circle =         3072
 ! average distance between points in degrees =   0.117187500
 ! average distance between points in km =    6.93359327
 ! average size of a spectral element in km =    27.7343731
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
 ! size of static arrays per slice =    153.73640000000000       MB
 !                                 =    146.61445617675781       MiB
 !                                 =   0.15373640000000000       GB
 !                                 =   0.14317817986011505       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    614.94560000000001       MB
 !                                      =    586.45782470703125       MiB
 !                                      =   0.61494559999999998       GB
 !                                      =   0.57271271944046021       GiB
 !                                      =    6.1494559999999998E-004  TB
 !                                      =    5.5928976507857442E-004  TiB
 !

 integer, parameter :: NEX_XI_VAL =           32
 integer, parameter :: NEX_ETA_VAL =           32

 integer, parameter :: NSPEC_CRUST_MANTLE =         3184
 integer, parameter :: NSPEC_OUTER_CORE =           84
 integer, parameter :: NSPEC_INNER_CORE =           12

 integer, parameter :: NGLOB_CRUST_MANTLE =       215769
 integer, parameter :: NGLOB_OUTER_CORE =         6429
 integer, parameter :: NGLOB_INNER_CORE =         1053

 integer, parameter :: NSPECMAX_ANISO_IC =            0

 integer, parameter :: NSPECMAX_ISO_MANTLE =         3184
 integer, parameter :: NSPECMAX_TISO_MANTLE =            0
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =         3184
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =           12

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =         3184
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =           12

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

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =         3184
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =           84

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

 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =           84

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.

 integer, parameter :: NPROC_XI_VAL =            2
 integer, parameter :: NPROC_ETA_VAL =            2
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =            4

 integer, parameter :: ATT1_VAL =            5
 integer, parameter :: ATT2_VAL =            5
 integer, parameter :: ATT3_VAL =            5
 integer, parameter :: ATT4_VAL =         3184
 integer, parameter :: ATT5_VAL =           12

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          300
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          300
 integer, parameter :: NSPEC2D_BOTTOM_CM =           16
 integer, parameter :: NSPEC2D_TOP_CM =          256
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =            6
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =            6
 integer, parameter :: NSPEC2D_BOTTOM_IC =            4
 integer, parameter :: NSPEC2D_TOP_IC =            4
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =           26
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =           26
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
 integer, parameter :: NGLOB_XY_CM =       215769
 integer, parameter :: NGLOB_XY_IC =            0

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 logical, parameter :: UNDO_ATTENUATION_VAL = .false.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =         1297

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    15.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    15.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =     4.500000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =   136.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

