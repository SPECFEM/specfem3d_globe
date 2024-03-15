
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
 ! maximum number of points per region =       750417
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      2251251
 !
 ! total elements per slice =        12432
 ! total points per slice =       832031
 !
 ! the time step of the solver will be DT =   0.125000000      (s)
 ! the (approximate) minimum period resolved will be =    15.4437399      (s)
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !    49728.000000000000
 ! approximate total number of points in entire mesh =
 !    3328124.0000000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    9361516.0000000000
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
 ! spectral elements along a great circle =         1152
 ! GLL points along a great circle =         4608
 ! average distance between points in degrees =    7.81250000E-02
 ! average distance between points in km =    8.68710327
 ! average size of a spectral element in km =    34.7484131
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
 ! size of static arrays per slice =    533.77149199999997       MB
 !                                 =    509.04416275024414       MiB
 !                                 =   0.53377149199999996       GB
 !                                 =   0.49711344018578529       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    2135.0859679999999       MB
 !                                      =    2036.1766510009766       MiB
 !                                      =    2.1350859679999998       GB
 !                                      =    1.9884537607431412       GiB
 !                                      =    2.1350859680000001E-003  TB
 !                                      =    1.9418493757257238E-003  TiB
 !

 integer, parameter :: NEX_XI_VAL =           64
 integer, parameter :: NEX_ETA_VAL =           64

 integer, parameter :: NSPEC_CRUST_MANTLE =        11264
 integer, parameter :: NSPEC_OUTER_CORE =         1120
 integer, parameter :: NSPEC_INNER_CORE =           48

 integer, parameter :: NGLOB_CRUST_MANTLE =       750417
 integer, parameter :: NGLOB_OUTER_CORE =        77857
 integer, parameter :: NGLOB_INNER_CORE =         3757

 integer, parameter :: NSPECMAX_ANISO_IC =            0

 integer, parameter :: NSPECMAX_ISO_MANTLE =        11264
 integer, parameter :: NSPECMAX_TISO_MANTLE =        11264
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            0

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =        11264
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =           48

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =        11264
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =           48

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

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =        11264
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =         1120

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =       750417

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

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

 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =         1120

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .true.

 integer, parameter :: NPROC_XI_VAL =            2
 integer, parameter :: NPROC_ETA_VAL =            2
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =            4

 integer, parameter :: ATT1_VAL =            5
 integer, parameter :: ATT2_VAL =            5
 integer, parameter :: ATT3_VAL =            5
 integer, parameter :: ATT4_VAL =        11264
 integer, parameter :: ATT5_VAL =           48

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          640
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          640
 integer, parameter :: NSPEC2D_BOTTOM_CM =           64
 integer, parameter :: NSPEC2D_TOP_CM =         1024
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           12
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           12
 integer, parameter :: NSPEC2D_BOTTOM_IC =           16
 integer, parameter :: NSPEC2D_TOP_IC =           16
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          168
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          168
 integer, parameter :: NSPEC2D_BOTTOM_OC =           16
 integer, parameter :: NSPEC2D_TOP_OC =           64
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1

 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 0
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 0

 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 0
 integer, parameter :: NGLOB_XY_CM =       750417
 integer, parameter :: NGLOB_XY_IC =            0

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 logical, parameter :: UNDO_ATTENUATION_VAL = .false.
 integer, parameter :: NT_DUMP_ATTENUATION_VAL =          301

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    40.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =    25.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

