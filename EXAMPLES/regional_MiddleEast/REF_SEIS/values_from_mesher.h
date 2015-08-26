
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 !
 ! number of chunks =  1
 !
 ! these statistics do not include the central cube
 !
 ! number of processors =  64
 !
 ! maximum number of points per region =  192897
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =  578691
 !
 ! total elements per slice =  3108
 ! total points per slice =  214855
 !
 ! the time step of the solver will be DT =  0.126843646
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !  198912.
 ! approximate total number of points in entire mesh =
 !  13750720.
 ! approximate total number of degrees of freedom in entire mesh =
 !  38576320.
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =  45.
 ! angular size in second direction in degrees =  40.
 !
 ! longitude of center in degrees =  55.
 ! latitude of center in degrees =  29.
 !
 ! angle of rotation of the first chunk =  0.
 !
 ! corner  1
 ! longitude in degrees =  33.491330699023123
 ! latitude in degrees =  8.437884305357894
 !
 ! corner  2
 ! longitude in degrees =  76.508669300976862
 ! latitude in degrees =  8.437884305357894
 !
 ! corner  3
 ! longitude in degrees =  24.3197569739134
 ! latitude in degrees =  44.884677931822168
 !
 ! corner  4
 ! longitude in degrees =  85.680243026086586
 ! latitude in degrees =  44.884677931822168
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =  1024
 ! GLL points along a great circle =  4096
 ! average distance between points in degrees =  1.533980831E-3
 ! average distance between points in km =  9.77299118
 ! average size of a spectral element in km =  39.0919647
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
 ! size of static arrays per slice =  144.08386400000001  MB
 !                                 =  137.40908050537109  MiB
 !                                 =  0.14408386400000001  GB
 !                                 =  0.13418855518102646  GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =  9221.3672960000004  MB
 !                                      =  8794.18115234375  MiB
 !                                      =  9.2213672960000004  GB
 !                                      =  8.5880675315856934  GiB
 !                                      =  9.22136729599999996E-3  TB
 !                                      =  8.38678469881415367E-3  TiB
 !

 integer, parameter :: NEX_XI_VAL =  128
 integer, parameter :: NEX_ETA_VAL =  128

 integer, parameter :: NSPEC_CRUST_MANTLE =  2816
 integer, parameter :: NSPEC_OUTER_CORE =  280
 integer, parameter :: NSPEC_INNER_CORE =  12

 integer, parameter :: NGLOB_CRUST_MANTLE =  192897
 integer, parameter :: NGLOB_OUTER_CORE =  20905
 integer, parameter :: NGLOB_INNER_CORE =  1053

 integer, parameter :: NSPECMAX_ANISO_IC =  1

 integer, parameter :: NSPECMAX_ISO_MANTLE =  2816
 integer, parameter :: NSPECMAX_TISO_MANTLE =  2816
 integer, parameter :: NSPECMAX_ANISO_MANTLE =  1

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =  2816
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =  12

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =  2816
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =  12

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

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =  2816
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =  280

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =  192897

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
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =  280

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .true.

 integer, parameter :: NPROC_XI_VAL =  8
 integer, parameter :: NPROC_ETA_VAL =  8
 integer, parameter :: NCHUNKS_VAL =  1
 integer, parameter :: NPROCTOT_VAL =  64

 integer, parameter :: ATT1_VAL =  5
 integer, parameter :: ATT2_VAL =  5
 integer, parameter :: ATT3_VAL =  5
 integer, parameter :: ATT4_VAL =  2816
 integer, parameter :: ATT5_VAL =  12

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =  320
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =  320
 integer, parameter :: NSPEC2D_BOTTOM_CM =  16
 integer, parameter :: NSPEC2D_TOP_CM =  256
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =  6
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =  6
 integer, parameter :: NSPEC2D_BOTTOM_IC =  4
 integer, parameter :: NSPEC2D_TOP_IC =  4
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =  84
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =  84
 integer, parameter :: NSPEC2D_BOTTOM_OC =  4
 integer, parameter :: NSPEC2D_TOP_OC =  16
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

 integer, parameter :: NGLOB_XY_CM =  192897
 integer, parameter :: NGLOB_XY_IC =  1

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 integer, parameter :: NT_DUMP_ATTENUATION =  100000000

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    40.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    45.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    29.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =    55.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

