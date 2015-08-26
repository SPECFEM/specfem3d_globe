
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
 ! number of processors =  4
 !
 ! maximum number of points per region =  750417
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =  2251251
 !
 ! total elements per slice =  12432
 ! total points per slice =  832031
 !
 ! the time step of the solver will be DT =  0.126792654
 !
 ! total for full 1-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !  49728.
 ! approximate total number of points in entire mesh =
 !  3328124.
 ! approximate total number of degrees of freedom in entire mesh =
 !  9361516.
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =  20.
 ! angular size in second direction in degrees =  20.
 !
 ! longitude of center in degrees =  25.
 ! latitude of center in degrees =  40.
 !
 ! angle of rotation of the first chunk =  0.
 !
 ! corner  1
 ! longitude in degrees =  13.661882747982686
 ! latitude in degrees =  29.677869656290316
 !
 ! corner  2
 ! longitude in degrees =  36.338117252017319
 ! latitude in degrees =  29.677869656290316
 !
 ! corner  3
 ! longitude in degrees =  9.8824891302835312
 ! latitude in degrees =  49.192981807843722
 !
 ! corner  4
 ! longitude in degrees =  40.117510869716469
 ! latitude in degrees =  49.192981807843715
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =  1152
 ! GLL points along a great circle =  4608
 ! average distance between points in degrees =  1.363538438E-3
 ! average distance between points in km =  8.68710327
 ! average size of a spectral element in km =  34.7484131
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
 ! size of static arrays per slice =  390.46826399999998  MB
 !                                 =  372.37955474853516  MiB
 !                                 =  0.39046826400000001  GB
 !                                 =  0.36365190893411636  GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =  1561.8730559999999  MB
 !                                      =  1489.5182189941406  MiB
 !                                      =  1.561873056  GB
 !                                      =  1.4546076357364655  GiB
 !                                      =  1.56187305600000007E-3  TB
 !                                      =  1.42051526927389205E-3  TiB
 !

 integer, parameter :: NEX_XI_VAL =  64
 integer, parameter :: NEX_ETA_VAL =  64

 integer, parameter :: NSPEC_CRUST_MANTLE =  11264
 integer, parameter :: NSPEC_OUTER_CORE =  1120
 integer, parameter :: NSPEC_INNER_CORE =  48

 integer, parameter :: NGLOB_CRUST_MANTLE =  750417
 integer, parameter :: NGLOB_OUTER_CORE =  77857
 integer, parameter :: NGLOB_INNER_CORE =  3757

 integer, parameter :: NSPECMAX_ANISO_IC =  1

 integer, parameter :: NSPECMAX_ISO_MANTLE =  11264
 integer, parameter :: NSPECMAX_TISO_MANTLE =  11264
 integer, parameter :: NSPECMAX_ANISO_MANTLE =  1

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =  11264
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =  48

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =  11264
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =  48

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

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =  11264
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =  1120

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =  750417

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
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =  1120

 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .true.

 integer, parameter :: NPROC_XI_VAL =  2
 integer, parameter :: NPROC_ETA_VAL =  2
 integer, parameter :: NCHUNKS_VAL =  1
 integer, parameter :: NPROCTOT_VAL =  4

 integer, parameter :: ATT1_VAL =  5
 integer, parameter :: ATT2_VAL =  5
 integer, parameter :: ATT3_VAL =  5
 integer, parameter :: ATT4_VAL =  11264
 integer, parameter :: ATT5_VAL =  48

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =  640
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =  640
 integer, parameter :: NSPEC2D_BOTTOM_CM =  64
 integer, parameter :: NSPEC2D_TOP_CM =  1024
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =  12
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =  12
 integer, parameter :: NSPEC2D_BOTTOM_IC =  16
 integer, parameter :: NSPEC2D_TOP_IC =  16
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =  168
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =  168
 integer, parameter :: NSPEC2D_BOTTOM_OC =  16
 integer, parameter :: NSPEC2D_TOP_OC =  64
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

 integer, parameter :: NGLOB_XY_CM =  750417
 integer, parameter :: NGLOB_XY_IC =  1

 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.

 logical, parameter :: FORCE_VECTORIZATION_VAL = .true.

 integer, parameter :: NT_DUMP_ATTENUATION =  100000000

 double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL =    20.000000
 double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL =    40.000000
 double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL =    25.000000
 double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL =     0.000000

