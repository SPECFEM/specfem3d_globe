 
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
 ! number of processors =           24
 !
 ! maximum number of points per region =       791825
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      2375475
 !
 ! total elements per slice =        13275
 ! total points per slice =       887483
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !    313600.000000000     
 ! approximate total number of points in entire mesh = 
 !    20954987.0000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    59310609.0000000     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          320
 ! GLL points along a great circle =         1280
 ! average distance between points in degrees =   0.2812500    
 ! average distance between points in km =    31.27357    
 ! average size of a spectral element in km =    125.0943    
 !
 ! number of time steps =         8800
 !
 ! number of seismic sources =            1
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
 ! size of static arrays per slice =    130.248624000000       MB
 !                                 =    124.214767456055       MiB
 !                                 =   0.130248624000000       GB
 !                                 =   0.121303483843803       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    3125.96697600000       MB
 !                                      =    2981.15441894531       MiB
 !                                      =    3.12596697600000       GB
 !                                      =    2.91128361225128       GiB
 !                                      =   3.125966976000000E-003  TB
 !                                      =   2.843050402589142E-003  TiB
 !
 
 integer, parameter :: NEX_XI_VAL =           80
 integer, parameter :: NEX_ETA_VAL =           80
 
 integer, parameter :: NSPEC_CRUST_MANTLE =        11900
 integer, parameter :: NSPEC_OUTER_CORE =         1075
 integer, parameter :: NSPEC_INNER_CORE =          300
 
 integer, parameter :: NGLOB_CRUST_MANTLE =       791825
 integer, parameter :: NGLOB_OUTER_CORE =        74049
 integer, parameter :: NGLOB_INNER_CORE =        21609
 
 integer, parameter :: NSPECMAX_ANISO_IC =            1
 
 integer, parameter :: NSPECMAX_ISO_MANTLE =        11900
 integer, parameter :: NSPECMAX_TISO_MANTLE =            1
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUAT =            1
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =            1
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =            1
 
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
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =            1
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =            1
 
 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =            1
 
 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.
 
 logical, parameter :: ATTENUATION_VAL = .false.
 
 logical, parameter :: ATTENUATION_3D_VAL = .false.
 
 logical, parameter :: ELLIPTICITY_VAL = .false.
 
 logical, parameter :: GRAVITY_VAL = .false.
 
 logical, parameter :: OCEANS_VAL = .false.
 
 integer, parameter :: NX_BATHY_VAL = 1
 integer, parameter :: NY_BATHY_VAL = 1
 
 logical, parameter :: ROTATION_VAL = .false.
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =            1
 
 integer, parameter :: NGLOB1D_RADIAL_CM =          105
 integer, parameter :: NGLOB1D_RADIAL_OC =           73
 integer, parameter :: NGLOB1D_RADIAL_IC =            9
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_CM =         9914
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_OC =         2666
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_IC =         1083
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_CM =         9914
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_OC =         2666
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_IC =         1083
 integer, parameter :: NPROC_XI_VAL =            2
 integer, parameter :: NPROC_ETA_VAL =            2
 integer, parameter :: NCHUNKS_VAL =            6
 integer, parameter :: NPROCTOT_VAL =           24
 integer, parameter :: NGLOB2DMAX_XY_CM_VAL =         9914
 integer, parameter :: NGLOB2DMAX_XY_OC_VAL =         2666
 integer, parameter :: NGLOB2DMAX_XY_IC_VAL =         1083
 integer, parameter :: NUMMSGS_FACES_VAL =           24
 integer, parameter :: NCORNERSCHUNKS_VAL =            8
 integer, parameter :: ATT1 =            1
 integer, parameter :: ATT2 =            1
 integer, parameter :: ATT3 =            1
 integer, parameter :: ATT4 =            1
 integer, parameter :: ATT5 =            1
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          510
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          510
 integer, parameter :: NSPEC2D_BOTTOM_CM =          100
 integer, parameter :: NSPEC2D_TOP_CM =         1600
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           60
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           60
 integer, parameter :: NSPEC2D_BOTTOM_IC =           25
 integer, parameter :: NSPEC2D_TOP_IC =           25
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          130
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          135
 integer, parameter :: NSPEC2D_BOTTOM_OC =           25
 integer, parameter :: NSPEC2D_TOP_OC =          100
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1
 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY = .false.
 logical, parameter :: COMPUTE_AND_STORE_STRAIN = .false.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 1
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 1
