
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
 ! maximum number of points per region =       576013
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      1728039
 !
 ! total elements per slice =         9360
 ! total points per slice =       626599
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !    224640.00000000000     
 ! approximate total number of points in entire mesh = 
 !    15038376.000000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    42811848.000000000     
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    90.000000    
 ! angular size in second direction in degrees =    90.000000    
 !
 ! longitude of center in degrees =    0.0000000    
 ! latitude of center in degrees =    90.000000    
 !
 ! angle of rotation of the first chunk =    0.0000000    
 !
 ! corner            1
 ! longitude in degrees =   -45.000000000000000     
 ! latitude in degrees =    35.446575249587262     
 !
 ! corner            2
 ! longitude in degrees =    45.000000000000000     
 ! latitude in degrees =    35.446575249587262     
 !
 ! corner            3
 ! longitude in degrees =   -135.00000000000000     
 ! latitude in degrees =    35.446575249587262     
 !
 ! corner            4
 ! longitude in degrees =    135.00000000000000     
 ! latitude in degrees =    35.446575249587262     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          256
 ! GLL points along a great circle =         1024
 ! average distance between points in degrees =   0.35156250    
 ! average distance between points in km =    39.091965    
 ! average size of a spectral element in km =    156.36786    
 !
 ! number of time steps =         7900
 !
 ! number of seismic sources =            1
 !

 ! approximate static memory needed by the solver:
 ! ----------------------------------------------
 !
 ! size of static arrays per slice =   8.13465118408203125E-002  GB
 !
 !   (should be below and typically equal to 80% of 1.5 GB = 1.2 GB on pangu
 !    at Caltech, and below and typically equal to 85% of 2 GB = 1.7 GB
 !    on Marenostrum in Barcelona)
 !   (if significantly more, the job will not run by lack of memory)
 !   (if significantly less, you waste a significant amount of memory)
 !
 ! size of static arrays for all slices =   0.32538604736328125       GB
 !                                      =   3.17759811878204346E-004  TB
 !

 integer, parameter :: NEX_XI_VAL =           64
 integer, parameter :: NEX_ETA_VAL =           64

 integer, parameter :: NSPEC_CRUST_MANTLE =         8640
 integer, parameter :: NSPEC_OUTER_CORE =          688
 integer, parameter :: NSPEC_INNER_CORE =           32

 integer, parameter :: NGLOB_CRUST_MANTLE =       576013
 integer, parameter :: NGLOB_OUTER_CORE =        47985
 integer, parameter :: NGLOB_INNER_CORE =         2601

 integer, parameter :: NSPECMAX_ANISO_IC =            1

 integer, parameter :: NSPECMAX_ISO_MANTLE =         8640
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

 logical, parameter :: ROTATION_VAL = .false.
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =            1

 integer, parameter :: NGLOB1D_RADIAL_CM =          109
 integer, parameter :: NGLOB1D_RADIAL_OC =           65
 integer, parameter :: NGLOB1D_RADIAL_IC =            9
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_CM =         8574
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_OC =         2134
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_IC =          178
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_CM =         8574
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_OC =         2134
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_IC =          178
 integer, parameter :: NPROC_XI_VAL =            2
 integer, parameter :: NPROC_ETA_VAL =            2
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =            4
 integer, parameter :: NGLOB2DMAX_XY_VAL_CM =         8574
 integer, parameter :: NGLOB2DMAX_XY_VAL_OC =         2134
 integer, parameter :: NGLOB2DMAX_XY_VAL_IC =          178
 integer, parameter :: NUMMSGS_FACES_VAL =            2
 integer, parameter :: NCORNERSCHUNKS_VAL =            1
 integer, parameter :: ATT1 =            1
 integer, parameter :: ATT2 =            1
 integer, parameter :: ATT3 =            1
 integer, parameter :: ATT4 =            1
 integer, parameter :: ATT5 =            1
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          440
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          440
 integer, parameter :: NSPEC2D_BOTTOM_CM =           64
 integer, parameter :: NSPEC2D_TOP_CM =         1024
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =            8
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =            8
 integer, parameter :: NSPEC2D_BOTTOM_IC =           16
 integer, parameter :: NSPEC2D_TOP_IC =           16
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          100
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          100
 integer, parameter :: NSPEC2D_BOTTOM_OC =           16
 integer, parameter :: NSPEC2D_TOP_OC =           64
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1
