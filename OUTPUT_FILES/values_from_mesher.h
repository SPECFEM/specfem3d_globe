
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 ! these statistics include the central cube
 !
 ! number of processors =           384
 !
 ! number of ES nodes =     48.00000    
 ! percentage of total 640 ES nodes =     7.500000      %
 ! total memory available on these ES nodes (Gb) =     768.0000    
 !
 ! max points in largest region = max vector length =       3061529
 ! min vector length =            25
 ! min critical vector length =            75
 !
 ! on ES and SX-5, make sure "loopcnt=" parameter
 ! in Makefile is greater than       3061529
 !
 ! total elements per AB slice =         54208
 ! total points per AB slice =       3584995
 !
 ! total elements per AC slice =         55552
 ! total points per AC slice =       3671907
 !
 ! total elements per BC slice =         58240
 ! total points per BC slice =       3844835
 !
 ! load balancing AB/BC for points =     93.24184      %
 ! load balancing AB/BC for elements =     93.07692      %
 !
 ! load balancing AC/BC for points =     95.50233      %
 ! load balancing AC/BC for elements =     95.38461      %
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !      20193280
 ! approximate total number of points in entire mesh = 
 !     1336149371.000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !     3643879793.000000     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          2048
 ! GLL points along a great circle =          8192
 ! average distance between points in degrees =    4.3945313E-02
 ! average distance between points in km =     4.886496    
 !

 integer, parameter :: NEX_XI_VAL =           512
 integer, parameter :: NEX_ETA_VAL =           512

 integer, parameter :: NSPEC_CRUST_MANTLE_AB =         42752
 integer, parameter :: NSPEC_CRUST_MANTLE_AC =         44032
 integer, parameter :: NSPEC_CRUST_MANTLE_BC =         46592
 integer, parameter :: NSPEC_OUTER_CORE_AB =          7040
 integer, parameter :: NSPEC_OUTER_CORE_AC =          7104
 integer, parameter :: NSPEC_OUTER_CORE_BC =          7232
 integer, parameter :: NSPEC_INNER_CORE =          4416

 integer, parameter :: NGLOB_CRUST_MANTLE_AB =       2814233
 integer, parameter :: NGLOB_CRUST_MANTLE_AC =       2896921
 integer, parameter :: NGLOB_CRUST_MANTLE_BC =       3061529
 integer, parameter :: NGLOB_OUTER_CORE_AB =        469109
 integer, parameter :: NGLOB_OUTER_CORE_AC =        473333
 integer, parameter :: NGLOB_OUTER_CORE_BC =        481653
 integer, parameter :: NGLOB_INNER_CORE =        301653

 integer, parameter :: NSPECMAX_CRUST_MANTLE = NSPEC_CRUST_MANTLE_BC
 integer, parameter :: NGLOBMAX_CRUST_MANTLE = NGLOB_CRUST_MANTLE_BC

 integer, parameter :: NSPECMAX_OUTER_CORE = NSPEC_OUTER_CORE_BC
 integer, parameter :: NGLOBMAX_OUTER_CORE = NGLOB_OUTER_CORE_BC

 integer, parameter :: NSPECMAX_ANISO_IC = 1
 integer, parameter :: NSPECMAX_ISO_MANTLE = NSPECMAX_CRUST_MANTLE
 integer, parameter :: NSPECMAX_TISO_MANTLE =         16384
 integer, parameter :: NSPECMAX_ANISO_MANTLE = 1

 integer, parameter :: NSPECMAX_CRUST_MANTLE_ATTENUAT = NSPECMAX_CRUST_MANTLE
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION = NSPEC_INNER_CORE

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

 logical, parameter :: ANISOTROPIC_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .true.

 logical, parameter :: ELLIPTICITY_VAL = .true.

 logical, parameter :: GRAVITY_VAL = .true.

 logical, parameter :: ROTATION_VAL = .true.
 integer, parameter :: NSPECMAX_OUTER_CORE_ROTATION = NSPECMAX_OUTER_CORE

