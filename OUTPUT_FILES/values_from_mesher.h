
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 ! these statistics include the central cube
 !
 ! number of processors =           150
 !
 ! number of ES nodes =     18.75000    
 ! percentage of total 640 ES nodes =     2.929688      %
 ! total memory available on these ES nodes (Gb) =     300.0000    
 !
 ! max points in largest region =        496525
 ! max vector length =       1489575
 ! min vector length =            25
 ! min critical vector length =            75
 !
 ! on ES and SX-5, make sure "loopcnt=" parameter
 ! in Makefile is greater than       1489575
 !
 ! total elements per AB slice =          7856
 ! total points per AB slice =        532047
 !
 ! total elements per AC slice =          8192
 ! total points per AC slice =        553999
 !
 ! total elements per BC slice =          8864
 ! total points per BC slice =        597455
 !
 ! load balancing AB/BC for points =     89.05223      %
 ! load balancing AB/BC for elements =     88.62816      %
 !
 ! load balancing AC/BC for points =     92.72648      %
 ! load balancing AC/BC for elements =     92.41877      %
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !       1205600
 ! approximate total number of points in entire mesh = 
 !     81517845.00000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !     222868835.0000000     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =           640
 ! GLL points along a great circle =          2560
 ! average distance between points in degrees =    0.1406250    
 ! average distance between points in km =     15.63679    
 !

 integer, parameter :: NEX_XI_VAL =           160
 integer, parameter :: NEX_ETA_VAL =           160

 integer, parameter :: NSPEC_CRUST_MANTLE_AB =          6464
 integer, parameter :: NSPEC_CRUST_MANTLE_AC =          6784
 integer, parameter :: NSPEC_CRUST_MANTLE_BC =          7424
 integer, parameter :: NSPEC_OUTER_CORE_AB =          1024
 integer, parameter :: NSPEC_OUTER_CORE_AC =          1040
 integer, parameter :: NSPEC_OUTER_CORE_BC =          1072
 integer, parameter :: NSPEC_INNER_CORE =           368

 integer, parameter :: NGLOB_CRUST_MANTLE_AB =        434317
 integer, parameter :: NGLOB_CRUST_MANTLE_AC =        455181
 integer, parameter :: NGLOB_CRUST_MANTLE_BC =        496525
 integer, parameter :: NGLOB_OUTER_CORE_AB =         70853
 integer, parameter :: NGLOB_OUTER_CORE_AC =         71941
 integer, parameter :: NGLOB_OUTER_CORE_BC =         74053
 integer, parameter :: NGLOB_INNER_CORE =         26877

 integer, parameter :: NSPECMAX_CRUST_MANTLE = NSPEC_CRUST_MANTLE_BC
 integer, parameter :: NGLOBMAX_CRUST_MANTLE = NGLOB_CRUST_MANTLE_BC

 integer, parameter :: NSPECMAX_OUTER_CORE = NSPEC_OUTER_CORE_BC
 integer, parameter :: NGLOBMAX_OUTER_CORE = NGLOB_OUTER_CORE_BC

 integer, parameter :: NSPECMAX_ANISO_IC = 1
 integer, parameter :: NSPECMAX_ISO_MANTLE = NSPECMAX_CRUST_MANTLE
 integer, parameter :: NSPECMAX_TISO_MANTLE =          3328
 integer, parameter :: NSPECMAX_ANISO_MANTLE = 1

 integer, parameter :: NSPECMAX_CRUST_MANTLE_ATTENUAT = NSPECMAX_CRUST_MANTLE
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION = NSPEC_INNER_CORE

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

 logical, parameter :: ANISOTROPIC_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .true.

 logical, parameter :: ATTENUATION_VAL_3D = .false.

 logical, parameter :: ELLIPTICITY_VAL = .true.

 logical, parameter :: GRAVITY_VAL = .true.

 logical, parameter :: ROTATION_VAL = .true.
 integer, parameter :: NSPECMAX_OUTER_CORE_ROTATION = NSPECMAX_OUTER_CORE

