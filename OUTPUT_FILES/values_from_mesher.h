 
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
 ! number of processors =          150
 !
 ! number of ES nodes =    18.75000    
 ! percentage of total 640 ES nodes =    2.929688      %
 ! total memory available on these ES nodes (Gb) =    300.0000    
 !
 ! max points in largest region =      1107613
 ! max vector length =      3322839
 ! min vector length =           25
 ! min critical vector length =           75
 !
 ! on ES and SX-5, make sure "loopcnt=" parameter
 ! in Makefile is greater than      3322839
 !
 ! total elements per AB slice =        18036
 ! total points per AB slice =      1207215
 !
 ! total elements per AC slice =        18792
 ! total points per AC slice =      1256271
 !
 ! total elements per BC slice =        20304
 ! total points per BC slice =      1353711
 !
 ! load balancing AB/BC for points =    89.17819      %
 ! load balancing AB/BC for elements =    88.82979      %
 !
 ! load balancing AC/BC for points =    92.80201      %
 ! load balancing AC/BC for elements =    92.55319      %
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !      2721600
 ! approximate total number of points in entire mesh = 
 !    182002045.000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    498295035.000000     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          960
 ! GLL points along a great circle =         3840
 ! average distance between points in degrees =   9.3750000E-02
 ! average distance between points in km =    10.42452    
 !
 ! number of time steps =        12100
 !
 ! number of seismic sources =            1
 !
 
 integer, parameter :: NEX_XI_VAL =          240
 integer, parameter :: NEX_ETA_VAL =          240
 
 integer, parameter :: NSPEC_CRUST_MANTLE_AB =        14544
 integer, parameter :: NSPEC_CRUST_MANTLE_AC =        15264
 integer, parameter :: NSPEC_CRUST_MANTLE_BC =        16704
 integer, parameter :: NSPEC_OUTER_CORE_AB =         2304
 integer, parameter :: NSPEC_OUTER_CORE_AC =         2340
 integer, parameter :: NSPEC_OUTER_CORE_BC =         2412
 integer, parameter :: NSPEC_INNER_CORE =         1188
 
 integer, parameter :: NGLOB_CRUST_MANTLE_AB =       968221
 integer, parameter :: NGLOB_CRUST_MANTLE_AC =      1014877
 integer, parameter :: NGLOB_CRUST_MANTLE_BC =      1107613
 integer, parameter :: NGLOB_OUTER_CORE_AB =       155869
 integer, parameter :: NGLOB_OUTER_CORE_AC =       158269
 integer, parameter :: NGLOB_OUTER_CORE_BC =       162973
 integer, parameter :: NGLOB_INNER_CORE =        83125
 
 integer, parameter :: NSPECMAX_CRUST_MANTLE = NSPEC_CRUST_MANTLE_BC
 integer, parameter :: NGLOBMAX_CRUST_MANTLE = NGLOB_CRUST_MANTLE_BC
 
 integer, parameter :: NSPECMAX_OUTER_CORE = NSPEC_OUTER_CORE_BC
 integer, parameter :: NGLOBMAX_OUTER_CORE = NGLOB_OUTER_CORE_BC
 
 integer, parameter :: NSPECMAX_ANISO_IC = 1
 integer, parameter :: NSPECMAX_ISO_MANTLE = NSPECMAX_CRUST_MANTLE
 integer, parameter :: NSPECMAX_TISO_MANTLE =            1
 integer, parameter :: NSPECMAX_ANISO_MANTLE = 1
 
 integer, parameter :: NSPECMAX_CRUST_MANTLE_ATTENUAT = 1
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION = 1
 
 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.
 
 logical, parameter :: ATTENUATION_VAL = .false.
 
 logical, parameter :: ATTENUATION_VAL_3D = .false.
 
 logical, parameter :: ELLIPTICITY_VAL = .false.
 
 logical, parameter :: GRAVITY_VAL = .false.
 
 logical, parameter :: ROTATION_VAL = .false.
 integer, parameter :: NSPECMAX_OUTER_CORE_ROTATION = 1
 
