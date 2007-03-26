 
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
 ! number of processors =          600
 !
 ! number of ES nodes =    75.00000    
 ! percentage of total 640 ES nodes =    11.71875      %
 ! total memory available on these ES nodes (Gb) =    1200.000    
 !
 ! max points in largest region =      1581573
 ! max vector length =      4744719
 !
 ! on ES and SX-5, make sure "loopcnt=" parameter
 ! in Makefile is greater than      4744719
 !
 ! total elements per AB slice =        28368
 ! total points per AB slice =      1893923
 !
 ! total elements per AC slice =        29124
 ! total points per AC slice =      1942979
 !
 ! total elements per BC slice =        30636
 ! total points per BC slice =      2040419
 !
 ! load balancing AB/BC for points =    92.82030      %
 ! load balancing AB/BC for elements =    92.59695      %
 !
 ! load balancing AC/BC for points =    95.22451      %
 ! load balancing AC/BC for elements =    95.06463      %
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !     16545600
 ! approximate total number of points in entire mesh = 
 !    1105476595.00000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    2969287785.00000     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =         1920
 ! GLL points along a great circle =         7680
 ! average distance between points in degrees =   4.6875000E-02
 ! average distance between points in km =    5.212262    
 !
 ! number of time steps =        72100
 !
 ! number of seismic sources =            1
 !
 
 integer, parameter :: NEX_XI_VAL =          480
 integer, parameter :: NEX_ETA_VAL =          480
 
 integer, parameter :: NSPEC_CRUST_MANTLE_AB =        21744
 integer, parameter :: NSPEC_CRUST_MANTLE_AC =        22464
 integer, parameter :: NSPEC_CRUST_MANTLE_BC =        23904
 integer, parameter :: NSPEC_OUTER_CORE_AB =         4248
 integer, parameter :: NSPEC_OUTER_CORE_AC =         4284
 integer, parameter :: NSPEC_OUTER_CORE_BC =         4356
 integer, parameter :: NSPEC_INNER_CORE =         2376
 
 integer, parameter :: NGLOB_CRUST_MANTLE_AB =      1442181
 integer, parameter :: NGLOB_CRUST_MANTLE_AC =      1488837
 integer, parameter :: NGLOB_CRUST_MANTLE_BC =      1581573
 integer, parameter :: NGLOB_OUTER_CORE_AB =       286117
 integer, parameter :: NGLOB_OUTER_CORE_AC =       288517
 integer, parameter :: NGLOB_OUTER_CORE_BC =       293221
 integer, parameter :: NGLOB_INNER_CORE =       165625
 
 integer, parameter :: NSPECMAX_CRUST_MANTLE = NSPEC_CRUST_MANTLE_BC
 integer, parameter :: NGLOBMAX_CRUST_MANTLE = NGLOB_CRUST_MANTLE_BC
 
 integer, parameter :: NSPECMAX_OUTER_CORE = NSPEC_OUTER_CORE_BC
 integer, parameter :: NGLOBMAX_OUTER_CORE = NGLOB_OUTER_CORE_BC
 
 integer, parameter :: NSPECMAX_ANISO_IC = 1
 integer, parameter :: NSPECMAX_ISO_MANTLE = NSPECMAX_CRUST_MANTLE
 integer, parameter :: NSPECMAX_TISO_MANTLE =         9216
 integer, parameter :: NSPECMAX_ANISO_MANTLE = 1
 
 integer, parameter :: NSPECMAX_CRUST_MANTLE_ATTENUAT = NSPECMAX_CRUST_MANTLE
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION = NSPEC_INNER_CORE
 
 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.
 
 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.
 
 logical, parameter :: ATTENUATION_VAL = .true.
 
 logical, parameter :: ATTENUATION_VAL_3D = .false.
 
 logical, parameter :: ELLIPTICITY_VAL = .false.
 
 logical, parameter :: GRAVITY_VAL = .false.
 
 logical, parameter :: ROTATION_VAL = .false.
 integer, parameter :: NSPECMAX_OUTER_CORE_ROTATION = 1
 
