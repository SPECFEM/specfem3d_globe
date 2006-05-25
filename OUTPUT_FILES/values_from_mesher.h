
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 !
 ! number of chunks =  6
 !
 ! these statistics include the central cube
 !
 ! number of processors =  6
 !
 ! number of ES nodes =  0.75
 ! percentage of total 640 ES nodes =  0.1171875  %
 ! total memory available on these ES nodes (Gb) =  12.
 !
 ! max points in largest region =  1960621
 ! max vector length =  5881863
 !
 ! on ES and SX-5, make sure "loopcnt=" parameter
 ! in Makefile is greater than  5881863
 !
 ! total elements per AB slice =  30656
 ! total points per AB slice =  2036303
 !
 ! total elements per AC slice =  32000
 ! total points per AC slice =  2123215
 !
 ! total elements per BC slice =  34688
 ! total points per BC slice =  2296143
 !
 ! load balancing AB/BC for points =  88.68363  %
 ! load balancing AB/BC for elements =  88.37638  %
 !
 ! load balancing AC/BC for points =  92.468765  %
 ! load balancing AC/BC for elements =  92.25092  %
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !  192128
 ! approximate total number of points in entire mesh = 
 !  1.2731637E+7
 ! approximate total number of degrees of freedom in entire mesh = 
 !  3.4840163E+7
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =  256
 ! GLL points along a great circle =  1024
 ! average distance between points in degrees =  0.3515625
 ! average distance between points in km =  39.091965
 !
 ! number of time steps =  100
 !
 ! number of seismic sources =  1
 !

 integer, parameter :: NEX_XI_VAL =  64
 integer, parameter :: NEX_ETA_VAL =  64

 integer, parameter :: NSPEC_CRUST_MANTLE_AB =  25856
 integer, parameter :: NSPEC_CRUST_MANTLE_AC =  27136
 integer, parameter :: NSPEC_CRUST_MANTLE_BC =  29696
 integer, parameter :: NSPEC_OUTER_CORE_AB =  4096
 integer, parameter :: NSPEC_OUTER_CORE_AC =  4160
 integer, parameter :: NSPEC_OUTER_CORE_BC =  4288
 integer, parameter :: NSPEC_INNER_CORE =  704

 integer, parameter :: NGLOB_CRUST_MANTLE_AB =  1713325
 integer, parameter :: NGLOB_CRUST_MANTLE_AC =  1796013
 integer, parameter :: NGLOB_CRUST_MANTLE_BC =  1960621
 integer, parameter :: NGLOB_OUTER_CORE_AB =  273973
 integer, parameter :: NGLOB_OUTER_CORE_AC =  278197
 integer, parameter :: NGLOB_OUTER_CORE_BC =  286517
 integer, parameter :: NGLOB_INNER_CORE =  49005

 integer, parameter :: NSPECMAX_CRUST_MANTLE = NSPEC_CRUST_MANTLE_BC
 integer, parameter :: NGLOBMAX_CRUST_MANTLE = NGLOB_CRUST_MANTLE_BC

 integer, parameter :: NSPECMAX_OUTER_CORE = NSPEC_OUTER_CORE_BC
 integer, parameter :: NGLOBMAX_OUTER_CORE = NGLOB_OUTER_CORE_BC

 integer, parameter :: NSPECMAX_ANISO_IC = 1
 integer, parameter :: NSPECMAX_ISO_MANTLE = NSPECMAX_CRUST_MANTLE
 integer, parameter :: NSPECMAX_TISO_MANTLE =  1
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

