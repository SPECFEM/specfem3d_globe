 
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
 ! maximum number of points per region =      1559049
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      4677147
 !
 ! total elements per slice =        26964
 ! total points per slice =      1788935
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !    3909600.00000000     
 ! approximate total number of points in entire mesh = 
 !    259482445.000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    734419035.000000     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          960
 ! GLL points along a great circle =         3840
 ! average distance between points in degrees =   9.3750000E-02
 ! average distance between points in km =    10.42452    
 ! average size of a spectral element in km =    41.69810    
 !
 ! number of time steps =         5600
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
 ! size of static arrays per slice =    636.357580000000       MB
 !                                 =    606.877880096436       MiB
 !                                 =   0.636357580000000       GB
 !                                 =   0.592654179781675       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    95.4536370000000       GB
 !                                      =    88.8981269672513       GiB
 !                                      =   9.545363699999999E-002  TB
 !                                      =   8.681457711645635E-002  TiB
 !
 
 integer, parameter :: NEX_XI_VAL =          240
 integer, parameter :: NEX_ETA_VAL =          240
 
 integer, parameter :: NSPEC_CRUST_MANTLE =        23616
 integer, parameter :: NSPEC_OUTER_CORE =         2160
 integer, parameter :: NSPEC_INNER_CORE =         1188
 
 integer, parameter :: NGLOB_CRUST_MANTLE =      1559049
 integer, parameter :: NGLOB_OUTER_CORE =       146761
 integer, parameter :: NGLOB_INNER_CORE =        83125
 
 integer, parameter :: NSPECMAX_ANISO_IC =            1
 
 integer, parameter :: NSPECMAX_ISO_MANTLE =        23616
 integer, parameter :: NSPECMAX_TISO_MANTLE =        23616
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUAT =        23616
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =         1188
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =        23616
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =         1188
 
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
 
 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =      1559049
 
 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.
 
 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.
 
 logical, parameter :: ATTENUATION_VAL = .true.
 
 logical, parameter :: ATTENUATION_3D_VAL = .false.
 
 logical, parameter :: ELLIPTICITY_VAL = .true.
 
 logical, parameter :: GRAVITY_VAL = .true.
 
 logical, parameter :: OCEANS_VAL = .true.
 
 integer, parameter :: NX_BATHY_VAL = NX_BATHY
 integer, parameter :: NY_BATHY_VAL = NY_BATHY
 
 logical, parameter :: ROTATION_VAL = .true.
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =         2160
 
 integer, parameter :: NGLOB1D_RADIAL_CM =          153
 integer, parameter :: NGLOB1D_RADIAL_OC =           97
 integer, parameter :: NGLOB1D_RADIAL_IC =           13
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_CM =        16934
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_OC =         4518
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_IC =         3411
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_CM =        16934
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_OC =         4518
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_IC =         3411
 integer, parameter :: NPROC_XI_VAL =            5
 integer, parameter :: NPROC_ETA_VAL =            5
 integer, parameter :: NCHUNKS_VAL =            6
 integer, parameter :: NPROCTOT_VAL =          150
 integer, parameter :: NGLOB2DMAX_XY_CM_VAL =        16934
 integer, parameter :: NGLOB2DMAX_XY_OC_VAL =         4518
 integer, parameter :: NGLOB2DMAX_XY_IC_VAL =         3411
 integer, parameter :: NUMMSGS_FACES_VAL =           60
 integer, parameter :: NCORNERSCHUNKS_VAL =            8
 integer, parameter :: ATT1 =            5
 integer, parameter :: ATT2 =            5
 integer, parameter :: ATT3 =            5
 integer, parameter :: ATT4 =        23616
 integer, parameter :: ATT5 =         1188
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          864
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          864
 integer, parameter :: NSPEC2D_BOTTOM_CM =          144
 integer, parameter :: NSPEC2D_TOP_CM =         2304
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =          198
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =          198
 integer, parameter :: NSPEC2D_BOTTOM_IC =           36
 integer, parameter :: NSPEC2D_TOP_IC =           36
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          216
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          216
 integer, parameter :: NSPEC2D_BOTTOM_OC =           36
 integer, parameter :: NSPEC2D_TOP_OC =          144
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1
 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY = .false.
 logical, parameter :: COMPUTE_AND_STORE_STRAIN = .true. 
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 1
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 1
