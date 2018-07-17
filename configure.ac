dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl autoreconf -i
dnl
dnl You must have recent versions of Autoconf and Automake installed.

############################################################

AC_PREREQ(2.61)
AC_INIT([Specfem3D_Globe], [7.0.2], [see the wiki], [Specfem3DGlobe])
AC_CONFIG_SRCDIR([src/specfem3D/specfem3D.F90])

AC_CONFIG_HEADER([setup/config.h])

AC_CONFIG_MACRO_DIR([m4])
m4_include(m4/cit_backports.m4)

############################################################

AS_BOX([setting up compilation flags])
AC_CANONICAL_HOST

# 'configure' options

###
### precision
###

AC_ARG_ENABLE([double-precision],
    [AS_HELP_STRING([--enable-double-precision],
        [solver in double precision @<:@default=no@:>@])],
    [want_double_precision="$enableval"],
    [want_double_precision=no])
if test x"$want_double_precision" = xno; then
    CUSTOM_REAL=SIZE_REAL
    CUSTOM_MPI_TYPE=MPI_REAL
else
    CUSTOM_REAL=SIZE_DOUBLE
    CUSTOM_MPI_TYPE=MPI_DOUBLE_PRECISION
fi
AC_SUBST([CUSTOM_REAL])
AC_SUBST([CUSTOM_MPI_TYPE])

###
### debugging
###

AC_ARG_ENABLE([debug],
    [AC_HELP_STRING([--enable-debug],
        [build with debugging options enabled @<:@default=no@:>@])],
    [want_debug="$enableval"],
    [want_debug=no])
AM_CONDITIONAL([COND_DEBUG], [test x"$want_debug" != xno])
export COND_DEBUG_FALSE
export COND_DEBUG_TRUE

###
### FORCE_VECTORIZATION
###

AC_ARG_ENABLE([vectorization],
    [AC_HELP_STRING([--enable-vectorization],
        [build FORCE_VECTORIZATION enabled version @<:@default=auto@:>@])],
    [want_vec="$enableval"],
    [want_vec=no])
if test x"$want_vec" == xauto; then
    if test x"$want_debug" != xno; then
        want_vec=no
    else
        want_vec=yes
    fi
else
    if test x"$want_vec" != xno && test x"$want_debug" != xno; then
        AC_MSG_ERROR([--enable-debug and --enable-vectorization cannot be specified simultaneously.])
    fi
fi
AM_CONDITIONAL([COND_VECTORIZATION], [test x"$want_vec" != xno])

###
### CUDA
###

AC_ARG_WITH([cuda],
    [AC_HELP_STRING([--with-cuda],
        [build CUDA GPU enabled version @<:@default=no@:>@])],
    [want_cuda="$withval"],
    [want_cuda=no])
AM_CONDITIONAL([COND_CUDA], [test x"$want_cuda" != xno])
AM_CONDITIONAL([COND_CUDA5], [test x"$want_cuda" = xcuda5])
AM_CONDITIONAL([COND_CUDA6], [test x"$want_cuda" = xcuda6])
AM_CONDITIONAL([COND_CUDA7], [test x"$want_cuda" = xcuda7])
AM_CONDITIONAL([COND_CUDA8], [test x"$want_cuda" = xcuda8])
# cuda linking for cuda 5x and 6x
AM_CONDITIONAL([COND_CUDA_PLUS], [test "$want_cuda" = cuda5 -o "$want_cuda" = cuda6 -o "$want_cuda" = cuda7 -o "$want_cuda" = cuda8])

AC_ARG_WITH([opencl],
 [AC_HELP_STRING([--with-opencl], [build OpenCL GPU enabled version @<:@default=no@:>@])],
    [want_opencl="$withval"],
    [want_opencl=no])
AM_CONDITIONAL([COND_OCL], [test "$want_opencl" != no])

###
### MIC (Xeon PHI)
###

AC_ARG_WITH([mic],
    [AC_HELP_STRING([--with-mic],
        [build native Xeon PHI enabled version @<:@default=no@:>@])],
    [want_mic="$withval"],
    [want_mic=no])
AM_CONDITIONAL([COND_MIC], [test "$want_mic" != no])

###
### LIBXSMM (Library for small matrix-matrix multiplications targeting Intel Architecture (x86)
###          see:https://github.com/hfp/libxsmm )
###
AC_ARG_WITH([xsmm],
    [AC_HELP_STRING([--with-xsmm],
        [build LIBXSMM (see:https://github.com/hfp/libxsmm) enabled version @<:@default=no@:>@])],
    [want_xsmm="$withval"],
    [want_xsmm=no])
AM_CONDITIONAL([COND_XSMM], [test "$want_xsmm" != no])

###
### OpenMP
###

AC_ARG_ENABLE([openmp],
    [AC_HELP_STRING([--enable-openmp],
        [build OpenMP enabled version @<:@default=no@:>@])],
    [want_omp="$enableval"],
    [want_omp=no])
AM_CONDITIONAL([COND_OMP], [test "$want_omp" != no])


###
### VTK
###

AC_ARG_ENABLE([vtk],
    [AC_HELP_STRING([--enable-vtk],
        [build VTK enabled version @<:@default=no@:>@])],
    [want_vtk="$enableval"],
    [want_vtk=no])
AM_CONDITIONAL([COND_VTK], [test x"$want_vtk" != xno])

###
### ADIOS
###

AC_ARG_WITH([adios],
    [AC_HELP_STRING([--with-adios],
        [build ADIOS enabled version @<:@default=no@:>@])],
    [want_adios="$withval"],
    [want_adios=no])
AM_CONDITIONAL([COND_ADIOS], [test x"$want_adios" != xno])

###
### ASDF
###

AC_ARG_WITH([asdf],
    [AC_HELP_STRING([--with-asdf],
        [build ASDF enabled version @<:@default=no@:>@])],
    [want_asdf="$withval"],
    [want_asdf=no])
AM_CONDITIONAL([COND_ASDF], [test x$"$want_asdf" != xno])

###
### CEM
###

AC_ARG_WITH([cem],
    [AC_HELP_STRING([--with-cem],
        [build CEM enabled version @<:@default=no@:>@])],
    [want_cem="$withval"],
    [want_cem=no])
AM_CONDITIONAL([COND_CEM], [test x"$want_cem" != xno])

############################################################

# Checks for programs.

# a courtesy to the installed base of users
if test x"$FC" = x && test x"$F90" != x; then
    FC="$F90"
fi
if test x"$MPIFC" = x && test x"$MPIF90" != x; then
    MPIFC="$MPIF90"
fi
# note: if not set, default FCFLAGS (and CFLAGS) will be set by autoconf/configure script (mostly -g -O2)
if test x"$FCFLAGS" = x && test x"$FLAGS_CHECK" != x ; then
    FCFLAGS=""
fi

AC_PROG_FC
export FC
export MPIFC
F77="$FC"
FFLAGS="$FCFLAGS"

AC_PROVIDE([AC_PROG_F77])
AC_SUBST([FCENV])
AC_SUBST(srcdir)

flags_guess="$SHELL $srcdir/flags.guess"
AC_MSG_NOTICE([running $flags_guess])
flags=`$flags_guess` ||
    AC_MSG_ERROR([$flags_guess failed])
eval $flags

AC_FC_WRAPPERS
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM()

AC_FC_PP_DEFINE()
AC_SUBST([FC_DEFINE])
AC_FC_PP_SRCEXT(F90)dnl Because AC_FC_PP_DEFINE messes with things.

AC_PROG_CC

############################################################
# influential environment variables

AC_ARG_VAR(MPIFC, [MPI Fortran compiler command])
AC_ARG_VAR(MPILIBS, [extra libraries for linking MPI programs])

AC_ARG_VAR(FLAGS_CHECK, [Fortran compiler flags])

# scratch disks
AC_ARG_VAR(LOCAL_PATH_IS_ALSO_GLOBAL, [files on a local path on each node are also seen as global with same path @<:@default=true@:>@])

############################################################

# tests
AS_IF([test x"$MPIFC" = x],[
    MPIFC=mpif90
])

AS_IF([test x"$LOCAL_PATH_IS_ALSO_GLOBAL" = x],[
    LOCAL_PATH_IS_ALSO_GLOBAL=true
])

# Checks for typedefs, structures, and compiler characteristics.
AS_BOX([module extensions])

AC_FC_MODULE_EXTENSION
AC_SUBST([FC_MODEXT])
AC_FC_MODULE_FLAG([], [])
AC_FC_MODULE_OUTPUT_FLAG([
	FCFLAGS_f90="$FC_MODOUT./obj $FC_MODINC./obj $FC_MODINC. $FCFLAGS_f90"
	FC_MODDIR=./obj
], [
	FC_MODDIR=.
])
AC_SUBST([FC_MODDIR])

# Checks for header files.

AC_LANG_PUSH(C)
AC_CHECK_HEADER(emmintrin.h,AC_DEFINE([HAVE_EMMINTRIN],[1],[Define if emmintrin.h]))
AC_CHECK_HEADER(xmmintrin.h,AC_DEFINE([HAVE_XMMINTRIN],[1],[Define if xmmintrin.h]))
AC_LANG_POP(C)

# CIT_FC_MPI_MODULE([mpif.h], [$MPIFC], [$FLAGS_CHECK], [:], [
#     AC_MSG_FAILURE([cannot compile a trivial MPI program using $MPIFC])
# ])

# Checks for libraries.

###
### ADIOS
###

AS_IF([test x"$want_adios" != xno], [
    AS_BOX([ADIOS])
    CIT_ADIOS_CONFIG
])

###
### ASDF
###

AS_IF([test x"$want_asdf" != xno], [
    AS_BOX([ASDF])

    AC_ARG_VAR(ASDF_LIBS, [ASDF libraries for linking programs])

])

#mpbl: adding flags to config.h.in through autoheader in order to save them into the adios file.
#configure_flags_str="FC=$FC FCFLAGS=$FCFLAGS $FLAGS_CHECK"
#AC_DEFINE_UNQUOTED([CONFIGURE_FLAGS],["${configure_flags_str}"],
 #         [configuration flags to be saved in the adios output file.])
CONFIGURE_FLAGS="FC=$FC FCFLAGS=$FCFLAGS"
AC_SUBST([CONFIGURE_FLAGS])

###
### GPU
###

AS_IF([test x"$want_cuda" != xno], [
    AS_BOX([CUDA])
    CIT_CUDA_CONFIG
])

AS_IF([test x"$want_opencl" != xno], [
    AS_BOX([OpenCL])
    AC_ARG_VAR(OCL_CPU_FLAGS, [OpenCL CPU compiler flags])
    AC_ARG_VAR(OCL_GPU_FLAGS, [OpenCL GPU compiler flags])
    CIT_OPENCL_CONFIG
])

###
### MIC
###

AS_IF([test x"$want_mic" != xno], [
  AS_BOX([MIC])
  AC_MSG_NOTICE([MIC native compilation])
  FLAGS_CHECK=$(echo $FLAGS_CHECK | sed -e 's/\-xHost -fpe0\>//g')
])
if test x"$want_vec" == xno && test x"$want_mic" == xyes; then
    AC_MSG_NOTICE([--enable-vectorization and --with-mic should be both specified for better performance.])
fi

###
### LIBXSMM
###

AS_IF([test x"$want_xsmm" != xno], [
  AS_BOX([LIBXSMM])
  AC_MSG_NOTICE([LIBXSMM compilation is enabled])

  AC_ARG_VAR(LIBXSMM_INC, [LIBXSMM include directory])
  AC_ARG_VAR(LIBXSMM_LIBS, [LIBXSMM library directory])
  AC_ARG_VAR(BLAS_LIBS, [LAPACK/BLAS library directory])
])


###
### OpenMP
###

AS_IF([test x"$want_omp" != xno], [
  AS_BOX([OpenMP])
  AC_MSG_NOTICE([OpenMP compilation is enabled])

  AC_ARG_VAR(OMP_FCFLAGS, [OpenMP Fortran compiler flags])

  # openmp checking
  AC_MSG_NOTICE([OpenMP flag $OMP_FCFLAGS])
  CIT_FC_OPENMP_MODULE([$FC], [$OMP_FCFLAGS])
])

###
### FORCE_VECTORIZATION
###

AS_IF([test x"$want_vec" != xno],[
  AS_BOX([FORCE VECTORIZATION])
  AC_MSG_NOTICE([FORCE_VECTORIZATION is enabled])

  FLAGS_CHECK="${FLAGS_CHECK} ${FC_DEFINE}FORCE_VECTORIZATION"
  CPPFLAGS="${CPPFLAGS} -DFORCE_VECTORIZATION"
])

###
### VTK
###

AS_IF([test x"$want_vtk" != xno], [
    AS_BOX([VTK])
    CIT_OPTIONS_VTK
    CIT_PATH_VTK
    # Add optional libraries we need:
    VTK_LIBS="-lvtkRendering -lvtkGraphics $VTK_LIBS"
])

###
### MPI
###

AS_IF([test x"$want_cuda" != xno -o x"$want_opencl" != xno -o x"$want_vtk" != xno], [
    AS_BOX([MPI])
    # checks MPI include directory
    CIT_MPI_INCDIR([$MPIFC])
])

###
### CEM
###

AS_IF([test x"$want_cem" != xno], [
    AS_BOX([Comprehensive Earth Model (CEM)])
    AC_MSG_NOTICE([CEM is enabled])

    AC_ARG_VAR(CEM_FCFLAGS, [CEM Fortran compiler flags])
    AC_ARG_VAR(CEM_LIBS, [CEM extra libraries for linking programs])

    # adds compiler flag
    FLAGS_CHECK="${FLAGS_CHECK} ${FC_DEFINE}CEM"

    # daniel: todo add NetCDF checking for Fortran, something like:
    # netcdf checking
    #AC_LANG_PUSH(Fortran)
    #AC_MSG_CHECKING([for NetCDF modules])
    #FCFLAGS="$CEM_FCFLAGS $FCFLAGS"
    #AC_COMPILE_IFELSE([
    #  AC_LANG_PROGRAM([], [[
    #  use netcdf
    #  ]])
    #], [
    #  AC_MSG_RESULT(yes)
    #], [
    #  AC_MSG_RESULT(no)
    #  AC_MSG_ERROR([NetCDF modules not found; is NetCDF built with Fortran support for this compiler?])
    #])
    #AC_LANG_POP([Fortran])
])


############################################################

AS_BOX([setting up default simulation setup])

# Output results.
AC_CONFIG_FILES([
    Makefile
    setup/constants.h
    setup/constants_tomography.h
    setup/precision.h
    setup/config.fh
])
AS_IF([test -d $srcdir/.git], [
    GIT_VERSION_DEPS="$srcdir/.git/logs/HEAD"
    AC_MSG_NOTICE([building from git repository])
    AC_SUBST(GIT_VERSION_DEPS)
], [
    dnl Only substitute if not in a git repository.
    AC_MSG_NOTICE([not a git repository])
    AC_CONFIG_FILES([setup/version.fh])
])
AC_CONFIG_FILES([
    DATA/Par_file:DATA/Par_file
    DATA/CMTSOLUTION:DATA/CMTSOLUTION
    DATA/STATIONS:DATA/STATIONS
])
AC_CONFIG_LINKS([
    DATA/crust2.0:DATA/crust2.0
    DATA/crustmap:DATA/crustmap
    DATA/epcrust:DATA/epcrust
    DATA/eucrust-07:DATA/eucrust-07
    DATA/heterogen:DATA/heterogen
    DATA/Lebedev_sea99:DATA/Lebedev_sea99
    DATA/Montagner_model:DATA/Montagner_model
    DATA/old:DATA/old
    DATA/PPM:DATA/PPM
    DATA/QRFSI12:DATA/QRFSI12
    DATA/s20rts:DATA/s20rts
    DATA/s362ani:DATA/s362ani
    DATA/s40rts:DATA/s40rts
    DATA/Simons_model:DATA/Simons_model
    DATA/topo_bathy:DATA/topo_bathy
    DATA/Zhao_JP_model:DATA/Zhao_JP_model
])
AC_CONFIG_COMMANDS([bin], [AS_MKDIR_P(bin)])
AC_CONFIG_COMMANDS([obj], [AS_MKDIR_P(obj)])
AC_CONFIG_COMMANDS([OUTPUT_FILES], [AS_MKDIR_P(OUTPUT_FILES)])

AC_OUTPUT

#daniel: adding custom lines to config.h.in through autoheader
AH_BOTTOM([
/* Define to select optimized file i/o for regional simulations */
/* map fails when output files are > 4GB, which is often the case for GPU simulations */
#undef USE_MAP_FUNCTION
])


dnl FYI, this is not a "result", but AS_ECHO is not available in 2.61, and it
dnl basically works like a portable echo.

AC_MSG_RESULT([])
AS_BOX([$PACKAGE_NAME $PACKAGE_VERSION])
AC_MSG_RESULT([])
AC_MSG_RESULT([$0 has completed and set up a default configuration to build.])
AC_MSG_RESULT([])
AC_MSG_RESULT([You may wish to modify the following files before building:])
AC_MSG_RESULT([AS_HELP_STRING([DATA/Par_file],
                              [Set parameters affecting the build and simulation before running make and compiling.])])
AC_MSG_RESULT([AS_HELP_STRING([DATA/CMTSOLUTION],
                              [Set the source parameters before running the solver.])])
AC_MSG_RESULT([AS_HELP_STRING([DATA/STATIONS],
                              [Set the receiver stations before running the solver.])])
AC_MSG_RESULT([])

dnl end of configure.ac
