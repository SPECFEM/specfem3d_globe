dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl autoreconf -i
dnl
dnl You must have recent versions of Autoconf and Automake installed.

############################################################

AC_PREREQ(2.61)
AC_INIT([Specfem3D_Globe], m4_normalize(m4_include([VERSION])), [see the wiki], [Specfem3DGlobe])
AC_CONFIG_SRCDIR([src/specfem3D/specfem3D.F90])

AC_CONFIG_HEADERS([setup/config.h])

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
    [AS_HELP_STRING([--enable-debug],
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
    [AS_HELP_STRING([--enable-vectorization],
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
### static compilation
###
AC_ARG_ENABLE([static-compilation],
    [AS_HELP_STRING([--enable-static-compilation],
        [build with static compilation option enabled @<:@default=no@:>@])],
    [want_static="$enableval"],
    [want_static=no])
AM_CONDITIONAL([COND_STATIC_COMPILATION], [test x"$want_static" != xno])

###
### GPU version
###

## CUDA
AC_ARG_WITH([cuda],
    [AS_HELP_STRING([--with-cuda],
        [build CUDA GPU enabled version @<:@default=no@:>@])],
    [want_cuda="$withval"],
    [want_cuda=no])

# sets default: --with-cuda
if test x"$want_cuda" = xyes; then
    want_cuda=cuda4
fi
AM_CONDITIONAL([COND_CUDA], [test x"$want_cuda" != xno])
AM_CONDITIONAL([COND_CUDA4], [test x"$want_cuda" = xcuda4])
AM_CONDITIONAL([COND_CUDA5], [test x"$want_cuda" = xcuda5])
AM_CONDITIONAL([COND_CUDA6], [test x"$want_cuda" = xcuda6])
AM_CONDITIONAL([COND_CUDA7], [test x"$want_cuda" = xcuda7])
AM_CONDITIONAL([COND_CUDA8], [test x"$want_cuda" = xcuda8])
AM_CONDITIONAL([COND_CUDA9], [test x"$want_cuda" = xcuda9])
AM_CONDITIONAL([COND_CUDA10], [test x"$want_cuda" = xcuda10])
AM_CONDITIONAL([COND_CUDA11], [test x"$want_cuda" = xcuda11])
AM_CONDITIONAL([COND_CUDA12], [test x"$want_cuda" = xcuda12])

# cuda linking for cuda 5x and 6x and 7x and 8x and ..
AM_CONDITIONAL([COND_CUDA_PLUS],
 [test "$want_cuda" = cuda4 \
    -o "$want_cuda" = cuda5 \
    -o "$want_cuda" = cuda6 \
    -o "$want_cuda" = cuda7 \
    -o "$want_cuda" = cuda8 \
    -o "$want_cuda" = cuda9 \
    -o "$want_cuda" = cuda10 \
    -o "$want_cuda" = cuda11 \
    -o "$want_cuda" = cuda12 \]
)

# CUDA-aware MPI setting
AC_ARG_ENABLE([cuda-aware-mpi],
    [AS_HELP_STRING([--enable-cuda-aware-mpi],
        [build CUDA-aware MPI enabled version @<:@default=no@:>@])],
    [enable_cuda_aware_mpi="$enableval"],
    [enable_cuda_aware_mpi=no])
AC_ARG_WITH([cuda-aware-mpi],
    [AS_HELP_STRING([--with-cuda-aware-mpi],
        [build CUDA-aware MPI enabled version @<:@default=no@:>@])],
    [with_cuda_aware_mpi="$withval"],
    [with_cuda_aware_mpi=no])
AS_IF([test x"$enable_cuda_aware_mpi" != xno || test x"$with_cuda_aware_mpi" != xno],
    [want_cuda_aware_mpi=yes],[want_cuda_aware_mpi=no])
AM_CONDITIONAL([COND_CUDA_AWARE_MPI], [test x"$want_cuda_aware_mpi" != xno])

## OpenCL
AC_ARG_WITH([opencl],
 [AS_HELP_STRING([--with-opencl], [build OpenCL GPU enabled version @<:@default=no@:>@])],
    [want_opencl="$withval"],
    [want_opencl=no])
AM_CONDITIONAL([COND_OCL], [test x"$want_opencl" != xno])

## HIP
AC_ARG_WITH([hip],
 [AS_HELP_STRING([--with-hip], [build HIP GPU enabled version @<:@default=no@:>@])],
    [want_hip="$withval"],
    [want_hip=no])
AM_CONDITIONAL([COND_HIP], [test x"$want_hip" != xno])
AM_CONDITIONAL([COND_HIP_MI8], [test x"$want_hip" = xMI8])
AM_CONDITIONAL([COND_HIP_MI25], [test x"$want_hip" = xMI25])
AM_CONDITIONAL([COND_HIP_MI50], [test x"$want_hip" = xMI50])
AM_CONDITIONAL([COND_HIP_MI100], [test x"$want_hip" = xMI100])
AM_CONDITIONAL([COND_HIP_MI250], [test x"$want_hip" = xMI250])

AM_CONDITIONAL([COND_HIP_CUDA5], [test x"$want_hip" = xcuda5])
AM_CONDITIONAL([COND_HIP_CUDA6], [test x"$want_hip" = xcuda6])
AM_CONDITIONAL([COND_HIP_CUDA7], [test x"$want_hip" = xcuda7])
AM_CONDITIONAL([COND_HIP_CUDA8], [test x"$want_hip" = xcuda8])
AM_CONDITIONAL([COND_HIP_CUDA9], [test x"$want_hip" = xcuda9])
AM_CONDITIONAL([COND_HIP_CUDA10], [test x"$want_hip" = xcuda10])
AM_CONDITIONAL([COND_HIP_CUDA11], [test x"$want_hip" = xcuda11])
AM_CONDITIONAL([COND_HIP_CUDA12], [test x"$want_hip" = xcuda12])

###
### MIC (Xeon PHI)
###

AC_ARG_WITH([mic],
    [AS_HELP_STRING([--with-mic],
        [build native Xeon PHI enabled version @<:@default=no@:>@])],
    [want_mic="$withval"],
    [want_mic=no])
AM_CONDITIONAL([COND_MIC], [test x"$want_mic" != xno])

###
### LIBXSMM (Library for small matrix-matrix multiplications targeting Intel Architecture (x86)
###          see:https://github.com/hfp/libxsmm )
###
AC_ARG_WITH([xsmm],
    [AS_HELP_STRING([--with-xsmm],
        [build LIBXSMM (see:https://github.com/libxsmm/libxsmm) enabled version @<:@default=no@:>@])],
    [want_xsmm="$withval"],
    [want_xsmm=no])
AM_CONDITIONAL([COND_XSMM], [test x"$want_xsmm" != xno])

##
## support C++ Parallel STL routines (for vectorization & multithreading)
##
AC_ARG_WITH([parallel_stl],
 [AS_HELP_STRING([--with-parallel-stl], [build mesh sorting with C++ Parallel STL routines @<:@default=no@:>@])],
    [want_parallel_stl="$withval"],
    [want_parallel_stl=no])
AM_CONDITIONAL([COND_PARALLEL_STL], [test x"$want_parallel_stl" != xno])


###
### OpenMP
###

AC_ARG_WITH([openmp],
    [AS_HELP_STRING([--with-openmp],
        [build OpenMP enabled version @<:@default=no@:>@])],
    [with_omp="$withval"],
    [with_omp=no])
AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--enable-openmp],
        [build OpenMP enabled version @<:@default=no@:>@])],
    [enable_omp="$enableval"],
    [enable_omp=no])
AS_IF([test x"$enable_omp" != xno || test x"$with_omp" != xno],
    [want_omp=yes],[want_omp=no])
AM_CONDITIONAL([COND_OMP], [test x"$want_omp" != xno])

###
### VTK
###

AC_ARG_ENABLE([vtk],
    [AS_HELP_STRING([--enable-vtk],
        [build VTK enabled version @<:@default=no@:>@])],
    [want_vtk="$enableval"],
    [want_vtk=no])
AM_CONDITIONAL([COND_VTK], [test x"$want_vtk" != xno])

###
### ADIOS
###

AC_ARG_WITH([adios],
    [AS_HELP_STRING([--with-adios],
        [build ADIOS enabled version @<:@default=no@:>@])],
    [want_adios="$withval"],
    [want_adios=no])
AM_CONDITIONAL([COND_ADIOS], [test x"$want_adios" != xno])

###
### ADIOS2
###

AC_ARG_WITH([adios2],
    [AS_HELP_STRING([--with-adios2],
        [build ADIOS enabled version @<:@default=no@:>@])],
    [want_adios2="$withval"],
    [want_adios2=no])
AM_CONDITIONAL([COND_ADIOS2], [test x"$want_adios2" != xno])

###
### ASDF
###

AC_ARG_WITH([asdf],
    [AS_HELP_STRING([--with-asdf],
        [build ASDF enabled version @<:@default=no@:>@])],
    [want_asdf="$withval"],
    [want_asdf=no])
AM_CONDITIONAL([COND_ASDF], [test x$"$want_asdf" != xno])

###
### CEM
###

AC_ARG_WITH([cem],
    [AS_HELP_STRING([--with-cem],
        [build CEM enabled version @<:@default=no@:>@])],
    [want_cem="$withval"],
    [want_cem=no])

###
### IRIS EMC models
###

AC_ARG_WITH([emc],
    [AS_HELP_STRING([--with-emc],
        [build IRIS EMC models enabled version @<:@default=no@:>@])],
    [want_emc="$withval"],
    [want_emc=no])

###
### NetCDF support
###

AC_ARG_WITH([netcdf],
    [AS_HELP_STRING([--with-netcdf],
        [build NetCDF enabled version (IRIS EMC & CEM models) @<:@default=no@:>@])],
    [with_netcdf="$withval"],
    [with_netcdf=no])
AC_ARG_ENABLE([netcdf],
    [AS_HELP_STRING([--enable-netcdf],
        [build NetCDF enabled version (IRIS EMC & CEM models) @<:@default=no@:>@])],
    [enable_netcdf="$enableval"],
    [enable_netcdf=no])
AS_IF([test x"$enable_netcdf" != xno || test x"$with_netcdf" != xno],
    [want_netcdf=yes],[want_netcdf=no])

## turns on NetCDF support for CEM & EMC version
AS_IF([test x"$want_cem" != xno || test x"$want_emc" != xno],
    [want_netcdf=yes],[])

## turns on CEM & EMC models for NetCDF version
AS_IF([test x"$want_netcdf" != xno],
    [want_cem=yes])
AS_IF([test x"$want_netcdf" != xno],
    [want_emc=yes])

AM_CONDITIONAL([COND_CEM], [test x"$want_cem" != xno])
AM_CONDITIONAL([COND_EMC], [test x"$want_emc" != xno])
AM_CONDITIONAL([COND_NETCDF], [test x"$want_netcdf" != xno])

###
### PETSc support
###

AC_ARG_WITH([petsc],
    [AS_HELP_STRING([--with-petsc],
        [build PETSc enabled version @<:@default=no@:>@])],
    [want_petsc="$withval"],
    [want_petsc=no])
AM_CONDITIONAL([COND_PETSC], [test x$"$want_petsc" != xno])


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
AC_ARG_VAR(MPICC, [MPI C compiler command])

AC_ARG_VAR(FLAGS_CHECK, [Fortran compiler flags])

# scratch disks
AC_ARG_VAR(LOCAL_PATH_IS_ALSO_GLOBAL, [files on a local path on each node are also seen as global with same path @<:@default=true@:>@])

############################################################

# tests
AS_IF([test x"$MPIFC" = x],[
    MPIFC=mpif90
])

AS_IF([test x"$MPICC" = x],[
    MPICC=mpicc
])

AS_IF([test x"$LOCAL_PATH_IS_ALSO_GLOBAL" = x],[
    LOCAL_PATH_IS_ALSO_GLOBAL=true
])

# Checks for typedefs, structures, and compiler characteristics.
# check fortran modules flag
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
# AIX doesn't have err.h so we need to check
AC_CHECK_HEADER(err.h,      AC_DEFINE([HAVE_ERR],      [1],[Define if err.h]))
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
### ADIOS2
###

AS_IF([test x"$want_adios2" != xno], [
    AS_BOX([ADIOS2])

    # for fortran compilation
    CIT_ADIOS2_CONFIG

    # for C file compilation
    AC_MSG_CHECKING([for ADIOS2 C flags and libs])
    AS_IF([test x"$ADIOS2_CONFIG" != x], [
      ADIOS2_CFLAGS=`$ADIOS2_CONFIG --c-flags`
      ADIOS2_CLIBS=`$ADIOS2_CONFIG --c-libs`
      AC_MSG_RESULT(done)
    ],[AC_MSG_RESULT(no adios2_config)])
    AC_SUBST([ADIOS2_CFLAGS])
    AC_SUBST([ADIOS2_CLIBS])

    # in case adios2 c file compilation needs mpi
    AC_SUBST([MPICC])
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

## CUDA
AS_IF([test x"$want_cuda" != xno], [
    AS_BOX([CUDA])
    CIT_CUDA_CONFIG
])

## OpenCL
AS_IF([test x"$want_opencl" != xno], [
    AS_BOX([OpenCL])
    AC_ARG_VAR(OCL_CPU_FLAGS, [OpenCL CPU compiler flags])
    AC_ARG_VAR(OCL_GPU_FLAGS, [OpenCL GPU compiler flags])
    CIT_OPENCL_CONFIG
])

## HIP
AS_IF([test x"$want_hip" != xno], [
    AS_BOX([HIP])
    # environment variables
    AC_ARG_VAR(HIPCC, [AMD HIP compiler command])
    AC_ARG_VAR(HIP_FLAGS, [HIP compiler flags])
    AC_ARG_VAR(HIP_INC, [Location of HIP include files])
    AC_ARG_VAR(HIP_LIB, [Location of HIP library libhip_hcc])

    # tests HIPCC variable
    AS_IF([test x"$HIPCC" = x],[
        HIPCC=hipcc
    ])

    # Check for compiler
    # checks if program in path
    AC_PATH_PROG(HIPCC_PROG, $HIPCC)
    if test -z "$HIPCC_PROG" ; then
        AC_MSG_ERROR([cannot find '$HIPCC' program, please check your PATH.])
    fi

    # sets default HIP path
    AC_PATH_PROG(HIPCONFIG_PROG, hipconfig)
    AS_IF([test x"$HIPCONFIG_PROG" = x],[
      AC_PATH_PROG(HIPCONFIG_PROG, hipconfig, path=/opt/rocm/bin)
    ])
    AS_IF([test x"$HIPCONFIG_PROG" != x],[
        HIP_PATH=`$HIPCONFIG_PROG --path`
    ])
    AC_MSG_NOTICE([hipconfig program: $HIPCONFIG_PROG])
    AC_MSG_NOTICE([HIP path: $HIP_PATH])

    #AS_IF([test x"$HIP_INC" = x],[
    #    HIP_INC="${HIP_PATH}/include"
    #])

    # adds default HIP library
    AS_IF([test x"$HIP_LIB" = x -a x"${HIP_PATH}" != x],[
        HIP_LIB="${HIP_PATH}/lib"
        AC_MSG_NOTICE([adding default HIP library path: ${HIP_LIB}])
    ])

    # checks platform
    if test x"$HIPCONFIG_PROG" != x ; then
      GPU_PLATFORM=`$HIPCONFIG_PROG --platform`
      AS_IF([test x"$GPU_PLATFORM" = xnvidia], [
          AC_MSG_NOTICE([HIP PLATFORM NVIDIA detected.])
          HIP_PLATFORM=nvidia
          FLAG_PLATFORM="${FC_DEFINE}__HIP_PLATFORM_NVIDIA__"
      ])

      AS_IF([test x"$GPU_PLATFORM" = xamd], [
          AC_MSG_NOTICE([HIP PLATFORM AMD detected.])
          HIP_PLATFORM=amd
          FLAG_PLATFORM="${FC_DEFINE}__HIP_PLATFORM_AMD__"
      ])
      # compiler might require a platform selection flag
      case "${HIP_FLAGS}" in
        *HIP_PLATFORM*) AC_MSG_NOTICE([hip flags contain platform specifier: ${HIP_FLAGS}]) ;;
        *) AC_MSG_NOTICE([adding flag for platform specifier ${FLAG_PLATFORM}]);
           HIP_FLAGS="${HIP_FLAGS} ${FLAG_PLATFORM}"
           AC_MSG_NOTICE([HIP_FLAGS: ${HIP_FLAGS}]) ;;
      esac
    fi
    # sets default platform to AMD if not set
    AS_IF([test x"$GPU_PLATFORM" = x],[
      GPU_PLATFORM=amd
    ])

    # for compilation errors like:
    #   /usr/bin/ld: obj/assemble_MPI_scalar_gpu.hip.o: relocation R_X86_64_32 against `.rodata.str1.1' can not be used
    #      when making a PIE object; recompile with -fPIC
    # try to add
    #   HIP_FLAGS="-fPIC"
    # or
    #   HIP_FLAGS="-fno-gpu-rdc"
    # compiler might require a -fPIC flag
    case "${HIP_FLAGS}" in
      *fPIC*) AC_MSG_NOTICE([HIP flags contain position independent code flag -fPIC: ${HIP_FLAGS}]);;
      *) AC_MSG_NOTICE([consider adding compilation flag for position independent code: HIP_FLAGS="-fPIC"]);;
    esac

    # test lib & include
    AS_IF([test x"$HIP_LIB" != x],[
        HIP_LDFLAGS="-L$HIP_LIB"
    ])
    AS_IF([test x"$HIP_INC" != x],[
        HIP_CPPFLAGS="-I$HIP_INC"
    ])
    HIP_LIBS="-lamdhip64"

    # checks header file
    AC_LANG_PUSH(C)
    AC_REQUIRE_CPP
    CFLAGS_save="$CFLAGS"
    LDFLAGS_save="$LDFLAGS"
    LIBS_save="$LIBS"

    # uses hipcc compiler
    CFLAGS="$HIP_FLAGS"
    if test "x$HIP_INC" != "x"; then
      HIP_CPPFLAGS="-I$HIP_INC"
      CFLAGS="$CFLAGS $HIP_CPPFLAGS"
    fi

    AC_MSG_CHECKING([for hip/hip_runtime.h])
    ac_compile='$HIPCC -c $CFLAGS conftest.$ac_ext >&5'
    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([[#include <hip/hip_runtime.h>]],[[void* ptr = 0;]])
    ], [
        AC_MSG_RESULT(yes)
    ], [
        AC_MSG_RESULT(no)
        AC_MSG_ERROR([HIP runtime header not found; try setting HIP_INC.])
    ])

    # checks for HIP library
    if test "x$HIP_LIB" != "x"; then
      HIP_LDFLAGS="-L$HIP_LIB"
      LDFLAGS="$HIP_LDFLAGS $LDFLAGS"
    fi
    HIP_LIBS="-lamdhip64"
    LIBS="$HIP_LIBS $LIBS"

    # runs compilation test with hipcc
    AC_MSG_CHECKING([hipcc compilation with hipMalloc in -lamdhip64])
    ac_compile='$HIPCC -c $CFLAGS conftest.$ac_ext >&5'
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([[#include <hip/hip_runtime.h>]],[[void* ptr = 0;hipMalloc(&ptr, 1);]])
    ], [
      AC_MSG_RESULT(yes)
    ], [
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([HIP library function with hipcc compilation failed; try setting HIP_INC.])
    ])

    # runs linking test with hipcc
    AC_MSG_CHECKING([hipcc linking with hipMalloc in -lamdhip64])
    ac_link='$HIPCC -o conftest$ac_exeext $CFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
    AC_LINK_IFELSE([
      AC_LANG_PROGRAM([[
#include <stdio.h>
#include <hip/hip_runtime.h>]],[[void* ptr = 0;hipMalloc(&ptr, 1);]])
    ],[
      AC_MSG_RESULT(yes)
    ],[
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([HIP library linking with HIP failed; try setting HIP_LIB.])
    ])

    #HIP_HEADER_H="hip/hip_runtime.h"
    #AC_CHECK_HEADER([${HIP_HEADER_H}], [], [
    #    AC_MSG_ERROR([HIP header not found; try setting HIP_INC.])
    #])

    CFLAGS="$CFLAGS_save"
    LDFLAGS="$LDFLAGS_save"
    LIBS="$LIBS_save"
    AC_LANG_POP(C)

    # export flags
    AC_SUBST([HIP_CPPFLAGS])
    AC_SUBST([HIP_LDFLAGS])
    AC_SUBST([HIP_LIBS])
])
AM_CONDITIONAL([COND_HIP_PLATFORM_AMD], [test x"$GPU_PLATFORM" = xamd])
AM_CONDITIONAL([COND_HIP_PLATFORM_NVIDIA], [test x"$GPU_PLATFORM" = xnvidia])

# CUDA-aware MPI setting
AS_IF([test x"$want_cuda_aware_mpi" != xno], [
    AC_MSG_NOTICE([enabling CUDA-aware MPI support])
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

##
## Parallel STL
##
AS_IF([test x"$want_parallel_stl" != xno], [
  AS_BOX([Parallel STL])
  AC_MSG_NOTICE([C++ Parallel STL compilation is enabled for mesh sorting])
  AC_ARG_VAR(PSTL_LIB, [Parallel STL library directory])
])


###
### OpenMP
###

AS_IF([test x"$want_omp" != xno], [
  AS_BOX([OpenMP])
  AC_MSG_NOTICE([OpenMP compilation is enabled])

  AC_ARG_VAR(OMP_FCFLAGS, [OpenMP Fortran compiler flags])
  AC_ARG_VAR(OMP_LIB, [Location of extra OpenMP libraries])

  # openmp checking
  AC_MSG_NOTICE([OpenMP flag $OMP_FCFLAGS])
  CIT_FC_OPENMP_MODULE([$FC], [$OMP_FCFLAGS])
])


###
### MPI
###

AS_IF([test x"$want_cuda" != xno -o x"$want_opencl" != xno -o x"$want_hip" != xno -o x"$want_vtk" != xno], [
    AS_BOX([MPI])
    # checks MPI include directory
    CIT_MPI_INCDIR([$MPIFC])
])

###
### NetCDF support
###

AS_IF([test x"$want_netcdf" != xno], [
    AS_BOX([NetCDF])
    AC_MSG_NOTICE([NetCDF is enabled])

    AC_ARG_VAR(NETCDF_INC, [NetCDF include directory])
    AC_ARG_VAR(NETCDF_FCFLAGS, [Fortran compiler flags])
    AC_ARG_VAR(NETCDF_LIBS, [extra libraries for linking programs])

    # NetCDF checking for Fortran
    AC_LANG_PUSH(Fortran)
    AC_MSG_CHECKING([for NetCDF modules])
    FCFLAGS="$NETCDF_FCFLAGS $FCFLAGS"
    # adds include directory
    if test "x$NETCDF_INC" != "x"; then
      NETCDF_INCLUDES="-I$NETCDF_INC"
      FCFLAGS+=" $NETCDF_INCLUDES"
    fi
    # adds CEM flags if given
    if test "x$CEM_FCFLAGS" != "x"; then
      FCFLAGS+=" $CEM_FCFLAGS"
    fi
    # adds CEM flags if given
    if test "x$EMC_FCFLAGS" != "x"; then
      FCFLAGS+=" $EMC_FCFLAGS"
    fi
    # compilation test
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([], [[
      use netcdf
      ]])
    ], [
      AC_MSG_RESULT(yes)
    ], [
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([NetCDF module not found; is NetCDF built with Fortran support for this compiler? try setting NETCDF_INC.])
    ])
    AC_LANG_POP([Fortran])
    AC_SUBST([NETCDF_INCLUDES])
])

###
### CEM
###

AS_IF([test x"$want_cem" != xno], [
    AS_BOX([Comprehensive Earth Model (CEM)])
    AC_MSG_NOTICE([CEM is enabled])

    AC_ARG_VAR(CEM_FCFLAGS, [CEM Fortran compiler flags])
    AC_ARG_VAR(CEM_LIBS, [CEM extra libraries for linking programs])

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

###
### IRIS EMC models
###

AS_IF([test x"$want_emc" != xno], [
    AS_BOX([IRIS EMC models])
    AC_MSG_NOTICE([IRIS EMC is enabled])

    AC_ARG_VAR(EMC_FCFLAGS, [IRIS EMC Fortran compiler flags])
    AC_ARG_VAR(EMC_LIBS, [IRIS EMC extra libraries for linking programs])

    # NetCDF checking for Fortran
    #AC_LANG_PUSH(Fortran)
    #AC_MSG_CHECKING([for NetCDF modules])
    #FCFLAGS="$EMC_FCFLAGS $FCFLAGS"
    #AC_COMPILE_IFELSE([
    #  AC_LANG_PROGRAM([], [[
    #  use netcdf
    #  ]])
    #], [
    #  AC_MSG_RESULT(yes)
    #], [
    #  AC_MSG_RESULT(no)
    #  AC_MSG_ERROR([NetCDF module not found; is NetCDF built with Fortran support for this compiler?])
    #])
    #AC_LANG_POP([Fortran])
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
### static compilation
###

AS_IF([test x"$want_static" != xno],[
  AS_BOX([STATIC COMPILATION])
  AC_MSG_NOTICE([static compilation is enabled])
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
### PETSC
###
AS_IF([test x"$want_petsc" != xno], [
    AS_BOX([PETSc])
    AC_MSG_NOTICE([PETSc is enabled])

    AC_ARG_VAR(PETSC_INC, [PETSC include directory])
    AC_ARG_VAR(PETSC_FCFLAGS, [Fortran compiler flags (PETSC specific flags)])
    AC_ARG_VAR(PETSC_LIB, [extra library directory for linking programs with PETSC])

    # test include & lib
    AS_IF([test x"$PETSC_INC" != x],[
        PETSC_FCFLAGS+=" -I$PETSC_INC"
    ])
    AS_IF([test x"$PETSC_LIB" != x],[
        PETSC_LDFLAGS="-L$PETSC_LIB"
    ])
    PETSC_LIBS="-lpetsc"

    # PETSC checking for Fortran compilation
    AC_LANG_PUSH(Fortran)
    FCFLAGS_save="$FCFLAGS"
    LDFLAGS_save="$LDFLAGS"
    LIBS_save="$LIBS"

    # compilation test
    AC_MSG_CHECKING([for PETSC modules])
    if test "x$PETSC_FCFLAGS" != "x"; then
      FCFLAGS="$PETSC_FCFLAGS $FCFLAGS"
    fi
    # we'll need PETSc's Krylov subspace method solver
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([], [[
! PETSc's Krylov subspace method (KSP) solver
#include "petsc/finclude/petscksp.h"
  use petscksp
      ]])
    ], [
      AC_MSG_RESULT(yes)
    ], [
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([PETSC module not found; is PETSC built with Fortran support for this compiler? try setting PETSC_INC.])
    ])

    # runs linking test
    AC_MSG_CHECKING([PETSc linking with PetscInitialize in -lpetsc])
    if test "x$PETSC_LDFLAGS" != "x"; then
      LDFLAGS="$PETSC_LDFLAGS $LDFLAGS"
    fi
    LIBS="$PETSC_LIBS $LIBS"
    AC_LINK_IFELSE([
      AC_LANG_PROGRAM([], [[
! PETSc's Krylov subspace method (KSP) solver
#include "petsc/finclude/petscksp.h"
  use petscksp
  implicit none
  PetscErrorCode :: ierr
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call PetscFinalize(ierr)
      ]])
    ], [
      AC_MSG_RESULT(yes)
    ],[
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([PETSc library linking failed; try setting PETSC_LIB.])
    ])

    # restore original flags
    FCFLAGS="$FCFLAGS_save"
    LDFLAGS="$LDFLAGS_save"
    LIBS="$LIBS_save"
    AC_LANG_POP([Fortran])

    # export PETSC flags
    AC_SUBST([PETSC_FCFLAGS])
    AC_SUBST([PETSC_LDFLAGS])
    AC_SUBST([PETSC_LIBS])
])

##
## Git version info
##

AS_BOX([git version info])
CIG_PKG_GIT(SPECFEM3D_GLOBE)

AC_DEFINE_UNQUOTED([SPECFEM3D_GLOBE_VERSION], ["$PACKAGE_VERSION"], [Define SPECFEM3D_GLOBE version])
#AC_DEFINE_UNQUOTED([SPECFEM3D_GLOBE_DOI], ["10.5281/zenodo.xxxxx"], [Define SPECFEM3D_GLOBE doi])

AS_IF([test x"$SPECFEM3D_GLOBE_RELEASE_VERSION" = x0],[
  AC_MSG_NOTICE([building from git repository])
  GIT_PACKAGE_VERSION=$SPECFEM3D_GLOBE_GIT_REVISION
  GIT_COMMIT_VERSION=$SPECFEM3D_GLOBE_GIT_HASH
  GIT_DATE_VERSION=$SPECFEM3D_GLOBE_GIT_DATE
],[
  AC_MSG_NOTICE([not a git repository])
  GIT_PACKAGE_VERSION=$PACKAGE_VERSION
  GIT_COMMIT_VERSION="unknown"
  GIT_DATE_VERSION="unknown"
])
AC_SUBST([GIT_PACKAGE_VERSION])
AC_SUBST([GIT_COMMIT_VERSION])
AC_SUBST([GIT_DATE_VERSION])

############################################################

AS_BOX([setting up default simulation setup])

# Output results.
AC_CONFIG_FILES([
    Makefile
    setup/constants.h
    setup/constants_tomography.h
    setup/precision.h
    setup/config.fh
    setup/version.fh
])

AC_CONFIG_FILES([
    DATA/Par_file:DATA/Par_file
    DATA/CMTSOLUTION:DATA/CMTSOLUTION
    DATA/STATIONS:DATA/STATIONS
])

AC_CONFIG_LINKS([
    DATA/CCREM:DATA/CCREM
    DATA/cemRequest:DATA/cemRequest
    DATA/crust1.0:DATA/crust1.0
    DATA/crust2.0:DATA/crust2.0
    DATA/crustmap:DATA/crustmap
    DATA/epcrust:DATA/epcrust
    DATA/eucrust-07:DATA/eucrust-07
    DATA/full_sphericalharmonic_model:DATA/full_sphericalharmonic_model
    DATA/gladm25:DATA/gladm25
    DATA/heterogen:DATA/heterogen
    DATA/IRIS_EMC:DATA/IRIS_EMC
    DATA/Lebedev_sea99:DATA/Lebedev_sea99
    DATA/Montagner_model:DATA/Montagner_model
    DATA/mars:DATA/mars
    DATA/moon:DATA/moon
    DATA/old:DATA/old
    DATA/PPM:DATA/PPM
    DATA/QRFSI12:DATA/QRFSI12
    DATA/s362ani:DATA/s362ani
    DATA/s20rts:DATA/s20rts
    DATA/s40rts:DATA/s40rts
    DATA/sglobe:DATA/sglobe
    DATA/spiral1.4:DATA/spiral1.4
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
