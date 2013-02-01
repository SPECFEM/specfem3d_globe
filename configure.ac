dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl autoreconf -I m4
dnl
dnl You must have recent versions of Autoconf and Automake installed.


AC_PREREQ(2.61)
AC_INIT([Specfem3D Globe], [6.0.0], [jtromp AT princeton.edu], [Specfem3DGlobe])
AC_CONFIG_SRCDIR([src/specfem3D/specfem3D.F90])
AC_CONFIG_HEADER([setup/config.h])
AC_CONFIG_MACRO_DIR([m4])


# 'configure' options

AC_ARG_ENABLE([double-precision],
    [AC_HELP_STRING([--enable-double-precision],
        [solver in double precision @<:@default=no@:>@])],
    [want_double_precision="$enableval"],
    [want_double_precision=no])
if test "$want_double_precision" = no; then
    CUSTOM_REAL=SIZE_REAL
    CUSTOM_MPI_TYPE=MPI_REAL
else
    CUSTOM_REAL=SIZE_DOUBLE
    CUSTOM_MPI_TYPE=MPI_DOUBLE_PRECISION
fi
AC_SUBST([CUSTOM_REAL])
AC_SUBST([CUSTOM_MPI_TYPE])

AC_ARG_WITH([cuda],
    [AC_HELP_STRING([--with-cuda],
        [build cuda GPU enabled version @<:@default=no@:>@])],
    [want_cuda="$withval"],
    [want_cuda=no])
#AM_CONDITIONAL([COND_CUDA], [test "$want_cuda" = yes])
AM_CONDITIONAL([COND_CUDA], [test "$want_cuda" != no])
AM_CONDITIONAL([COND_CUDA5], [test "$want_cuda" = cuda5])

AC_ARG_WITH([vtk],
    [AC_HELP_STRING([--with-vtk],
        [build VTK enabled version @<:@default=no@:>@])],
    [want_vtk="$withval"],
    [want_vtk=no])
AM_CONDITIONAL([COND_VTK], [test "$want_vtk" != no])

# Checks for programs.

# a courtesy to the installed base of users
if test x"$FC" = x && test x"$F90" != x; then
    FC="$F90"
fi
if test x"$MPIFC" = x && test x"$MPIF90" != x; then
    MPIFC="$MPIF90"
fi

AC_PROG_FC
export FC
export MPIFC
F77="$FC"
FFLAGS="$FCFLAGS"
AC_PROVIDE([AC_PROG_F77])
AC_SUBST([FCENV])

flags_guess="$SHELL $srcdir/flags.guess"
AC_MSG_NOTICE([running $flags_guess])
flags=`$flags_guess` ||
    AC_MSG_ERROR([$flags_guess failed])
eval $flags

AC_FC_WRAPPERS
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM()

AC_PROG_CC

AC_LANG_PUSH(C)
AC_CHECK_HEADER(emmintrin.h,AC_DEFINE([HAVE_EMMINTRIN],[1],[Define if emmintrin.h]))
AC_CHECK_HEADER(xmmintrin.h,AC_DEFINE([HAVE_XMMINTRIN],[1],[Define if xmmintrin.h]))
AC_LANG_POP(C)

# influential environment variables

AC_ARG_VAR(MPIFC, [MPI Fortran compiler command])
AC_ARG_VAR(MPILIBS, [extra libraries for linking MPI programs])
AC_ARG_VAR(FLAGS_CHECK, [Fortran compiler flags for non-critical subroutines])
AC_ARG_VAR(FLAGS_NO_CHECK, [Fortran compiler flags for creating fast, production-run code for critical subroutines])

AC_ARG_VAR(CUDA_LIB,[Location of CUDA library libcudart])
AC_ARG_VAR(CUDA_INC,[Location of CUDA include files])
AC_ARG_VAR(MPI_INC,[Location of MPI include mpi.h, which is needed by nvcc when compiling cuda files])

AC_ARG_VAR(AR, [ar library creation])
AC_ARG_VAR(ARFLAGS, [ar flag library creation])
AC_ARG_VAR(RANLIB, [ranlib library creation])

if test x"$MPIFC" = x; then
    MPIFC=mpif90
fi

AC_ARG_VAR(LOCAL_PATH_IS_ALSO_GLOBAL, [files on a local path on each node are also seen as global with same path @<:@default=true@:>@])
if test x"$LOCAL_PATH_IS_ALSO_GLOBAL" = x; then
    LOCAL_PATH_IS_ALSO_GLOBAL=true
fi



# Checks for libraries.


# Checks for header files.

CIT_FC_MPI_HEADER([$MPIFC], [$FLAGS_NO_CHECK])

if test "$want_cuda" = yes; then
    CIT_CUDA_COMPILER
fi

# Checks for typedefs, structures, and compiler characteristics.


# Checks for library functions.

AC_CONFIG_FILES([
  Makefile 
  src/auxiliaries/Makefile 
  src/meshfem3D/Makefile 
  src/specfem3D/Makefile 
  src/create_header_file/Makefile 
  setup/constants.h setup/precision.h])
AC_OUTPUT

#daniel: adding custom lines to config.h.in through autoheader
AH_BOTTOM([
/* Uncomment and define to select optimized file i/o for regional simulations */
/* map fails when output files are > 4GB, which is often the case for GPU simulations */
// #define USE_MAP_FUNCTION
])

dnl end of configure.ac
