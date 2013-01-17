dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl autoreconf -I m4
dnl
dnl You must have recent versions of Autoconf and Automake installed.


AC_PREREQ(2.61)
AC_INIT([Specfem3D Globe], [5.1.5], [jtromp AT princeton.edu], [Specfem3DGlobe])
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
AC_ARG_VAR(FLAGS_CHECK, [Fortran compiler flags])
if test x"$MPIFC" = x; then
    MPIFC=mpif90
fi

AC_ARG_VAR(LOCAL_PATH_IS_ALSO_GLOBAL, [files on a local path on each node are also seen as global with same path @<:@default=true@:>@])
if test x"$LOCAL_PATH_IS_ALSO_GLOBAL" = x; then
    LOCAL_PATH_IS_ALSO_GLOBAL=true
fi



# Checks for libraries.


# Checks for header files.

CIT_FC_MPI_HEADER([$MPIFC], [$FLAGS_CHECK])


# Checks for typedefs, structures, and compiler characteristics.

CIT_FC_MODULE_EXTENSION
AC_SUBST([FC_MODEXT])
CIT_FC_MODULE_FLAG([], [])
CIT_FC_MODULE_OUTPUT_FLAG([
	FCFLAGS_f90="$FC_MODOUT./obj $FC_MODINC./obj $FC_MODINC. $FCFLAGS_f90"
	FC_MODDIR=./obj
], [
	FC_MODDIR=.
])
AC_SUBST([FC_MODDIR])

# Checks for library functions.


AC_CONFIG_FILES([
    Makefile
    setup/constants.h
    setup/precision.h
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
