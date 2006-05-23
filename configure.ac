dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl ACLOCAL="aclocal -I m4" autoreconf
dnl
dnl You must have recent versions of Autoconf and Automake installed.

# $Id: configure.ac,v 1.1 2005/08/24 22:45:51 leif Exp $


AC_PREREQ(2.59)
AC_INIT([Specfem3D Globe], [3.5.0], [jtromp AT caltech.edu], [Specfem3DGlobe])
AC_CONFIG_SRCDIR([specfem3D.f90])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])


# 'configure' options

AC_ARG_WITH([pyre],
    [AC_HELP_STRING([--with-pyre],
        [build Pyrized version @<:@default=no@:>@])],
    [want_pyre="$withval"],
    [want_pyre=no])
AM_CONDITIONAL([COND_PYRE], [test "$want_pyre" = yes])

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

if test "$want_pyre" = yes; then
    AM_PATH_PYTHON([2.3])
    CIT_PYTHON_SYSCONFIG
fi

# a courtesy to the installed base of users
if test x"$FC" = x && test x"$F90" != x; then
    FC="$F90"
fi
if test x"$MPIFC" = x && test x"$MPIF90" != x; then
    MPIFC="$MPIF90"
fi

AC_PROG_FC
F77="$FC"
FFLAGS="$FCFLAGS"
AC_PROVIDE([AC_PROG_F77])

dnl Use AC_CANONICAL_BUILD (and package config.guess, etc.) in the future?
if test x"$UNAME_MS" = x; then
    UNAME_MACHINE=`(uname -m) 2>/dev/null` || UNAME_MACHINE=unknown
    UNAME_SYSTEM=`(uname -s) 2>/dev/null`  || UNAME_SYSTEM=unknown
    UNAME_MS="${UNAME_MACHINE}:${UNAME_SYSTEM}"
fi
AC_SUBST([FCENV])
case $UNAME_MS in
    i*86:Linux | x86_64:Linux)
        ################ PC Linux #################
        case $FC in
            pgf90|*/pgf90)
		#
		# Beowulf Portland pgf90
		#
                if test x"$FLAGS_CHECK" = x; then
                    FLAGS_CHECK="-O2 -Mbounds -Mneginfo -Mdclchk -Knoieee"
                fi
                if test x"$FLAGS_NO_CHECK" = x; then
                    FLAGS_NO_CHECK="-O2 -Mnobounds -Mneginfo -Mdclchk -Knoieee"
                fi
                ;;
            ifort|*/ifort)
		#
		# Intel ifort Fortran90 for Linux
		#
                if test x"$FLAGS_CHECK" = x; then
                    FLAGS_CHECK="\$(FLAGS_NO_CHECK)"
                fi
                if test x"$FLAGS_NO_CHECK" = x; then
                    # standard options
                    FLAGS_NO_CHECK="-O3 -e95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds"
                    # old hrothgar cluster at Caltech
                    #FLAGS_NO_CHECK = -O3 -tpp6 -xK -e95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds
                    # Pangu at Caltech
                    #FLAGS_NO_CHECK = $(IFORT_PROF) -vec_report0 -O2 -static -ip -xP -Wl,--allow-multiple-definition -L $$IFORT_ROOT/lib -limf -lirc
                    # debug with range checking
                    #FLAGS_NO_CHECK = -O0 -static -e95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check bounds
                fi
                #MPI_LIBS = -Vaxlib
                ;;
            g95|*/g95)
		#
		# g95 (free f95 compiler from http://www.g95.org, still under development, but works)
		#
                if test x"$FLAGS_CHECK" = x; then
                    FLAGS_CHECK="-O"
                fi
                if test x"$FLAGS_NO_CHECK" = x; then
                    FLAGS_NO_CHECK="-O"
                fi
                ;;
            f90|*/f90)
		#
		# AbSoft
		#
                if test x"$FLAGS_CHECK" = x; then
                    FLAGS_CHECK="-s -O2 -cpu:p7 -v -YDEALLOC=ALL"
                fi
                if test x"$FLAGS_NO_CHECK" = x; then
                    FLAGS_NO_CHECK="\$(FLAGS_CHECK)"
                fi
                ;;
            f95|*/f95)
		#
		# NAG compiler for Linux
		#
                if test x"$FLAGS_CHECK" = x; then
                    FLAGS_CHECK="-O -u -strict95 -C=all"
                fi
                if test x"$FLAGS_NO_CHECK" = x; then
                    FLAGS_NO_CHECK="-O -u -strict95"
                fi
                ;;
            lf95|*/lf95)
		#
		# Lahey f90
		#
                if test x"$FLAGS_CHECK" = x; then
                    FLAGS_CHECK="--warn --wo --tpp --f95 --dal -O --chk"
                fi
                if test x"$FLAGS_NO_CHECK" = x; then
                    FLAGS_NO_CHECK="--warn --wo --tpp --f95 --dal -O"
                fi
                ;;
        esac
        ;;
    *:IRIX*)
	################ SGI Irix #################
        case $FC in
            f90|*/f90)
                if test x"$MPIFC" = x; then
                    MPIFC=$FC
		    if test x"$MPILIBS" = x; then
                        MPILIBS="-lmpi -lfastm -lfpe"
                    fi
                fi
		if test x"$FLAGS_CHECK" = x; then
		    FLAGS_CHECK="\$(FLAGS_NO_CHECK) -check_bounds"
		fi
		if test x"$FLAGS_NO_CHECK" = x; then
		    FLAGS_NO_CHECK="-ansi -u -64 -O3 -OPT:Olimit=0 -OPT:roundoff=3 -OPT:IEEE_arithmetic=3 -r10000 -mips4"
		fi
                ;;
        esac
	##
	##  CAUTION: always define setenv TRAP_FPE OFF on SGI before compiling
	##
	FCENV="TRAP_FPE=OFF"
        ;;
    alpha:OSF1)
        ################## Compaq Dec Alpha #################
        case $FC in
            f90|*/f90)
                if test x"$MPIFC" = x; then
                    MPIFC=$FC
		    if test x"$MPILIBS" = x; then
                        MPILIBS="-lfmpi -lmpi"
                    fi
                fi
		if test x"$FLAGS_CHECK" = x; then
		    FLAGS_CHECK="\$(FLAGS_NO_CHECK) -check bounds"
		fi
		if test x"$FLAGS_NO_CHECK" = x; then
		    FLAGS_NO_CHECK="-fast -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nounderflow"
		fi
                ;;
        esac
        ;;
    SX-*:SUPER-UX | ES:ESOS)
        ################## Earth Simulator and NEC SX-5 ##################
        case $FC in
            f90|*/f90)
                if test x"$MPIFC" = x; then
                    MPIFC=$FC
                fi
		if test x"$FLAGS_CHECK" = x; then
		    FLAGS_CHECK="-C hopt -R2 -Wf\" -L nostdout noinclist mrgmsg noeject -msg b -pvctl loopcnt=14000000 expand=10 fullmsg vecthreshold=20 -s\" -pi auto line=100 exp=swap_all,rank"
		fi
		if test x"$FLAGS_NO_CHECK" = x; then
		    FLAGS_NO_CHECK="\$(FLAGS_CHECK)"
		fi
                ;;
        esac
        ;;
    *:AIX)
        ######## IBM ######
        case $FC in
            mpxlf*|*/mpxlf*)
                if test x"$MPIFC" = x; then
                    MPIFC=$FC
                fi
                ;;
        esac
        case $FC in
            *xlf*|*/*xlf*)
		if test x"$FLAGS_CHECK" = x; then
		    FLAGS_CHECK="-O3 -Q -qarch=auto -qcache=auto -qtune=auto -qlanglvl=95pure -qmaxmem=65536 -qflag=L:L -qhalt=L -qsuffix=f=f90"
                    #FLAGS_CHECK = -q32 -O3 -qarch=pwr4 -qlanglvl=95pure -qflag=L:L -qhalt=L -qsuffix=f=f90
		fi
		if test x"$FLAGS_NO_CHECK" = x; then
		    FLAGS_NO_CHECK="\$(FLAGS_CHECK)"
		fi
                # use this on IDRIS machines, www.idris.fr
                #FLAGS_CHECK = -q64 -O4 -qfree=f90 -qsuffix=f=f90
                # use this on MareNostrum at the Barcelona SuperComputing Center (Spain)
                #FLAGS_CHECK = -O3 -qstrict -qtune=ppc970 -qarch=ppc970 -qfree=f90 -qsuffix=f=f90
                ;;
        esac
        ;;
esac

AC_FC_WRAPPERS
AC_PROG_CC


# influential environment variables

AC_ARG_VAR(MPIFC, [MPI Fortran compiler command])
AC_ARG_VAR(MPILIBS, [extra libraries for linking MPI programs])
AC_ARG_VAR(MPICC, [MPI C compiler command])
AC_ARG_VAR(FLAGS_CHECK, [Fortran compiler flags for non-critical subroutines])
AC_ARG_VAR(FLAGS_NO_CHECK, [Fortran compiler flags for creating fast, production-run code for critical subroutines])
if test x"$MPIFC" = x; then
    MPIFC=mpif90
fi
if test x"$MPICC" = x; then
    MPICC=mpicc
fi

AC_ARG_VAR(LOCAL_PATH_IS_ALSO_GLOBAL, [files on a local path on each node are also seen as global with same path @<:@default=false@:>@])
if test x"$LOCAL_PATH_IS_ALSO_GLOBAL" = x; then
    LOCAL_PATH_IS_ALSO_GLOBAL=false
fi

AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])



# Checks for Python modules and packages.

if test "$want_pyre" = yes; then
    PYTHONPATH="${pythondir}${PYTHONPATH:+:${PYTHONPATH}}"; export PYTHONPATH
    CIT_PATH_PYTHIA([0.8], [mpi])
    AC_SUBST([PYTHONPATH])
fi


# Checks for libraries.


# Checks for header files.


# Checks for typedefs, structures, and compiler characteristics.

if test "$want_pyre" = yes; then
    CIT_FC_MAIN
fi


# Checks for library functions.


AC_CONFIG_FILES([Makefile constants.h precision.h])
AC_OUTPUT


dnl end of configure.ac
