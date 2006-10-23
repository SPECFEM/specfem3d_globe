# CIT_PYTHON_INCDIR
# -----------------
# Determine the directory containing <Python.h> using distutils.
AC_DEFUN([CIT_PYTHON_INCDIR], [
# $Id: cit_python_incdir.m4 2367 2005-09-09 16:46:52Z leif $
AC_REQUIRE([AM_PATH_PYTHON])
AC_CACHE_CHECK([for $am_display_PYTHON include directory],
    [PYTHON_INCDIR],
    [PYTHON_INCDIR=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_inc()" 2>/dev/null ||
     echo "$PYTHON_PREFIX/include/python$PYTHON_VERSION"`])
AC_SUBST([PYTHON_INCDIR], [$PYTHON_INCDIR])
])dnl CIT_PYTHON_INCDIR
dnl end of file
