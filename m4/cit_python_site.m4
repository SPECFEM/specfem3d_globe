# CIT_PYTHON_SITE
# ---------------
AC_DEFUN([CIT_PYTHON_SITE], [
# $Id: cit_python_site.m4 2659 2006-04-01 01:41:01Z leif $
AC_REQUIRE([AM_PATH_PYTHON])
AC_MSG_CHECKING([whether we are installing to Python's prefix])
cit_python_prefix=`$PYTHON -c "import sys; print sys.prefix"`
if test "$cit_python_prefix" = "$prefix"; then
    AC_MSG_RESULT(yes)
    cit_cond_python_site=true
else
    AC_MSG_RESULT(no)
    cit_cond_python_site=false
fi
AC_MSG_CHECKING([whether we are installing to Python's exec prefix])
cit_python_exec_prefix=`$PYTHON -c "import sys; print sys.exec_prefix"`
cit_exec_prefix=$exec_prefix
test "x$cit_exec_prefix" = xNONE && cit_exec_prefix=$prefix
if test "$cit_python_exec_prefix" = "$cit_exec_prefix"; then
    AC_MSG_RESULT(yes)
    cit_cond_pyexec_site=true
else
    AC_MSG_RESULT(no)
    cit_cond_pyexec_site=false
fi
AM_CONDITIONAL([COND_PYTHON_SITE], [$cit_cond_python_site])
AM_CONDITIONAL([COND_PYEXEC_SITE], [$cit_cond_pyexec_site])
])dnl CIT_PYTHON_SITE
dnl end of file
