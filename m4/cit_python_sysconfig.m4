# CIT_PYTHON_SYSCONFIG
# --------------------
AC_DEFUN([CIT_PYTHON_SYSCONFIG], [
# $Id: cit_python_sysconfig.m4 2367 2005-09-09 16:46:52Z leif $
AC_REQUIRE([AM_PATH_PYTHON])
AC_MSG_CHECKING([$am_display_PYTHON sysconfig])
cat >sysconfig.py <<END_OF_PYTHON
[from distutils import sysconfig
print 'PYTHON_INCDIR="%s"' % sysconfig.get_python_inc()
keys = (
    'BLDLIBRARY',
    'LDFLAGS',
    'LDLAST',
    'LDLIBRARY',
    'LIBDIR',
    'LIBP',
    'LIBPL',
    'LIBS',
    'LINKFORSHARED',
    'MODLIBS',
    'SYSLIBS',
)
vars = sysconfig.get_config_vars()
# transform AIX's python.exp
vars['LINKFORSHARED'] = vars['LINKFORSHARED'].replace('Modules',vars['LIBPL'])
if vars['LDLIBRARY'] == vars['LIBRARY']:
    # "On systems without shared libraries, LDLIBRARY is the same as LIBRARY"
    vars['BLDLIBRARY'] = "-L%(LIBPL)s -lpython%(VERSION)s" % vars
elif vars['BLDLIBRARY']:
    # "On Mac OS X frameworks, BLDLIBRARY is blank"
    vars['BLDLIBRARY'] = "-L%(LIBDIR)s -lpython%(VERSION)s" % vars
for key in keys:
    print 'PYTHON_%s="%s"' % (key, vars.get(key, ''))
]
END_OF_PYTHON
eval `$PYTHON sysconfig.py 2>/dev/null`
if test -n "$PYTHON_INCDIR"; then
    AC_MSG_RESULT(ok)
else
    AC_MSG_ERROR(["failed

Run '$PYTHON sysconfig.py' to see what wrong.
"])
fi
rm -f sysconfig.py
AC_SUBST([PYTHON_INCDIR], [$PYTHON_INCDIR])
AC_SUBST([PYTHON_BLDLIBRARY], [$PYTHON_BLDLIBRARY])
AC_SUBST([PYTHON_LDFLAGS], [$PYTHON_LDFLAGS])
AC_SUBST([PYTHON_LDLAST], [$PYTHON_LDLAST])
AC_SUBST([PYTHON_LDLIBRARY], [$PYTHON_LDLIBRARY])
AC_SUBST([PYTHON_LIBDIR], [$PYTHON_LIBDIR])
AC_SUBST([PYTHON_LIBP], [$PYTHON_LIBP])
AC_SUBST([PYTHON_LIBPL], [$PYTHON_LIBPL])
AC_SUBST([PYTHON_LIBS], [$PYTHON_LIBS])
AC_SUBST([PYTHON_LINKFORSHARED], [$PYTHON_LINKFORSHARED])
AC_SUBST([PYTHON_MODLIBS], [$PYTHON_MODLIBS])
AC_SUBST([PYTHON_SYSLIBS], [$PYTHON_SYSLIBS])
])dnl CIT_PYTHON_SYSCONFIG
dnl end of file
