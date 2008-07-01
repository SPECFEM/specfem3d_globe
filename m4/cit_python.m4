# -*- Autoconf -*-


## --------------------------- ##
## Autoconf macros for Python. ##
## --------------------------- ##


# CIT_PYTHON_INCDIR
# -----------------
# Determine the directory containing <Python.h> using distutils.
AC_DEFUN([CIT_PYTHON_INCDIR], [
# $Id$
AC_REQUIRE([AM_PATH_PYTHON])
AC_CACHE_CHECK([for $am_display_PYTHON include directory],
    [PYTHON_INCDIR],
    [PYTHON_INCDIR=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_inc()" 2>/dev/null ||
     echo "$PYTHON_PREFIX/include/python$PYTHON_VERSION"`])
AC_SUBST([PYTHON_INCDIR], [$PYTHON_INCDIR])
])dnl CIT_PYTHON_INCDIR


# CIT_PYTHON_SYSCONFIG
# --------------------
AC_DEFUN([CIT_PYTHON_SYSCONFIG], [
# $Id$
AC_REQUIRE([AM_PATH_PYTHON])
AC_MSG_CHECKING([$am_display_PYTHON sysconfig])
cat >sysconfig.py <<END_OF_PYTHON
[import os, sys
from distutils import sysconfig
def cygpath(wpath):
    s = os.popen('cygpath -u "%s"' % wpath)
    path = s.read().strip()
    s.close()
    return path
incdir = sysconfig.get_python_inc()
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
    'LA_LDFLAGS',
)
if os.name == "nt":
    # We are running under Python for Windows (the real one...
    # not Cygwin Python, under which 'os.name' is 'posix').
    # We assume that we are still in the Cygwin POSIX environment,
    # however (this is 'configure', after all); so we convert
    # all Windows pathnames to POSIX pathnames using 'cygpath'.
    incdir = cygpath(incdir)
    vars = {}
    libs = os.path.join(sys.prefix, "libs")
    libs = cygpath(libs)
    version = sysconfig.get_python_version()
    version = version.replace('.', '')
    vars['BLDLIBRARY'] = "-L%s -lpython%s" % (libs, version)
else:
    vars = sysconfig.get_config_vars()
    # transform AIX's python.exp
    vars['LINKFORSHARED'] = vars['LINKFORSHARED'].replace('Modules',vars['LIBPL'])
    if vars['LDLIBRARY'] == vars['LIBRARY']:
        # "On systems without shared libraries, LDLIBRARY is the same as LIBRARY"
        vars['BLDLIBRARY'] = "-L%(LIBPL)s -lpython%(VERSION)s" % vars
    elif vars['BLDLIBRARY']:
        #     The 'mpicc' wrapper for LAM/MPI isn't very smart about "-L"
        # options.  Adding "-L/usr/lib" can cause "-lmpi" to be found in /usr/lib
        # instead of LAM's 'lib' directory.  Of course, avoiding "-L/usr/lib"
        # doesn't really fix the problem, but it does make it much less likely;
        # and "-L/usr/lib" is redundant and potentially problematic anyway.
        #     Python 2.4 and later puts a symlink to libpython.so in LIBPL
        # (/usr/lib/python2.x/config), which makes adding "-L$LIBDIR"
        # (in addition to "-L$LIBPL") completely redundant.
        #     But we still support Python 2.3, and we prefer shared to static,
        # so we still add "-L$LIBDIR" when Python is installed in a non-standard
        # location.  Note that the linker will still prefer shared over static
        # with only "-L/usr/lib/python2.3/config" on the link line.
        libdir = ""
        if vars['LIBDIR'] != "/usr/lib":
            libdir = "-L%(LIBDIR)s "
        # Important: on Cygwin, the import library for libpython.dll is
        # nested inside Python's 'config' directory (see Issue39).  This means
        # that the linker always needs help finding "-lpython2.x" (in the form
        # of "-L$LIBPL"), even for the "system" Python installed under /usr.
        vars['BLDLIBRARY'] = (libdir + "-L%(LIBPL)s -lpython%(VERSION)s") % vars
    else:
        # "On Mac OS X frameworks, BLDLIBRARY is blank"
        # See also Issue39.
        framework = "%(PYTHONFRAMEWORKDIR)s/Versions/%(VERSION)s/%(PYTHONFRAMEWORK)s" % vars
        PYTHONFRAMEWORK = vars.get('PYTHONFRAMEWORK', 'Python')
        vars['LINKFORSHARED'] = vars['LINKFORSHARED'].replace(framework, "-framework " + PYTHONFRAMEWORK)
        vars['LA_LDFLAGS'] = "-Wl,-framework,%s" % PYTHONFRAMEWORK
print 'PYTHON_INCDIR="%s"' % incdir
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
AC_SUBST([PYTHON_LA_LDFLAGS], [$PYTHON_LA_LDFLAGS])
])dnl CIT_PYTHON_SYSCONFIG


# CIT_PYTHON_SITE
# ---------------
AC_DEFUN([CIT_PYTHON_SITE], [
# $Id$
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


# CIT_CHECK_PYTHON_EGG(REQUIREMENT,
#                      [ACTION-IF-FOUND, [ACTION-IF-NOT-FOUND]])
# --------------------------------------------------------------

# Check for REQUIREMENT using pkg_resources.require().  If the
# corresponding distribution is found, execute ACTION-IF-FOUND.
# Otherwise, execute ACTION-IF-NOT-FOUND.

AC_DEFUN([CIT_CHECK_PYTHON_EGG], [
# $Id$

AC_MSG_CHECKING([for "$1"])

cat >check_python_egg.py <<END_OF_PYTHON
[
import sys
try:
    from pkg_resources import require
    require("$1")
except Exception, e:
    print >>sys.stderr, e
    print "cit_egg_status=1"
else:
    print "cit_egg_status=0"
]
END_OF_PYTHON

AS_IF([AC_TRY_COMMAND([$PYTHON check_python_egg.py >conftest.sh 2>&AS_MESSAGE_LOG_FD])],
      [],
      [AC_MSG_RESULT(failed)
      AC_MSG_FAILURE([cannot check for Python eggs])])
eval `cat conftest.sh`
rm -f conftest.sh check_python_egg.py

if test "$cit_egg_status" == 0; then
    AC_MSG_RESULT(yes)
    $2
else
    AC_MSG_RESULT(no)
    m4_default([$3], [AC_MSG_ERROR([required Python package not found: $1])])
fi

])dnl CIT_CHECK_PYTHON_EGG


# CIT_PYTHON_EGG_SETUP
# --------------------

AC_DEFUN([CIT_PYTHON_EGG_SETUP], [
# $Id$
AC_REQUIRE([AM_PATH_PYTHON])

cit_builddir=`pwd`
cit_save_PYTHONPATH="$PYTHONPATH"
PYTHONPATH="$cit_builddir/python:$PYTHONPATH"; export PYTHONPATH
cd $srcdir

AC_MSG_NOTICE([downloading missing Python dependencies])
AS_IF([AC_TRY_COMMAND([$PYTHON setup.py install_deps -f $cit_builddir/deps -zmxd $cit_builddir/deps >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD])],
      [],
      [AC_MSG_FAILURE([cannot download missing Python dependencies])])

AC_MSG_NOTICE([building Python dependencies])
AS_IF([AC_TRY_COMMAND([$PYTHON setup.py develop -H None -f $cit_builddir/deps -x -d $cit_builddir/python >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD])],
      [],
      [AC_MSG_FAILURE([building Python dependencies])])

AC_MSG_CHECKING([for egg-related flags])
AS_IF([AC_TRY_COMMAND([$PYTHON setup.py egg_flags >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD])],
      [AC_MSG_RESULT(ok)
       . egg-flags.sh
       rm -f egg-flags.sh
      ],
      [AC_MSG_RESULT(failed)
      AC_MSG_FAILURE([cannot scan Python eggs for flags])])

cd $cit_builddir
PYTHONPATH="$cit_save_PYTHONPATH"
PYTHONPATH="${pythondir}:${pyexecdir}${cit_save_PYTHONPATH:+:${cit_save_PYTHONPATH}}"

AC_SUBST(PYTHONPATH)
AC_SUBST(PYTHON_EGG_CFLAGS)
AC_SUBST(PYTHON_EGG_CPPFLAGS)
AC_SUBST(PYTHON_EGG_LDFLAGS)
AC_SUBST(PYTHON_EGG_LIBS)
AC_SUBST(PYTHON_EGG_PYXFLAGS)

])dnl CIT_PYTHON_EGG_SETUP


# CIT_PROG_PYCONFIG
# -----------------
# Provide a simple Python script which generates a Python module to
# expose our package configuration, similar to Python's
# distutils.sysconfig.
AC_DEFUN([CIT_PROG_PYCONFIG], [
# $Id$
PYCONFIG='$(top_builddir)/pyconfig'
AC_SUBST(PYCONFIG)
ofile=pyconfig
cfgfile="${ofile}T"
trap "rm \"$cfgfile\"; exit 1" 1 2 15
rm -f "$cfgfile"
AC_MSG_NOTICE([creating $ofile])
cat >"$cfgfile" <<END_OF_PYTHON
[#!/usr/bin/env python

from getopt import getopt, GetoptError
from sys import argv, exit
from getopt import getopt
from distutils.sysconfig import parse_config_h, parse_makefile, expand_makefile_vars

def printUsage():
    print "Usage: %s -h HEADER -m MAKEFILE -o OUTPUT" % argv[0]

try:
    (opts, args) = getopt(argv[1:], "h:m:o:")
except GetoptError, error:
    print "%s: %s" % (argv[0], error)
    printUsage()
    exit(1)

header = '';
makefile = '';
output = '';
for option, parameter in opts:
    if option == '-h':
        header = parameter
    elif option == '-m':
        makefile = parameter
    elif option == '-o':
        output = parameter
if not (header and makefile and output):
    printUsage()
    exit(1)

f = open(header)
config_vars = parse_config_h(f)
f.close()

makefile_vars = parse_makefile(makefile)
keys = makefile_vars.keys()
for key in keys:
    makefile_vars[key] = expand_makefile_vars(makefile_vars[key], makefile_vars)

f = open(output, 'w')
print >>f, "#!/usr/bin/env python"
print >>f
print >>f, "config =", config_vars
print >>f
print >>f, "makefile =", makefile_vars
print >>f
print >>f, "# end of file"
f.close()

# end of file]
END_OF_PYTHON
mv -f "$cfgfile" "$ofile" || \
    (rm -f "$ofile" && cp "$cfgfile" "$ofile" && rm -f "$cfgfile")
chmod +x "$ofile"
])dnl CIT_PROG_PYCONFIG


# CIT_PATH_NEMESIS
# -----------------
AC_DEFUN([CIT_PATH_NEMESIS], [
# $Id$
AC_BEFORE([$0], [AM_PATH_PYTHON])
AC_PATH_PROG(PYTHON, nemesis, no)
if test "$PYTHON" = no; then
    AC_MSG_ERROR([program 'nemesis' not found])
fi
])dnl CIT_PATH_NEMESIS


dnl end of file
