# -*- Autoconf -*-


## -------------------------------- ##
## Autoconf macros for Python eggs. ##
## -------------------------------- ##


# CIT_CHECK_PYTHON_EGG(REQUIREMENT,
#                      [ACTION-IF-FOUND, [ACTION-IF-NOT-FOUND]])
# --------------------------------------------------------------

# Check for REQUIREMENT using pkg_resources.require().  If the
# corresponding distribution is found, append it to the list of
# requirements and execute ACTION-IF-FOUND.  Otherwise, execute
# ACTION-IF-NOT-FOUND.

AC_DEFUN([CIT_CHECK_PYTHON_EGG], [
# $Id: cit_python_eggs.m4 4616 2006-09-25 23:41:07Z leif $

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
    cit_egg_requirements="$1:$cit_egg_requirements"
    $2
else
    AC_MSG_RESULT(no)
    m4_default([$3], [AC_MSG_ERROR([required Python package not found; try running "$PYTHON setup.py"])])
fi

])dnl CIT_CHECK_PYTHON_EGG


# CIT_PYTHON_EGG_FLAGS
# --------------------

# Perform a breadth-first traversal of Python dependencies (as
# indicated by the requirements accumulated by CIT_CHECK_PYTHON_EGG).
# Set PYTHON_EGG_CFLAGS, PYTHON_EGG_CPPFLAGS, and PYTHON_EGG_LDFLAGS
# according to each dependency's "config.cfg" metadata, if present.

# Loosely inspired by PKG_CHECK_MODULES.  See pkg-config(1).

AC_DEFUN([CIT_PYTHON_EGG_FLAGS], [
# $Id: cit_python_eggs.m4 4616 2006-09-25 23:41:07Z leif $

AC_MSG_CHECKING([for egg-related flags])

cat >check_python_egg.py <<END_OF_PYTHON
[
try:
    from pkg_resources import require
except Exception, e:
    print >>sys.stderr, e
    sys.exit(0)

import sys
from ConfigParser import ConfigParser, NoOptionError
from StringIO import StringIO

flags = dict(
    CFLAGS = [],
    CPPFLAGS = [],
    LDFLAGS = [],
)

cit_egg_requirements = "$cit_egg_requirements"
requirements = cit_egg_requirements.split(':')

deps = require(*requirements)
deps.reverse()
dependencies = []
processed = {}
for dist in deps:
    if dist in processed:
        continue
    dependencies.insert(0, dist)
    processed[dist] = True
for dist in dependencies:
    if dist.has_metadata('config.cfg'):
        parser = ConfigParser({'location': dist.location})
        config = dist.get_metadata('config.cfg')
        fp = StringIO(config)
        parser.readfp(fp, 'config.cfg')
        for k,v in flags.iteritems():
            try:
                v.append(parser.get('flags', k))
            except NoOptionError:
                pass

for k,v in flags.iteritems():
    print 'PYTHON_EGG_%s="%s"' % (k, ' '.join(v))
]
END_OF_PYTHON

AS_IF([AC_TRY_COMMAND([$PYTHON check_python_egg.py >conftest.sh 2>&AS_MESSAGE_LOG_FD])],
      [AC_MSG_RESULT(ok)],
      [AC_MSG_RESULT(failed)
      AC_MSG_FAILURE([cannot scan Python eggs for flags])])
eval `cat conftest.sh`
rm -f conftest.sh check_python_egg.py

AC_SUBST(PYTHON_EGG_CFLAGS)
AC_SUBST(PYTHON_EGG_CPPFLAGS)
AC_SUBST(PYTHON_EGG_LDFLAGS)

])dnl CIT_PYTHON_EGG_FLAGS


# CIT_PYTHON_EGG_REQUIRES
# -----------------------

# Dump Python egg requirements (accumulated by CIT_CHECK_PYTHON_EGG)
# to 'requires.txt'.

AC_DEFUN([CIT_PYTHON_EGG_REQUIRES], [
# $Id: cit_python_eggs.m4 4616 2006-09-25 23:41:07Z leif $

ofile=requires.txt
requiresfile="${ofile}T"
trap "rm \"$requiresfile\"; exit 1" 1 2 15
rm -f "$requiresfile"

AC_MSG_NOTICE([creating $ofile])

cit_save_IFS=$IFS; IFS=:
for cit_egg_requirement in $cit_egg_requirements
do
    IFS=$cit_save_IFS
    echo $cit_egg_requirement >>$requiresfile
done

mv -f "$requiresfile" "$ofile" || \
    (rm -f "$ofile" && cp "$requiresfile" "$ofile" && rm -f "$requiresfile")

AC_SUBST([pythoneggdir], [\${pythondir}/$PACKAGE-$PACKAGE_VERSION.egg])
AC_SUBST([pythonegginfodir], [\${pythoneggdir}/EGG-INFO])

])dnl CIT_PYTHON_EGG_REQUIRES


dnl end of file
