# CIT_PROG_PYCONFIG
# -----------------
# Provide a simple Python script which generates a Python module to
# expose our package configuration, similar to Python's
# distutils.sysconfig.
AC_DEFUN([CIT_PROG_PYCONFIG], [
# $Id: cit_prog_pyconfig.m4 2659 2006-04-01 01:41:01Z leif $
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
dnl end of file
