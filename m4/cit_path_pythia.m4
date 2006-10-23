# CIT_PATH_PYTHIA([VERSION], [SUBPACKAGES],
#                 [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------
# Check for the Pythia package.  If SUBPACKAGES is
# specified, check for each whitespace-separated subpackage
# listed (useful for optional subpackages such as 'mpi'
# and 'acis').
AC_DEFUN([CIT_PATH_PYTHIA], [
# $Id: cit_path_pythia.m4 2721 2006-04-11 02:45:29Z leif $
AC_REQUIRE([AM_PATH_PYTHON])
# undocumented configure arg --with-pythia=[auto|prepackaged|VERSION|VERSION-SUBPACKAGES]
if test "${with_pythia+set}" = set; then
    case "$with_pythia" in
        yes | no) want_pythia="auto" ;;
        auto | prepackaged | *.*) want_pythia="$with_pythia" ;;
        * ) want_pythia="auto" ;;
    esac
else
    want_pythia="auto"
fi
pythia_mpi="no"
for pythia_subpackage in : $2; do
    test "x$pythia_subpackage" = x: && continue
    if test "$pythia_subpackage" = "mpi"; then
        pythia_mpi="yes"
    fi
done
AC_MSG_CHECKING([for Pythia v$1])
if test "$want_pythia" = "prepackaged"; then
    if test -d $srcdir/pythia-$1; then
        AC_MSG_RESULT([(prepackaged) yes])
        MAYBE_PYTHIA=pythia-$1
        # Override these tests in any subpackages.
        if test -n "$2"; then
            pythia_version=$1-`echo $2 | sed 's/ /-/g'`
        else
            pythia_version="$1"
        fi
        ac_configure_args="$ac_configure_args --with-pythia=$pythia_version"
        # Find Pythia headers and libraries in the build directory.
        pythia_builddir=`pwd`/pythia-$1
        pythia_pkgdir=$pythia_builddir/packages
        CPPFLAGS="-I$pythia_builddir/include $CPPFLAGS"; export CPPFLAGS
        LDFLAGS="-L$pythia_pkgdir/journal/libjournal -L$pythia_pkgdir/mpi $LDFLAGS"; export LDFLAGS
        if test "$pythia_mpi" = "yes"; then
            AC_SUBST([PYTHIA_MPIPYTHON], ["\${bindir}/mpipython.exe"])
        fi
        $3
    else
        AC_MSG_RESULT(no)
        m4_default([$4], [AC_MSG_ERROR([prepackaged Pythia not found])])
        :
    fi
elif test "$want_pythia" != "auto"; then
    # Override the tests.
    pythia_version=`echo $want_pythia | sed 's/-/ /' | sed 's/ .*//'`
    pythia_subpackages=,`echo $want_pythia | sed 's/-/ /' | sed 's/^.* //' | sed 's/-/,/g'`,
    if test "$pythia_version" = $1; then
        AC_MSG_RESULT([(prepackaged) yes])
        pythia_found="yes"
        for pythia_subpackage in : $2; do
           test "x$pythia_subpackage" = x: && continue
            AC_MSG_CHECKING([for subpackage '$pythia_subpackage' in Pythia])
            if test `echo $pythia_subpackages | grep ,$pythia_subpackage,`; then
                AC_MSG_RESULT([(prepackaged) yes])
            else
                AC_MSG_RESULT([(prepackaged) no])
                pythia_found="no"
            fi
        done
        if test "$pythia_found" = "yes"; then
            if test "$pythia_mpi" = "yes"; then
                AC_SUBST([PYTHIA_MPIPYTHON], ["\${bindir}/mpipython.exe"])
            fi
            $3
        else
            m4_default([$4], [AC_MSG_ERROR([prepackaged Pythia is unsuitable; need subpackages: $2])])
            :
        fi
    else
        AC_MSG_RESULT([(prepackaged) no])
        m4_default([$4], [AC_MSG_ERROR([prepackaged Pythia v$pythia_version is unsuitable; need v$1])])
    fi
else
    # It is common practice to create a 'pyre' project subdirectory, which
    # Python will search instead of the installed Pyre!
    test -d empty || mkdir empty
    pythia_version=`cd empty && $PYTHON -c "import pyre; print pyre.__version__" 2>/dev/null`
    rmdir empty
    if test "$pythia_version" = $1; then
        AC_MSG_RESULT(yes)
        pythia_found="yes"
        for pythia_subpackage in : $2; do
            test "x$pythia_subpackage" = x: && continue
            AC_MSG_CHECKING([for subpackage '$pythia_subpackage' in Pythia])
            test -d empty || mkdir empty
            pythia_subversion=`cd empty && $PYTHON -c "import $pythia_subpackage; print $pythia_subpackage.__version__" 2>/dev/null`
            rmdir empty
            if test "$pythia_subversion" = $1; then
                AC_MSG_RESULT(yes)
            else
                AC_MSG_RESULT(no)
                pythia_found="no"
            fi
        done
        if test "$pythia_found" = "yes"; then
            AC_MSG_CHECKING([Pythia include directory])
            test -d empty || mkdir empty
            [pythia_pkgincludedir=`cd empty && $PYTHON -c "from pyre.config import makefile; print makefile['pkgincludedir']" 2>/dev/null`]
            rmdir empty
            if test -d "$pythia_pkgincludedir"; then
                AC_MSG_RESULT([$pythia_pkgincludedir])
                CPPFLAGS="-I$pythia_pkgincludedir $CPPFLAGS"; export CPPFLAGS
            else
                AC_MSG_RESULT(no)
            fi
            AC_MSG_CHECKING([Pythia lib directory])
            test -d empty || mkdir empty
            [pythia_libdir=`cd empty && $PYTHON -c "from pyre.config import makefile; print makefile['libdir']" 2>/dev/null`]
            rmdir empty
            if test -d "$pythia_libdir"; then
                AC_MSG_RESULT([$pythia_libdir])
                LDFLAGS="-L$pythia_libdir $LDFLAGS"; export LDFLAGS
            else
                AC_MSG_RESULT(no)
            fi
            AC_MSG_CHECKING([Pythia bin directory])
            test -d empty || mkdir empty
            [pythia_bindir=`cd empty && $PYTHON -c "from pyre.config import makefile; print makefile['bindir']" 2>/dev/null`]
            rmdir empty
            if test -d "$pythia_bindir"; then
                AC_MSG_RESULT([$pythia_bindir])
            else
                AC_MSG_RESULT(no)
            fi
            AC_CHECK_LIB(journal, firewall_hit, [
                AC_LANG_PUSH(C++)
                AC_CHECK_HEADER([journal/diagnostics.h], [
                    if test "$pythia_mpi" = "no"; then
                        :
                    elif test -n "$pythia_bindir"; then
                        AC_MSG_CHECKING([for mpipython.exe])
                        if test -x "$pythia_bindir/mpipython.exe"; then
                            AC_SUBST([PYTHIA_MPIPYTHON], ["$pythia_bindir/mpipython.exe"])
                            AC_MSG_RESULT([$PYTHIA_MPIPYTHON])
                            $3
                        else
                            AC_MSG_RESULT(no)
                            m4_default([$4], [AC_MSG_ERROR([Pythia program 'mpipython.exe' not found])])
                        fi
                    else
                        AC_PATH_PROG([PYTHIA_MPIPYTHON], [mpipython.exe], [no])
                        if test "$PYTHIA_MPIPYTHON" != "no"; then
                            $3
                            :
                        else
                            m4_default([$4], [AC_MSG_ERROR([Pythia program 'mpipython.exe' not found])])
                            :
                        fi
                    fi
                ], [
                    m4_default([$4], [AC_MSG_ERROR([Pythia headers not found; try CPPFLAGS="-I<pythia-$1 include dir>"])])
                    :
                ])
                AC_LANG_POP(C++)
            ], [
                m4_default([$4], [AC_MSG_ERROR([Pythia libraries not found; try LDFLAGS="-L<Pythia lib dir>"])])
                :
            ])
        else
            m4_default([$4], [AC_MSG_ERROR([required Pythia subpackages not found: $2])])
            :
        fi
    else
        AC_MSG_RESULT(no)
        AC_MSG_CHECKING([for prepackaged Pythia])
        if test -d $srcdir/pythia-$1; then
            AC_MSG_RESULT(yes)
            MAYBE_PYTHIA=pythia-$1
            # Override the above tests in any subpackages.
            if test -n "$2"; then
                pythia_version=$1-`echo $2 | sed 's/ /-/g'`
            else
                pythia_version="$1"
            fi
            ac_configure_args="$ac_configure_args --with-pythia=$pythia_version"
            # Find Pythia headers and libraries in the build directory.
            pythia_builddir=`pwd`/pythia-$1
            pythia_pkgdir=$pythia_builddir/packages
            CPPFLAGS="-I$pythia_builddir/include $CPPFLAGS"; export CPPFLAGS
            LDFLAGS="-L$pythia_pkgdir/journal/libjournal -L$pythia_pkgdir/mpi $LDFLAGS"; export LDFLAGS
            if test "$pythia_mpi" = "yes"; then
                AC_SUBST([PYTHIA_MPIPYTHON], ["\${bindir}/mpipython.exe"])
            fi
            $3
        else
            AC_MSG_RESULT(no)
            m4_default([$4], [AC_MSG_ERROR([no suitable Pythia package found; check PYTHONPATH])])
        fi
    fi
fi
if test -d $srcdir/pythia-$1; then
    MAYBE_DIST_PYTHIA=pythia-$1
fi
AC_SUBST([MAYBE_PYTHIA])
AC_SUBST([MAYBE_DIST_PYTHIA])
])dnl CIT_PATH_PYTHIA
dnl end of file
