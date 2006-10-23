# _CIT_FC_MAIN
# ------------
# Define {F77,FC}_MAIN to the name of the alternate main() function
# for use with the Fortran libraries (i.e., MAIN__ or whatever), or
# 'main' if no such alternate name is found.
#
# As of Autoconf 2.59, the macro AC_FC_MAIN does not work with ifort
# v9, because the macro assumes that 'main' will be resolved by
# FCLIBS, but FCLIBS does not include Intel's 'for_main.o'.  This
# macro simply links with the Fortran compiler instead.
#
AC_DEFUN([_CIT_FC_MAIN],
[_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([for alternate main to link with Fortran libraries],
               ac_cv_[]_AC_LANG_ABBREV[]_main,
[ac_[]_AC_LANG_ABBREV[]_m_save_LIBS=$LIBS
 LIBS="cfortran_test.$ac_objext $LIBS"
 ac_fortran_dm_var=[]_AC_FC[]_DUMMY_MAIN
 ac_cv_fortran_main="main" # default entry point name
 for ac_func in MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ main__ _main; do
   AC_LANG_PUSH(C)
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([@%:@ifdef FC_DUMMY_MAIN_EQ_F77
@%:@  undef F77_DUMMY_MAIN
@%:@  undef FC_DUMMY_MAIN
@%:@else
@%:@  undef $ac_fortran_dm_var
@%:@endif
@%:@define main $ac_func])],
                  [mv conftest.$ac_objext cfortran_test.$ac_objext],
                  [AC_MSG_FAILURE([cannot compile a simple C program])])
   AC_LANG_POP(C)
   AC_LINK_IFELSE([AC_LANG_SOURCE(
[      subroutine foobar()
      return
      end])], [ac_cv_fortran_main=$ac_func; break])
   rm -f cfortran_test* conftest*
 done
 ac_cv_[]_AC_LANG_ABBREV[]_main=$ac_cv_fortran_main
 rm -f cfortran_test* conftest*
 LIBS=$ac_[]_AC_LANG_ABBREV[]_m_save_LIBS
])
AC_DEFINE_UNQUOTED([]_AC_FC[]_MAIN, $ac_cv_[]_AC_LANG_ABBREV[]_main,
                   [Define to alternate name for `main' routine that is
                    called from a `main' in the Fortran libraries.])
])# _CIT_FC_MAIN


# CIT_F77_MAIN
# ------------
AC_DEFUN([CIT_F77_MAIN],
[AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])dnl
AC_LANG_PUSH(Fortran 77)dnl
_AC_FC_MAIN
AC_LANG_POP(Fortran 77)dnl
])# CIT_F77_MAIN


# CIT_FC_MAIN
# -----------
AC_DEFUN([CIT_FC_MAIN],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
AC_LANG_PUSH(Fortran)dnl
_CIT_FC_MAIN
AC_LANG_POP(Fortran)dnl
])# CIT_FC_MAIN


dnl end of file
