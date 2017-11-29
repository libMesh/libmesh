# SYNOPSIS
#
#   Add code coverage support with gcov/lcov.
#
#   AX_CODE_COVERAGE()
#
# DESCRIPTION
#
#   Provides a --enable-coverage option which checks for available
#   gcov/lcov binaries and provides ENABLE_CODE_COVERAGE conditional.
#
# LAST MODIFICATION
#
#   $Id: coverage.m4 40881 2013-08-20 17:54:39Z damon $
#
# COPYLEFT
#
#   Copyright (c) 2012 Roy H. Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_CODE_COVERAGE],
[

AC_ARG_ENABLE(coverage, AC_HELP_STRING([--enable-coverage],[configure code coverage analysis tools]))

HAVE_GCOV_TOOLS=0

GCOV_FLAGS=""

if test "x$enable_coverage" = "xyes"; then

   # ----------------------------
   # Check for gcov/lcov binaries
   # ----------------------------
   
   AC_CHECK_PROG(have_gcov,gcov, yes, no)

   if test "x$have_gcov" = "xno"; then
      AC_MSG_ERROR([

      gcov coverage testing tool not found. Please install or update
      your PATH accordingly prior to enabling code coverage features.

      	   ])
   fi

   # ----------------------------------
   # include coverage compiler options
   # ----------------------------------

   HAVE_GCOV_TOOLS=1
   GCOV_FLAGS="-fprofile-arcs -ftest-coverage --coverage"
   GCOV_LDFLAGS="--coverage -lgcov"

   ac_coverage_save_LDFLAGS="$LDFLAGS"
   LDFLAGS="${LDFLAGS} ${GCOV_LDFLAGS}"

   # Test for C...

   ac_coverage_save_CFLAGS="$CFLAGS"

   AC_LANG_PUSH([C])
   CFLAGS="${GCOV_FLAGS} ${CFLAGS}"
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],[], 
     [AC_MSG_ERROR([unable to compile with code coverage ($CC --coverage)])])
   AC_LANG_POP([C])

   # Test for C++...

   ac_coverage_save_CXXFLAGS="$CXXFLAGS"

   AC_LANG_PUSH([C++])
   CXXFLAGS="${GCOV_FLAGS} ${CXXFLAGS}"
   LIBS="-lgcov ${LIBS}"
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],[], 
     [AC_MSG_ERROR([unable to compile with code coverage ($CXX --coverage)])])
   AC_LANG_POP([C++])

   # Test for Fortran...

   ac_coverage_save_FCFLAGS="$FCFLAGS"

   AC_LANG_PUSH([Fortran])
   FCFLAGS="${GCOV_FLAGS} ${FCFLAGS}"
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],[],
     [AC_MSG_ERROR([unable to compile with code coverage ($FC --coverage)])])
   AC_LANG_POP([Fortran])

fi

AC_SUBST(GCOV_FLAGS)
AC_SUBST(GCOV_LDFLAGS)
AC_SUBST(HAVE_GCOV_TOOLS)
AM_CONDITIONAL(CODE_COVERAGE_ENABLED,test x$HAVE_GCOV_TOOLS = x1)

])


