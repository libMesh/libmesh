# SYNOPSIS
#
#   Test for MASA Library
#
#   AM_PATH_MASA( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-masa=DIR option. Searches --with-masa,
#   $MASA_DIR, and the usual places for MASA headers and libraries.
#
#   On success, sets MASA_CXXFLAGS, MASA_LIBS, MASA_FC_LIBS (for
#   Fortran) and #defines HAVE_MASA.  Assumes package is optional
#   unless overridden with $2=yes
#
# LAST MODIFICATION
#
#   $Id: masa.m4 38825 2013-04-23 17:53:45Z benkirk $
#
# COPYLEFT
#
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.
#

AC_DEFUN([AX_PATH_MASA],
[

AC_ARG_VAR(MASA_DIR,[root directory of MASA installation])

AC_ARG_WITH(masa,
  [AS_HELP_STRING([--with-masa[=DIR]],[root directory of MASA installation (default = MASA_DIR)])],
  [with_masa=$withval
if test "${with_masa}" != yes; then
    MASA_PREFIX=$withval
fi
],[
# assume a sensible default of --with-masa=yes
with_masa=yes
if test "x${MASA_DIR}" != "x"; then
   MASA_PREFIX=${MASA_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_MASA=0

if test "${with_masa}" != no ; then

    if test -d "${MASA_PREFIX}/lib" ; then
       MASA_LIBS="-L${MASA_PREFIX}/lib -lmasa -Wl,-rpath,${MASA_PREFIX}/lib"
       MASA_FC_LIBS="-L${MASA_PREFIX}/lib -lfmasa -lmasa -Wl,-rpath,${MASA_PREFIX}/lib"
       MASA_FCFLAGS="-I${MASA_PREFIX}/lib -DHAVE_MASA"
    fi

    if test -d "${MASA_PREFIX}/include" ; then
        MASA_CXXFLAGS="-I${MASA_PREFIX}/include"
    fi

    ac_MASA_save_CXXFLAGS="$CXXFLAGS"
    ac_MASA_save_CPPFLAGS="$CPPFLAGS"
    ac_MASA_save_LDFLAGS="$LDFLAGS"
    ac_MASA_save_LIBS="$LIBS"

    CXXFLAGS="${MASA_CXXFLAGS} ${CXXFLAGS}"
    CPPFLAGS="${MASA_CXXFLAGS} ${CPPFLAGS}"
    LDFLAGS="${MASA_LIBS} ${LDFLAGS}"

    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([masa.h],[found_header=yes],[found_header=no])

    #-----------------------
    # Minimum version check
    #----------------------

    min_masa_version=ifelse([$1], ,0.10, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_masa_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_masa_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_masa_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    dnl begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for masa - version >= $min_masa_version)
        version_succeeded=no

    	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
       	@%:@include <masa.h>
            ]], [[
            #if MASA_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (MASA_MAJOR_VERSION >= $MAJOR_VER) && (MASA_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (MASA_MAJOR_VERSION >= $MAJOR_VER) && (MASA_MINOR_VERSION >= $MINOR_VER) && (MASA_MICRO_VERSION >= $MICRO_VER)
            /* Winner winner, chicken dinner */
            #else
            #  error version is too old
            #endif
        ]])],[
            AC_MSG_RESULT(yes)
            version_succeeded=yes
        ],[
            AC_MSG_RESULT(no)
        ])

    if test "$version_succeeded" != "yes";then
       if test "$is_package_required" = yes; then
       	  AC_MSG_ERROR([

   Your MASA library version does not meet the minimum versioning
   requirements ($min_masa_version).  Please use --with-masa to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Library availability

    AC_MSG_CHECKING([for -lmasa linkage])

    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <masa.h>],[MASA::masa_version_stdout])],
    [TEST_LIBS="$TEST_LIBS -lmasa"] [
    AC_MSG_RESULT(yes)
    found_library=yes ],[AC_MSG_RESULT(no)])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CXXFLAGS="$ac_MASA_save_CXXFLAGS"
    CPPFLAGS="$ac_MASA_save_CPPFLAGS"
    LDFLAGS="$ac_MASA_save_LDFLAGS"
    LIBS="$ac_MASA_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$version_succeeded" = yes; then
	    if test "$found_library" = yes; then
               succeeded=yes
	    fi
        fi
    fi

    if test "$succeeded" = no; then
       if test "$is_package_required" = yes; then
       	  AC_MSG_ERROR([MASA not found.  Try either --with-masa or setting MASA_DIR.])
       else
          AC_MSG_NOTICE([optional MASA library not found])
       fi
    else
        HAVE_MASA=1
        AC_DEFINE(HAVE_MASA,1,[Define if MASA is available])
        AC_SUBST(MASA_CXXFLAGS)
	AC_SUBST(MASA_FCFLAGS)
        AC_SUBST(MASA_LIBS)
        AC_SUBST(MASA_FC_LIBS)
	AC_SUBST(MASA_PREFIX)
    fi

    AC_SUBST(HAVE_MASA)
fi

AM_CONDITIONAL(MASA_ENABLED,test x$HAVE_MASA = x1)

])
