# SYNOPSIS
#
#   Test for MetaPhysicL
#
#   AX_PATH_METAPHYSICL( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-metaphysicl=DIR option. Searches
#   --with-metaphysicl, $METAPHYSICL_DIR, and the usual places for
#   METAPHYSICL headers and libraries.
#
#   On success, sets METAPHYSICL_CPPFLAGS, METAPHYSICL_LIBS, and
#   #defines HAVE_METAPHYSICL.
#   Also defines automake conditional METAPHYSICL_ENABLED.  Assumes
#   package is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: metaphysicl.m4 37037 2013-02-16 01:03:09Z pbauman $
#
# COPYLEFT
#
#   Copyright (c) 2013 Roy H. Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2013 Paul T. Bauman <pbauman@ices.utexas.edu>
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

AC_DEFUN([AX_PATH_METAPHYSICL],
[

AC_ARG_VAR(METAPHYSICL_DIR,[root directory of MetaPhysicL installation])

AC_ARG_WITH(metaphysicl,
  [AS_HELP_STRING([--with-metaphysicl[=DIR]],[root directory of MetaPhysicL installation (default = METAPHYSICL_DIR)])],
  [with_metaphysicl=$withval
if test "${with_metaphysicl}" != yes; then
    METAPHYSICL_PREFIX=$withval
fi
],[
with_metaphysicl=$withval
if test "x${METAPHYSICL_DIR}" != "x"; then
   METAPHYSICL_PREFIX=${METAPHYSICL_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_METAPHYSICL=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_metaphysicl}" != no ; then

    if test -d "${METAPHYSICL_PREFIX}/lib" ; then
       METAPHYSICL_LDFLAGS="-L${METAPHYSICL_PREFIX}/lib -Wl,-rpath,${METAPHYSICL_PREFIX}/lib"
       METAPHYSICL_LIBS="-lmetaphysicl"
    fi

    if test -d "${METAPHYSICL_PREFIX}/include" ; then
       METAPHYSICL_CPPFLAGS="-I${METAPHYSICL_PREFIX}/include -I${METAPHYSICL_PREFIX}/src"
    fi

    ac_METAPHYSICL_save_CPPFLAGS="$CPPFLAGS"
    ac_METAPHYSICL_save_LDFLAGS="$LDFLAGS"
    ac_METAPHYSICL_save_LIBS="$LIBS"

    CPPFLAGS="${METAPHYSICL_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${METAPHYSICL_LDFLAGS} ${LDFLAGS}"
    LIBS="${METAPHYSICL_LIBS} ${LDFLAGS}"

    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([metaphysicl/metaphysicl_version.h],[found_header=yes],[found_header=no])
    AC_LANG_POP([C++])

    #-----------------------
    # Minimum version check
    #----------------------

    min_metaphysicl_version=ifelse([$1], ,0.0.0, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_metaphysicl_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_metaphysicl_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_metaphysicl_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for metaphysicl - version >= $min_metaphysicl_version)
        version_succeeded=no

	AC_LANG_PUSH([C++])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include "metaphysicl/metaphysicl_version.h"
            ]], [[
            #if METAPHYSICL_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (METAPHYSICL_MAJOR_VERSION >= $MAJOR_VER) && (METAPHYSICL_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (METAPHYSICL_MAJOR_VERSION >= $MAJOR_VER) && (METAPHYSICL_MINOR_VERSION >= $MINOR_VER) && (METAPHYSICL_MICRO_VERSION >= $MICRO_VER)
            /* I feel like chicken tonight, like chicken tonight? */
            #else
            #  error version is too old
            #endif
        ]])],[
            AC_MSG_RESULT(yes)
            version_succeeded=yes
        ],[
            AC_MSG_RESULT(no)
        ])
	AC_LANG_POP([C++])

    if test "$version_succeeded" != "yes";then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([

   Your METAPHYSICL version does not meet the minimum versioning
   requirements ($min_metaphysicl_version).  Please use --with-metaphysicl to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Library availability

    AC_MSG_CHECKING([for -lmetaphysicl linkage])

    AC_LANG_PUSH([C++])

    AC_LINK_IFELSE(
                  [AC_LANG_PROGRAM([#include "metaphysicl/metaphysicl_version.h"],
                                   [MetaPhysicL::get_metaphysicl_version()])],
                  [AC_MSG_RESULT(yes)
                   found_library=yes],
                  [AC_MSG_RESULT(no) 
                   found_library=no])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CPPFLAGS="$ac_METAPHYSICL_save_CPPFLAGS"
    LDFLAGS="$ac_METAPHYSICL_save_LDFLAGS"
    LIBS="$ac_METAPHYSICL_save_LIBS"

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
          AC_MSG_ERROR([MetaPhysicL not found.  Try either --with-metaphysicl or setting METAPHYSICL_DIR.])
       else
          AC_MSG_NOTICE([optional MetaPhysicL library not found])
          METAPHYSICL_CPPFLAGS=""   # METAPHYSICL_CFLAGS empty on failure
          METAPHYSICL_LDFLAGS=""    # METAPHYSICL_LDFLAGS empty on failure
          METAPHYSICL_LIBS=""       # METAPHYSICL_LIBS empty on failure
       fi
    else
        HAVE_METAPHYSICL=1
        AC_DEFINE(HAVE_METAPHYSICL,1,[Define if MetaPhysicL is available])
        AC_SUBST(METAPHYSICL_CPPFLAGS)
        AC_SUBST(METAPHYSICL_LDFLAGS)
        AC_SUBST(METAPHYSICL_LIBS)
        AC_SUBST(METAPHYSICL_PREFIX)
    fi

    AC_SUBST(HAVE_METAPHYSICL)

# fi

AM_CONDITIONAL(METAPHYSICL_ENABLED,test x$HAVE_METAPHYSICL = x1)

])
