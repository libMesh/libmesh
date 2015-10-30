# -------------------------------------------------------------
# Hdf5 library
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_HDF5],
[
  AC_ARG_ENABLE(hdf5,
                AS_HELP_STRING([--disable-hdf5],
                               [build without HDF5 support]),
                [case "${enableval}" in
                  yes)  enablehdf5=yes ;;
                  no)  enablehdf5=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-hdf5) ;;
                esac],
                [enablehdf5=$enableoptional])

  if (test $enablehdf5 = yes); then
    AX_PATH_HDF5(1.8.0,no)
    if (test "x$HAVE_HDF5" = "x0"); then
      enablehdf5=no
    else
       AC_MSG_RESULT(<<< Configuring library with HDF5 support >>>)
    fi
  fi
])




# SYNOPSIS
#
#   Test for HDF5
#
#   AX_PATH_HDF5( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-hdf5=DIR option and minimum version check for
#   the HDF I/O library. Searches --with-hdf5, $HDF5_DIR, and the
#   usual places for HDF5 headers and libraries.
#
#   On success, sets HDF5_CFLAGS, HDF5_LIBS, and #defines HAVE_HDF5.
#   Assumes package is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: hdf5_new.m4 31520 2012-06-28 14:53:30Z mpanesi $
#
# COPYLEFT
#
#   Copyright (c) 2013 Roy H. Stogner <roystgnr@ices.utexas.edu>
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

AC_DEFUN([AX_PATH_HDF5],
[

HAVE_HDF5=0

AC_ARG_VAR(HDF5_DIR,[root directory of HDF5 installation])

AC_ARG_WITH(hdf5,
  [AS_HELP_STRING([--with-hdf5[=DIR]],[root directory of HDF5 installation (default = HDF5_DIR)])],
  [with_hdf5=$withval
if test "${with_hdf5}" != yes; then
    HDF5_PREFIX=$withval
fi
],[
with_hdf5=$withval
if test "x${HDF5_DIR}" != "x"; then
   HDF5_PREFIX=${HDF5_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

if test "${with_hdf5}" != no ; then

    if test -d "${HDF5_PREFIX}/lib" ; then
       HDF5_LIBS="-L${HDF5_PREFIX}/lib -lhdf5 -Wl,-rpath,${HDF5_PREFIX}/lib"
       HDF5_FLIBS="-L${HDF5_PREFIX}/lib -lhdf5_fortran -Wl,-rpath,${HDF5_PREFIX}/lib"
       HDF5_CXXLIBS="-L${HDF5_PREFIX}/lib -lhdf5_cpp -Wl,-rpath,${HDF5_PREFIX}/lib"
    fi

    if test -d "${HDF5_PREFIX}/include" ; then
        HDF5_CPPFLAGS="-I${HDF5_PREFIX}/include"
    fi

    ac_HDF5_save_CFLAGS="$CFLAGS"
    ac_HDF5_save_CPPFLAGS="$CPPFLAGS"
    ac_HDF5_save_LDFLAGS="$LDFLAGS"
    ac_HDF5_save_LIBS="$LIBS"

    CFLAGS="${HDF5_CPPFLAGS} ${CFLAGS}"
    CPPFLAGS="${HDF5_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${HDF5_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([hdf5.h],[found_header=yes],[found_header=no])

    #-----------------------
    # Minimum version check
    #----------------------

    min_hdf5_version=ifelse([$1], ,1.8.0, $1)

    AC_MSG_CHECKING(for hdf5 - version >= $min_hdf5_version)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_hdf5_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_hdf5_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_hdf5_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    succeeded=no
    AC_LANG_PUSH([C])

    if test "x${found_header}" = "xyes" ; then
      version_succeeded=no

      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <hdf5.h>
            ]], [[
            #if H5_VERS_MAJOR > $MAJOR_VER
            /* Sweet nibblets */
            #elif (H5_VERS_MAJOR >= $MAJOR_VER) && (H5_VERS_MINOR >= $MINOR_VER) && (H5_VERS_RELEASE >= $MICRO_VER)
            /* Winner winner, chicken dinner */
            #else
            #  error HDF5 version is too old
            #endif
        ]])],[
            AC_MSG_RESULT(yes)
            version_succeeded=yes
        ],[
            AC_MSG_RESULT(no)
        ])

      AC_LANG_POP([C])

      if test "$version_succeeded" != "yes";then
        if test "$is_package_required" = yes; then
          AC_MSG_ERROR([Your HDF5 library version does not meet the minimum versioning
                        requirements ($min_hdf5_version).  Please use --with-hdf5 to specify the location
                        of an updated installation or consider upgrading the system version.])
        fi
      fi

      # Library availability
      AC_CHECK_LIB([hdf5],[H5Fopen],[found_library=yes],[found_library=no])
      AC_LANG_POP([C])

      succeeded=no
      if test "$found_header" = yes; then
        if test "$version_succeeded" = yes; then
          if test "$found_library" = yes; then
            succeeded=yes
          fi
        fi
      fi
    fi dnl end test if header if available

    CFLAGS="$ac_HDF5_save_CFLAGS"
    CPPFLAGS="$ac_HDF5_save_CPPFLAGS"
    LDFLAGS="$ac_HDF5_save_LDFLAGS"
    LIBS="$ac_HDF5_save_LIBS"

    if test "$succeeded" = no; then
      if test "$is_package_required" = yes; then
        AC_MSG_ERROR([HDF5 not found.  Try either --with-hdf5 or setting HDF5_DIR.])
      else
         AC_MSG_NOTICE([optional HDF5 library not found, or does not meet version requirements])
      fi

      HDF5_CFLAGS=""
      HDF5_CPPFLAGS=""
      HDF5_LIBS=""
      HDF5_FLIBS=""
      HDF5_CXXLIBS=""
      HDF5_PREFIX=""
    else
      HAVE_HDF5=1
      AC_DEFINE(HAVE_HDF5,1,[Define if HDF5 is available])
      AC_SUBST(HDF5_CFLAGS)
      AC_SUBST(HDF5_CPPFLAGS)
      AC_SUBST(HDF5_LIBS)
      AC_SUBST(HDF5_FLIBS)
      AC_SUBST(HDF5_CXXLIBS)
      AC_SUBST(HDF5_PREFIX)
    fi
    #AC_SUBST(HAVE_HDF5)
fi

#AM_CONDITIONAL(HDF5_ENABLED,test x$HAVE_HDF5 = x1)

])
