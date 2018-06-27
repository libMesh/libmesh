# -------------------------------------------------------------
# Hdf5 library
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_HDF5],
[
  AC_ARG_ENABLE(hdf5,
                AS_HELP_STRING([--enable-hdf5],
                               [build libmesh with HDF5 support, the selected HDF5 must be compatible with contrib/netcdf]),
                [AS_CASE("${enableval}",
                         [yes], [enablehdf5=yes],
                         [no],  [enablehdf5=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-hdf5)])],
                [enablehdf5=no])

  AS_IF([test "x$enablehdf5" = "xyes"],
        [
          AX_PATH_HDF5(1.8.0,no)
          AS_IF([test "x$HAVE_HDF5" = "x0"],
                [
                  enablehdf5=no
                  AC_MSG_RESULT(<<< HDF5 support not found or disabled >>>)
                ],
                [AC_MSG_RESULT(<<< Configuring library with HDF5 support >>>)])
        ])
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
  [AS_HELP_STRING([--with-hdf5=DIR],[root directory of HDF5 installation (default = HDF5_DIR)])],
  dnl action-if-given
  [
    with_hdf5=$withval
    AS_IF([test "x${with_hdf5}" != "xyes"], [HDF5_PREFIX=$withval])
  ],
  dnl action-if-not-given
  [
    dnl This is "no" if the user did not specify --with-hdf5=foo
    with_hdf5=$withval

    dnl If $HDF5_DIR is set in the user's environment, then treat that
    dnl as though they had said --with-hdf5=$HDF5_DIR.
    AS_IF([test "x${HDF5_DIR}" != "x"],
          [
            HDF5_PREFIX=${HDF5_DIR}
            with_hdf5=yes
          ])
  ])

dnl package requirement; if not specified, the default is to assume that
dnl the package is optional

dnl GNU-m4 ifelse documentation:
dnl ifelse (string-1, string-2, equal, [not-equal])
dnl If string-1 and string-2 are equal (character for character),
dnl expands to the string in 'equal', otherwise to the string in
dnl 'not-equal'.
is_package_required=ifelse([$2], ,no, $2)

AS_IF([test "x${with_hdf5}" != "xno"],
      [
    AS_IF([test -d "${HDF5_PREFIX}/lib"],
          [
            HDF5_LIBS="-L${HDF5_PREFIX}/lib -lhdf5"
            HDF5_FLIBS="-L${HDF5_PREFIX}/lib -lhdf5_fortran"
            HDF5_CXXLIBS="-L${HDF5_PREFIX}/lib -lhdf5_cpp"
          ])

    dnl If there is an "rpath" flag detected, append it to the various
    dnl LIBS vars.  This avoids hard-coding -Wl,-rpath, in case that is
    dnl not the right approach for some compilers.
    AS_IF([test "x$RPATHFLAG" != "x" && test -d "${HDF5_PREFIX}/lib"],
          [
            HDF5_LIBS="${HDF5_LIBS} ${RPATHFLAG}${HDF5_PREFIX}/lib"
            HDF5_FLIBS="${HDF5_FLIBS} ${RPATHFLAG}${HDF5_PREFIX}/lib"
            HDF5_CXXLIBS="${HDF5_CXXLIBS} ${RPATHFLAG}${HDF5_PREFIX}/lib"
          ])

    AS_IF([test -d "${HDF5_PREFIX}/include"], [HDF5_CPPFLAGS="-I${HDF5_PREFIX}/include"])

    ac_HDF5_save_CFLAGS="$CFLAGS"
    ac_HDF5_save_CPPFLAGS="$CPPFLAGS"
    ac_HDF5_save_LDFLAGS="$LDFLAGS"
    ac_HDF5_save_LIBS="$LIBS"

    CFLAGS="${HDF5_CPPFLAGS} ${CFLAGS}"
    CPPFLAGS="${HDF5_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${HDF5_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([hdf5.h],[found_header=yes],[found_header=no])

    dnl ----------------------
    dnl  Minimum version check
    dnl ----------------------

    min_hdf5_version=ifelse([$1], ,1.8.0, $1)

    AC_MSG_CHECKING([for HDF5 version >= $min_hdf5_version])

    dnl Strip the major.minor.micro version numbers out of the min version string
    MAJOR_VER=`echo $min_hdf5_version | sed 's/^\([[0-9]]*\).*/\1/'`
    AS_IF([test "x${MAJOR_VER}" = "x"], [MAJOR_VER=0])

    MINOR_VER=`echo $min_hdf5_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    AS_IF([test "x${MINOR_VER}" = "x"], [MINOR_VER=0])

    MICRO_VER=`echo $min_hdf5_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    AS_IF([test "x${MICRO_VER}" = "x"], [MICRO_VER=0])

    dnl begin additional test(s) if header if available
    succeeded=no
    AC_LANG_PUSH([C])

    AS_IF([test "x${found_header}" = "xyes"],
          [
            min_version_succeeded=no
            hdf5_has_cxx=no

            dnl Test that HDF5 version is greater than or equal to the required min version.
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
              @%:@include <hdf5.h>
                  ]], [[
                  @%:@if H5_VERS_MAJOR > $MAJOR_VER
                  /* Sweet nibblets */
                  @%:@elif (H5_VERS_MAJOR >= $MAJOR_VER) && (H5_VERS_MINOR >= $MINOR_VER) && (H5_VERS_RELEASE >= $MICRO_VER)
                  /* Winner winner, chicken dinner */
                  @%:@else
                  @%:@  error HDF5 version is too old
                  @%:@endif
              ]])],[
                  min_version_succeeded=yes
              ],[
                  min_version_succeeded=no
              ])

            AC_LANG_POP([C])

            AS_IF([test "x$min_version_succeeded" = "xno"],
                  [
                    AC_MSG_RESULT(no)
                    AS_IF([test "x$is_package_required" = "xyes"],
                          [AC_MSG_ERROR([Your HDF5 library version does not meet the minimum version requirement (HDF5 >= $min_hdf5_version). Please use --with-hdf5 to specify the location of a valid installation.])])
                  ],
                  [AC_MSG_RESULT(yes)])

            dnl Check for -lhdf5
            AC_CHECK_LIB([hdf5],[H5Fopen],[found_library=yes],[found_library=no])

            dnl Test for the HDF5 C++ interface by trying to link a test code.
            AC_LANG_PUSH([C++])

            AC_MSG_CHECKING([If HDF5 C++ interface is present])

            dnl Using the C++ interface requires linking against both the C
            dnl and C++ libs.
            LIBS="${HDF5_LIBS} ${HDF5_CXXLIBS}"

            AC_LINK_IFELSE([AC_LANG_PROGRAM([[
              @%:@include <H5Cpp.h>
              @%:@ifndef H5_NO_NAMESPACE
              using namespace H5;
              @%:@endif
                  ]], [[
              H5std_string  fname("test.h5");
              H5File file (fname, H5F_ACC_TRUNC);
              ]])],[
                  hdf5_has_cxx=yes
              ],[
                  hdf5_has_cxx=no
              ])

            AC_LANG_POP([C++])

            dnl Not having the C++ interface doesn't disqualify us from using
            dnl the C interface.  We'll set a define if C++ is available, so
            dnl code can conditionally make use of it.
            AS_IF([test "x$hdf5_has_cxx" = "xyes"],
                  [
                    AC_MSG_RESULT(yes)
                    AC_DEFINE(HAVE_HDF5_CXX, 1, [Define if the HDF5 C++ interface is available])
                  ],
                  [AC_MSG_RESULT(no)])

            succeeded=no
            AS_IF([test "x$found_header" = "xyes" && test "x$min_version_succeeded" = "xyes" && test "x$found_library" = "xyes"],
                  [succeeded=yes])
          ])

    dnl Reset variables used by configure tests.
    CFLAGS="$ac_HDF5_save_CFLAGS"
    CPPFLAGS="$ac_HDF5_save_CPPFLAGS"
    LDFLAGS="$ac_HDF5_save_LDFLAGS"
    LIBS="$ac_HDF5_save_LIBS"

    AS_IF([test "x$succeeded" = "xno"],
          [
            AS_IF([test "x$is_package_required" = "xyes"],
                  [AC_MSG_ERROR([HDF5 not found.  Try either --with-hdf5 or setting HDF5_DIR.])],
                  [AC_MSG_NOTICE([optional HDF5 library not found, or does not meet version requirements])])

            HDF5_CFLAGS=""
            HDF5_CPPFLAGS=""
            HDF5_LIBS=""
            HDF5_FLIBS=""
            HDF5_CXXLIBS=""
            HDF5_PREFIX=""
          ],
          [
            HAVE_HDF5=1
            AC_DEFINE(HAVE_HDF5,1,[Define if HDF5 is available])
            AC_SUBST(HDF5_CFLAGS)
            AC_SUBST(HDF5_CPPFLAGS)
            AC_SUBST(HDF5_LIBS)
            AC_SUBST(HDF5_FLIBS)
            AC_SUBST(HDF5_CXXLIBS)
            AC_SUBST(HDF5_PREFIX)
          ])
      ])

#AM_CONDITIONAL(HDF5_ENABLED,test x$HAVE_HDF5 = x1)

])
