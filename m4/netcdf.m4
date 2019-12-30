dnl -------------------------------------------------------------
dnl netCDF
dnl -------------------------------------------------------------
dnl Declare compilation test function preamble to
dnl be used later with AC_LANG_PROGRAM
m4_define([_AX_CXX_COMPILE_NETCDF_preamble],
          [[
            @%:@include "netcdf.h"
          ]])

dnl Declare compilation test function body to
dnl be used later with AC_LANG_PROGRAM
m4_define([_AX_CXX_COMPILE_NETCDF_body],
          [[
          static int ncid;
          nc_create("test.nc", NC_NETCDF4|NC_CLOBBER, &ncid);
          ]])

dnl Declare compilation test function body to
dnl be used later with AC_LANG_PROGRAM
m4_define([_AX_CXX_COMPILE_NETCDF_body],
          [[
              static int ncid;
              nc_create("test.nc", NC_NETCDF4|NC_CLOBBER, &ncid);
          ]])

AC_DEFUN([CONFIGURE_NETCDF],
[
  dnl User-specific include path
  AC_ARG_WITH(netcdf,
              AS_HELP_STRING([--with-netcdf-include=PATH],[Specify the path for NETCDF header files]),
              withnetcdfinc=$withval,
              withnetcdfinc=no)

  dnl User-specific library path
  AC_ARG_WITH(netcdf-lib,
              AS_HELP_STRING([--with-netcdf-lib=PATH],[Specify the path for NETCDF libs]),
              withnetcdflib=$withval,
              withnetcdflib=no)

  AC_ARG_VAR([NETCDF_INCLUDE], [path to NETCDF header files])
  AC_ARG_VAR([NETCDF_DIR],     [path to NETCDF installation])

  AC_ARG_ENABLE(netcdf,
                AS_HELP_STRING([--disable-netcdf],
                               [build without netCDF binary I/O]),
                [AS_CASE("${enableval}",
                  [yes|new|v4], [enablenetcdf=yes
                                 netcdfmajor=4],
                  [old|v3],     [enablenetcdf=yes
                                 netcdfmajor=3],
                  [no],         [enablenetcdf=no
                                 netcdfmajor=no],
                  [AC_MSG_ERROR(bad value ${enableval} for --enable-netcdf)])],
                [enablenetcdf=$enableoptional; netcdfmajor=4])

  dnl Setting --enable-netcdf-required causes an error to be emitted during
  dnl configure if NETCDF is not successfully detected. This is useful for app
  dnl codes which require NETCDF (like MOOSE-based apps), since it prevents
  dnl situations where libmesh is accidentally built without NETCDF support
  dnl (which may take a very long time), and then the app fails to compile,
  dnl requiring you to redo everything.
  AC_ARG_ENABLE(netcdf-required,
                AC_HELP_STRING([--enable-netcdf-required],
                               [Error if NETCDF is not detected by configure]),
                [AS_CASE("${enableval}",
                         [yes], [netcdfrequired=yes],
                         [no],  [netcdfrequired=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-netcdf-required)])],
                [netcdfrequired=no])

  AS_IF([test "$enablenetcdf" = "yes"],
        [
        AS_IF([test "$withnetcdfinc" != "no"], [NETCDF_INC="$withnetcdfinc"])
        AS_IF([test "$withnetcdflib" != "no"], [NETCDF_LIB="$withnetcdflib"])

        dnl Honor NETCDF_DIR if it is set (or take the provided contrib package)
        AS_IF([test "x$NETCDF_DIR" = "x"], [NETCDF_DIR=${top_srcdir}/contrib/netcdf/v4])

        dnl Look for NETCDF location in the environment, then default paths
        NETCDF_LS_CHECK=$(dirname $(ls -d $NETCDF_DIR/include/netcdf.h 2>/dev/null | tail -n 1) 2>/dev/null)
        AS_IF([test "x$NETCDF_INC" = "x"],
              [
                AS_IF([test "x$NETCDF_INCLUDE" != "x"], [AS_IF([test -d $NETCDF_INCLUDE], [NETCDF_INC=$NETCDF_INCLUDE])],
                      [test -d $NETCDF_LS_CHECK],       [NETCDF_INC=$NETCDF_LS_CHECK],
                      [test "x$NETCDF_DIR" != "x"],     [NETCDF_INC=$NETCDF_DIR/include],
                      [NETCDF_INC="/usr/include"])
              ])

        AS_IF([test "x$NETCDF_LIB" = "x"],
              [
                AS_IF([test "x$NETCDF_DIR" != x], [NETCDF_LIB=$NETCDF_DIR/lib],
                      [NETCDF_LIB="/usr/lib"])
              ])

        dnl Properly let the substitution variables

        dnl Check for existence of a header file in the specified location
        netcdfincFound=no;

        AC_CHECK_HEADERS($NETCDF_INC/netcdf.h, netcdfincFound=yes)

        AS_IF([test "$netcdfincFound" = "no"],
              [
                AC_MSG_RESULT(NETCDF header files not found!)
                enablenetcdf=no;
              ])

        dnl Discover the major, minor, and build versions of NETCDF by looking in
        dnl netcdf_meta.h
        AS_IF([test "x$enablenetcdf" = "xyes"],
              [
                dnl If we have the netcdf_meta.h (NETCDF 4.x), find the version there
                AS_IF([test -r $NETCDF_INC/netcdf_meta.h],
                      [
                      netcdfmajor=`grep "define NC_VERSION_MAJOR" /usr/include/netcdf_meta.h | sed -e "s/.*define NC_VERSION_MAJOR[ ]*\(.\) .*/\1/g"`
                      netcdfminor=`grep "define NC_VERSION_MINOR" /usr/include/netcdf_meta.h | sed -e "s/.*define NC_VERSION_MINOR[ ]*\(.\) .*/\1/g"`
                      netcdfpatch=`grep "define NC_VERSION_PATCH" /usr/include/netcdf_meta.h | sed -e "s/.*define NC_VERSION_PATCH[ ]*\(.\) .*/\1/g"`
                      ])
                netcdfversion=$netcdfmajor.$netcdfminor.$netcdfpatch
                netcdfmajorminor=$netcdfmajor.$netcdfminor
              ])

        dnl fix for --disable-optional
        AS_IF([test "x$enablenetcdf" = "xno"], [netcdfmajor=no])

        dnl netCDF3 is no longer distributed with libMesh
        AS_IF([test "x$netcdfmajor" = "x3"],
              [
                AC_MSG_RESULT([<<< Using netCDF 3.x is no longer supported, using version 4.x instead. >>>])
                netcdfmajor=4
              ])

        AS_CASE("${netcdfmajor}",
                [3], [
                       dnl We shouldn't get here, see if test above.
                       AC_MSG_ERROR([>>> Error: netCDF3 is no longer distributed with libMesh <<<])
                     ],
                [4], [NETCDF_INCLUDE=$NETCDF_INC

                      dnl Save original value of LIBS, then append $VTK_LIB
                      old_LIBS="$LIBS"
                      old_CPPFLAGS="$CPPFLAGS"

                      NETCDF_INCLUDE="-I$NETCDF_INC"

                      # If this compiler supports -rpath commands, create a
                      # variable for them now that can be used in $LIBS below.  We
                      # ran across an issue where GCC's linker actually needed
                      # -rpath flags in order to *link* a test program.  From the
                      # man page for GNU ld:
                      #   The -rpath option is also used when locating shared objects
                      #   which are needed by shared objects explicitly included in
                      #   the link; see the description of the -rpath-link option.
                      AS_IF([test "x$RPATHFLAG" != "x" && test -d $NETCDF_LIB],
                            [NETCDF_RPATH_FLAGS="${RPATHFLAG}${NETCDF_LIB}"
                             NETCDF_LIBRARY="${RPATHFLAG}${NETCDF_LIB} $NETCDF_LIBRARY"
                            ],
                            [NETCDF_LIBRARY="-L${NETCDF_LIB} -lnetcdf $NETCDF_LIBRARY"])
                      dnl LIBS="$LIBS $NETCDF_LIBRARY"
                      dnl CPPFLAGS="$CPPFLAGS $NETCDF_INCLUDE"
                      LIBS="$LIBS -lnetcdf -L/usr/lib"
                      CPPFLAGS="$CPPFLAGS -I/usr/include"

                      dnl test compiling and linking a test program.
                      dnl AC_LINK_IFELSE (input, [action-if-true], [action-if-false])z
                      AC_LINK_IFELSE([AC_LANG_PROGRAM([_AX_CXX_COMPILE_NETCDF_preamble],[_AX_CXX_COMPILE_NETCDF_body])],
                                     [AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with Netcdf support])],
                                     [AC_MSG_ERROR(*** Linking a test program against the NETCDF libraries failed)])

                      dnl building netcdf-4 requires that we support nested subpackages
                      AS_IF([test "$NETCDF_INCLUDE" = "xno"], [AC_MSG_ERROR([NetCDF v4 requires nested subpackages, try --enable-nested])])

                      dnl pass --disable-netcdf-4 to the subpackage so that we do not require HDF-5
                      AS_IF([test "x$enablehdf5" = "xno"], [libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-netcdf-4"])

                      dnl netcdf will install its own pkgconfig script, use this to get proper static linking
                      libmesh_pkgconfig_requires="netcdf >= 4.2 $libmesh_pkgconfig_requires"
                      AC_MSG_RESULT(<<< Configuring library with NetCDF version $netcdfversion support >>>)],
                [
                  NETCDF_INCLUDE=""
                  enablenetcdf=no
                ])

        dnl allow opt-out for nested subpackages
        AS_IF([test "x$enablenested" = "xyes"],
              [
               dnl ensure that the configuration is consistent
               AS_IF([test "x$enablecurl" = "xyes"],
                     [libmesh_subpackage_arguments="$libmesh_subpackage_arguments --enable-dap" ],
                     [libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-dap --disable-curl"])
               dnl pass --disable-testsets to the netcdf subpackage to disable the most rigorous tests
               libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-testsets"
               AS_IF([test "$NETCDF_INC" = contrib/netcdf/v4/include],
                     [AC_CONFIG_SUBDIRS(contrib/netcdf/v4)])
              ])

  ])

  AC_SUBST(NETCDF_INCLUDE)
  AC_SUBST(NETCDF_LIBRARY)

  dnl If NETCDF is not enabled, but it *was* required, error out now
  dnl instead of compiling libmesh in an invalid configuration.
  AS_IF([test "$enablenetcdf" = "no" && test "$netcdfrequired" = "yes"],
        dnl We return error code 4 here, since 0 means success and 1 is
        dnl indistinguishable from other errors.  Ideally, all of the
        dnl AC_MSG_ERROR calls in our m4 files would return a different
        dnl error code, but currently this is not enforced.
        [AC_MSG_ERROR([*** NETCDF was not found, but --enable-netcdf-required was specified.], 4)])
  ])
