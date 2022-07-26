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

AC_DEFUN([CONFIGURE_NETCDF],
[
  AC_ARG_ENABLE([netcdf],
                [AS_HELP_STRING([--disable-netcdf],
                                [build without netCDF binary I/O])],
                [AS_CASE(["${enableval}"],
                         [yes], [enablenetcdf=yes],
                         [no],  [enablenetcdf=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-netcdf)])],
                [enablenetcdf=$enableoptional])

  AS_IF([test "x$enablenetcdf" = "xyes"],
        [
          AC_ARG_WITH([netcdf-include],
                [AS_HELP_STRING([--with-netcdf-include=PATH],
                                [Specify the path for NetCDF header])],
                [withnetcdfinc=$withval],
                [withnetcdfinc=no])

          AC_ARG_WITH([netcdf-lib],
                      [AS_HELP_STRING([--with-netcdf-lib=PATH],
                                      [Specify the path for NetCDF libs])],
                      [withnetcdflib=$withval],
                      [withnetcdflib=no])

          dnl first check for external netcdf library
          install_internal_netcdf=no

          dnl Start with reasonable external defaults.
          NETCDF_INC="/usr/include/netcdf"
          NETCDF_LIB="/usr/lib"
          AS_IF([test "x$withnetcdfinc" != "xno"], [NETCDF_INC="$withnetcdfinc"])
          AS_IF([test "x$withnetcdflib" != "xno"], [NETCDF_LIB="$withnetcdflib"])

          # Initialize eventual Makefile/config.h substitution variables
          NETCDF_INCLUDE=""
          NETCDF_LIBRARY=""

          external_netcdf_inc_found=no;
          ac_netcdf_save_CPPFLAGS="$CPPFLAGS"
          CPPFLAGS="-I$NETCDF_INC $CPPFLAGS"
          AC_CHECK_HEADERS([$NETCDF_INC/netcdf.h], [external_netcdf_inc_found=yes])

          AS_IF([test "x$external_netcdf_inc_found" = "xyes"], [NETCDF_INCLUDE="-I$NETCDF_INC"
                                                                NETCDF_LIBRARY="-L$NETCDF_LIB -lnetcdf"],

                [test -d $top_srcdir/contrib/netcdf/v4],  [AC_MSG_RESULT([<<< external NetCDF files not found, using NetCDF from ./contrib >>>])
                                                           NETCDF_INC=$top_srcdir/contrib/netcdf/v4/include
                                                           NETCDF_LIB=/usr/lib
                                                           NETCDF_INCLUDE="-I\$(top_srcdir)/contrib/netcdf/v4/include"
                                                           NETCDF_LIBRARY="-L$NETCDF_LIB"
                                                           install_internal_netcdf=yes],
                [enablenetcdf=no])

          CPPFLAGS="${ac_netcdf_save_CPPFLAGS}"
          AS_IF([test "x$install_internal_netcdf" = "xyes"],
                [
                   CPPFLAGS="$NETCDF_INCLUDE $CPPFLAGS"
                   dnl Do not use cached results for header checks.
                   AS_UNSET([ac_cv_header_netcdf_h])
                   AC_CHECK_HEADERS([netcdf.h], [internal_netcdf_found=yes])
                   AS_IF([test "x$internal_netcdf_found" = "xno"],
                         [AC_MSG_RESULT([<<< internal NetCDF files not found >>>])
                          enablenetcdf=no])
                   CPPFLAGS="${ac_netcdf_save_CPPFLAGS}"
                ])

          dnl Check to make sure the NetCDF version is supported.
          AS_IF([test "x$enablenetcdf" = "xyes"],
                [
                  AS_IF([test -r "$NETCDF_INC/netcdf_meta.h"],
                        [
                          netcdfmajor=`grep "define NC_VERSION_MAJOR" $NETCDF_INC/netcdf_meta.h | sed -e "s/.*define NC_VERSION_MAJOR[ ]*\(.\) .*/\1/g"`
                          netcdfminor=`grep "define NC_VERSION_MINOR" $NETCDF_INC/netcdf_meta.h | sed -e "s/.*define NC_VERSION_MINOR[ ]*\(.\) .*/\1/g"`
                          netcdfpatch=`grep "define NC_VERSION_PATCH" $NETCDF_INC/netcdf_meta.h | sed -e "s/.*define NC_VERSION_PATCH[ ]*\(.\) .*/\1/g"`
                        ])

                  netcdfversion=$netcdfmajor

                  AS_IF([test "x$netcdfversion" = "x3"],
                        [AC_MSG_RESULT([<<< Using netCDF 3.x is no longer supported >>>])
                         enablenetcdf=no
                        ])
                ])

          dnl If we are using an external library, lets try to link it.
          AS_IF([test "x$enablenetcdf" = "xyes" && test "x$install_internal_netcdf" = "xno"],
                [
                  ac_netcdf_save_CPPFLAGS="$CPPFLAGS"
                  ac_netcdf_save_LIBS="$LIBS"
                  CPPFLAGS="-I$NETCDF_INC $CPPFLAGS"

                  AS_IF([test "x$RPATHFLAG" != "x" && test -d $NETCDF_LIB],
                        [
                          NETCDF_RPATH_FLAGS="${RPATHFLAG}${NETCDF_LIB}"
                          NETCDF_LIBRARY="${RPATHFLAG}${NETCDF_LIB} ${NETCDF_LIBRARY}"
                        ],
                        [
                          NETCDF_LIBRARY="-L${NETCDF_LIB} -lnetcdf"
                        ])
                  LIBS="$LIBS -lnetcdf -L/usr/lib"
                  AC_LINK_IFELSE([AC_LANG_PROGRAM([_AX_CXX_COMPILE_NETCDF_preamble],
                                                  [_AX_CXX_COMPILE_NETCDF_body])],
                                 [AC_MSG_RESULT(yes)],
                                 [AC_MSG_RESULT(no); enablenetcdf=no])
                  CPPFLAGS="${ac_netcdf_save_CPPFLAGS}"
                  LIBS="${ac_netcdf_save_LIBS}"
                ])

          AS_IF([test "x$enablenetcdf" = "xyes"],
                [
                  HAVE_NETCDF=1
                  AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with NetCDF support])

                  AS_IF([test "x$enablenested" = "xno"],
                        [AC_MSG_ERROR([NetCDF v4 requires nested subpackages, try --enable-nested])])

                  dnl pass --disable-netcdf-4 to the subpackage so that we do not require HDF-5
                  AS_IF([test "x$enablehdf5" = "xno"],
                        [libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-netcdf-4"])

                  dnl netcdf will install its own pkgconfig script, use this to get proper static linking
                  libmesh_pkgconfig_requires="netcdf >= 4.2 $libmesh_pkgconfig_requires"
                  AC_MSG_RESULT([<<< Configuring library with NetCDF version $netcdfversion support >>>])
                ],
                [AC_MSG_ERROR([<<< NetCDF could not be configured >>>])])
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
          AS_IF([test "$NETCDF_INC" = $top_srcdir/contrib/netcdf/v4/include],
                [AC_CONFIG_SUBDIRS([contrib/netcdf/v4])])
        ])

  AC_SUBST(NETCDF_INCLUDE)
  AC_SUBST(NETCDF_LIBRARY)
])
