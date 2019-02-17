dnl -------------------------------------------------------------
dnl netCDF
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NETCDF],
[
  AC_ARG_ENABLE(netcdf,
                AS_HELP_STRING([--disable-netcdf],
                               [build without netCDF binary I/O]),
                [AS_CASE("${enableval}",
                  [yes|new|v4], [enablenetcdf=yes
                                 netcdfversion=4],
                  [old|v3],     [enablenetcdf=yes
                                 netcdfversion=3],
                  [no],         [enablenetcdf=no
                                 netcdfversion=no],
                  [AC_MSG_ERROR(bad value ${enableval} for --enable-netcdf)])],
                [enablenetcdf=$enableoptional; netcdfversion=4])

  dnl fix for --disable-optional
  AS_IF([test "x$enablenetcdf" = "xno"], [netcdfversion=no])

  dnl netCDF3 is no longer distributed with libMesh
  AS_IF([test "x$netcdfversion" = "x3"],
        [
          AC_MSG_RESULT([<<< Using netCDF 3.x is no longer supported, using version 4.x instead. >>>])
          netcdfversion=4
        ])

  AS_CASE("${netcdfversion}",
          [3], [
                 dnl We shouldn't get here, see if test above.
                 AC_MSG_ERROR([>>> Error: netCDF3 is no longer distributed with libMesh <<<])
               ],
          [4], [NETCDF_INCLUDE="-I\$(top_srcdir)/contrib/netcdf/v4/include"
                AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with Netcdf support])
                dnl building netcdf-4 requires that we support nested subpackages
                AS_IF([test "x$enablenested" = "xno"], [AC_MSG_ERROR([NetCDF v4 requires nested subpackages, try --enable-nested])])

                dnl pass --disable-netcdf-4 to the subpackage so that we do not require HDF-5
                AS_IF([test "x$enablehdf5" = "xno"], [libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-netcdf-4"])

                dnl netcdf will install its own pkgconfig script, use this to get proper static linking
                libmesh_pkgconfig_requires="netcdf >= 4.2 $libmesh_pkgconfig_requires"
                AC_MSG_RESULT(<<< Configuring library with NetCDF version 4 support >>>)],
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
          AC_CONFIG_SUBDIRS([contrib/netcdf/v4])
        ])

  AC_SUBST(NETCDF_INCLUDE)
])
