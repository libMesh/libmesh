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

  AS_CASE("${netcdfversion}",
          [3], [dnl The NETCDF API is distributed with libmesh, so we don't have to guess
                dnl where it might be installed...
                NETCDF_INCLUDE="-I\$(top_srcdir)/contrib/netcdf/v3"
                AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with Netcdf support])
                AC_MSG_RESULT(<<< Configuring library with NetCDF version 3 support >>>)
                dnl pass --disable-netcdf-4 to the subpackage so that we do not require HDF-5
                dnl note this is madness - we will run configure in the subdirectory v4, but not use it,
                dnl so this is nothing more than a hedge against that failing.  we need it to work for
                dnl 'make dist' to work
                libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-netcdf-4"],
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

  AC_CONFIG_FILES([contrib/netcdf/v3/Makefile])

  dnl allow opt-out for nested subpackages
  AS_IF([test "x$enablenested" = "xyes"],
        [
          dnl pass --disable-testsets to the netcdf subpackage to disable the most rigorous tests
          libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-testsets"
          AC_CONFIG_SUBDIRS([contrib/netcdf/v4])
        ])

  AC_SUBST(NETCDF_INCLUDE)
])
