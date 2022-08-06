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

  netcdf_v4_arg=''
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
                AS_IF([test "x$enablehdf5" = "xno"], [netcdf_v4_arg=--disable-netcdf-4])

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
          dnl We need our netcdf subpackage to be able to use whatever HDF5 we've
          dnl detected.  That means using whatever $HDF5_CPPFLAGS and $HDF5_LIBS
          dnl we've determined work - but if we're using PETSc it also means
          dnl passing along PETSc's directories in case we want to link to a
          dnl --download-hdf5 build there.
          SUB_CPPFLAGS="$CPPFLAGS"
          SUB_LIBS="$LIBS"
          AS_IF([test $enablehdf5 = yes],
                [
                 SUB_CPPFLAGS="$HDF5_CPPFLAGS $SUB_CPPFLAGS"
                 SUB_LIBS="$HDF5_LIBS $SUB_LIBS"
                 AS_IF([test "x$enablepetsc" != "xno"],
                       [
                        SUB_CPPFLAGS="$PETSCINCLUDEDIRS $SUB_CPPFLAGS"
                        SUB_LIBS="$PETSCLINKLIBS $SUB_LIBS"
                       ])])

          dnl ensure that the configuration is consistent
          dnl pass --disable-testsets to the netcdf subpackage to disable the most rigorous tests
          AS_IF([test "x$enablecurl" = "xyes"],
                [netcdf_dap_arg=--enable-dap
                 netcdf_curl_arg=''],
                [netcdf_dap_arg=--disable-dap
                 netcdf_curl_arg=--disable-curl])
          AX_SUBDIRS_CONFIGURE([contrib/netcdf/v4],[[CXX=$CXX],[CC=$CC],[F77=$F77],[FC=$FC],[CPPFLAGS=$SUB_CPPFLAGS],[LIBS=$SUB_LIBS],[$netcdf_dap_arg],[$netcdf_curl_arg],[--disable-testsets],[$netcdf_v4_arg]])
        ])

  AC_SUBST(NETCDF_INCLUDE)
])
