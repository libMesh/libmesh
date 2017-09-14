dnl -------------------------------------------------------------
dnl netCDF
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NETCDF],
[
  AC_ARG_ENABLE(netcdf,
                AS_HELP_STRING([--disable-netcdf],
                               [build without netCDF binary I/O]),
                [case "${enableval}" in
                  yes|new|v4) enablenetcdf=yes; netcdfversion=4  ;;
                  old|v3) enablenetcdf=yes; netcdfversion=3  ;;
                  no) enablenetcdf=no;  netcdfversion=no ;;
                  *) AC_MSG_ERROR(bad value ${enableval} for --enable-netcdf) ;;
                esac],
                [enablenetcdf=$enableoptional; netcdfversion=4])

  # fix for --disable-optional
  if (test "x$enablenetcdf" = "xno"); then
    netcdfversion=no
  fi

  case "${netcdfversion}" in
      3)
        # The NETCDF API is distributed with libmesh, so we don't have to guess
        # where it might be installed...
        NETCDF_INCLUDE="-I\$(top_srcdir)/contrib/netcdf/v3"
        AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with Netcdf support])
        AC_MSG_RESULT(<<< Configuring library with NetCDF version 3 support >>>)

        # pass --disable-netcdf-4 to the subpackage so that we do not require HDF-5
        # note this is madness - we will run configure in the subdirectory v4, but not use it,
        # so this is nothing more than a hedge against that failing.  we need it to work for
        # 'make dist' to work
        libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-netcdf-4"
        ;;

      4)
        NETCDF_INCLUDE="-I\$(top_srcdir)/contrib/netcdf/v4/include"
        AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with Netcdf support])
        # building netcdf-4 requires that we support nested subpackages
        if (test "x$enablenested" = "xno"); then
          AC_MSG_ERROR([NetCDF v4 requires nested subpackages, try --enable-nested])
        fi

        if (test "x$enablehdf5" = "xno"); then
          # pass --disable-netcdf-4 to the subpackage so that we do not require HDF-5
          libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-netcdf-4"
        fi

        # netcdf will install its own pkgconfig script, use this to get proper static linking
        libmesh_pkgconfig_requires="netcdf >= 4.2 $libmesh_pkgconfig_requires"
        AC_MSG_RESULT(<<< Configuring library with NetCDF version 4 support >>>)
        ;;

      *)
        NETCDF_INCLUDE=""
        enablenetcdf=no
        ;;
  esac

  AC_CONFIG_FILES([contrib/netcdf/v3/Makefile])

  # allow opt-out for nested subpackages
  if (test "x$enablenested" = "xyes"); then
      #  pass --disable-testsets to the netcdf subpackage to disable the most rigorous tests
      libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-testsets"
      AC_CONFIG_SUBDIRS([contrib/netcdf/v4])
  fi

  AC_SUBST(NETCDF_INCLUDE)
])
