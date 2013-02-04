dnl -------------------------------------------------------------
dnl netCDF
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NETCDF], 
[
  AC_ARG_ENABLE(netcdf,
                AC_HELP_STRING([--enable-netcdf],
                               [build with netCDF binary I/O]),
		[case "${enableval}" in
		  yes|v3) enablenetcdf=yes; netcdfversion=3  ;;
                  new|v4) enablenetcdf=yes; netcdfversion=4  ;;
		   no)    enablenetcdf=no;  netcdfversion=no ;;
 		    *)    AC_MSG_ERROR(bad value ${enableval} for --enable-netcdf) ;;
		 esac],
		 [enablenetcdf=$enableoptional])


				
  if (test "x$netcdfversion" = "x3"); then
     # The NETCDF API is distributed with libmesh, so we don't have to guess
     # where it might be installed...
     NETCDF_INCLUDE="-I\$(top_srcdir)/contrib/netcdf/v3"
     AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with Netcdf support])
     AC_MSG_RESULT(<<< Configuring library with NetCDF version 3 support >>>)
     AC_CONFIG_FILES([contrib/netcdf/v3/Makefile])

  elif (test "x$netcdfversion" = "x4"); then
     NETCDF_INCLUDE="-I\$(top_srcdir)/contrib/netcdf/v4/include"
     AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with Netcdf support])
     AC_MSG_RESULT(<<< Configuring library with NetCDF version 4 support >>>)
     AC_CONFIG_SUBDIRS([contrib/netcdf/v4])

  else
     NETCDF_INCLUDE=""
     enablenetcdf=no
  fi
 
  AC_SUBST(NETCDF_INCLUDE)
])
