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
		 [enablenetcdf=$enableoptional; netcdfversion=3])				
		 		
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
	  ;;

      4)
	  NETCDF_INCLUDE="-I\$(top_srcdir)/contrib/netcdf/v4/include"
	  AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with Netcdf support])
	  # building netcdf-4 requires that we support nested subpackages
	  if (test "x$enablenested" = "xno"); then
	      AC_MSG_ERROR([NetCDF v4 requres nested subpackages, try --enable-nested])
	  fi								
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
      AC_CONFIG_SUBDIRS([contrib/netcdf/v4])

      if (test "x$enablehdf5" = "xno"); then	    
  	  #  pass --disable-netcdf-4 to the subpackage so that we do not require HDF-5
	  libmesh_subpackage_arguments="$libmesh_subpackage_arguments --disable-netcdf-4"
      fi
  fi

  AC_SUBST(NETCDF_INCLUDE)
])
