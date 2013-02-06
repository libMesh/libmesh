dnl -------------------------------------------------------------
dnl exodus
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_EXODUS], 
[
  AC_ARG_ENABLE(exodus,
                AC_HELP_STRING([--enable-exodus],
                               [build with ExodusII API support]),
		[case "${enableval}" in
		  yes|v509)  enableexodus=yes; exodusversion="v5.09" ;;
		  new|v522)  enableexodus=yes; exodusversion="v5.22" ;;
		   no)       enableexodus=no;  exodusversion=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-exodus) ;;
		 esac],
		 [enableexodus=$enablenetcdf ; exodusversion="v5.09"]) # if unspecified, depend on netcdf


  # fix for --disable-optional
  if (test "x$enableexodus" = "xno"); then
    exodusversion=no
  fi   
		
  case "${exodusversion}" in
      
      "v5.09")
          # The EXODUS API is distributed with libmesh, so we don't have to guess
          # where it might be installed...
	  EXODUS_INCLUDE="-I\$(top_srcdir)/contrib/exodusii/$exodusversion/include"
	  AC_DEFINE(HAVE_EXODUS_API, 1, [Flag indicating whether the library will be compiled with Exodus support])
	  AC_MSG_RESULT(<<< Configuring library with Exodus version $exodusversion support >>>)
	  ;;

      "v5.22")
	  EXODUS_INCLUDE="-I\$(top_srcdir)/contrib/exodusii/$exodusversion/exodus/cbind/include"
	  AC_DEFINE(HAVE_EXODUS_API, 1, [Flag indicating whether the library will be compiled with Exodus support])
	  AC_MSG_RESULT(<<< Configuring library with Exodus version $exodusversion support >>>)
	  ;;

      *)
	  EXODUS_INCLUDE=""
	  enableexodus=no
	  ;;
  esac

  AC_CONFIG_FILES([contrib/exodusii/v5.09/Makefile])
  AC_CONFIG_FILES([contrib/exodusii/v5.22/exodus/cbind/Makefile])
  AC_SUBST(EXODUS_INCLUDE)
])
