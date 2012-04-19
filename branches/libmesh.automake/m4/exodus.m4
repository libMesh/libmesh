dnl -------------------------------------------------------------
dnl exodus
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_EXODUS], 
[
  AC_ARG_ENABLE(exodus,
                AC_HELP_STRING([--enable-exodus],
                               [build with ExodusII API support]),
		[case "${enableval}" in
		  yes)  enableexodus=yes ;;
		   no)  enableexodus=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-exodus) ;;
		 esac],
		 [enableexodus=$enablenetcdf]) # if unspecified, depend on netcdf


		
  dnl The EXODUS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enableexodus = yes); then
     EXODUS_INCLUDE="-I\$(top_srcdir)/contrib/exodusii/Lib/include"
     AC_DEFINE(HAVE_EXODUS_API, 1, [Flag indicating whether the library will be compiled with Exodus support])
     AC_MSG_RESULT(<<< Configuring library with Exodus support >>>)
     libmesh_contrib_INCLUDES="$EXODUS_INCLUDE $libmesh_contrib_INCLUDES"
     have_exodus=yes
  else
     EXODUS_INCLUDE=""
     enableexodus=no
     have_exodus=no
  fi

  AC_SUBST(EXODUS_INCLUDE)
  AC_SUBST(enableexodus)
		 
  AM_CONDITIONAL(LIBMESH_ENABLE_EXODUS, test x$enableexodus = xyes)
])
