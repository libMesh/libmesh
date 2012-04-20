dnl -------------------------------------------------------------
dnl Metis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_METIS], 
[
  AC_ARG_ENABLE(metis,
                AC_HELP_STRING([--enable-metis],
                               [build with Metis graph partitioning suppport]),
		[case "${enableval}" in
		  yes)  enablemetis=yes ;;
		   no)  enablemetis=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-metis) ;;
		 esac],
		 [enablemetis=$enableoptional])



  dnl The METIS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablemetis = yes); then
     METIS_INCLUDE="-I\$(top_srcdir)/contrib/metis/Lib"
     METIS_LIB="\$(EXTERNAL_LIBDIR)/libmetis\$(libext)"
     AC_DEFINE(HAVE_METIS, 1, [Flag indicating whether the library will be compiled with Metis support])
     AC_MSG_RESULT(<<< Configuring library with Metis support >>>)
     libmesh_contrib_INCLUDES="$METIS_INCLUDE $libmesh_contrib_INCLUDES"
  else
     METIS_INCLUDE=""
     METIS_LIB=""
     enablemetis=no
  fi

  AC_SUBST(METIS_INCLUDE)
  AC_SUBST(METIS_LIB)	
  AC_SUBST(enablemetis)

  AM_CONDITIONAL(LIBMESH_ENABLE_METIS, test x$enablemetis = xyes)		 
])
