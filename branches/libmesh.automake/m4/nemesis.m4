dnl -------------------------------------------------------------
dnl nemesis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NEMESIS], 
[
  AC_ARG_ENABLE(nemesis,
                AC_HELP_STRING([--enable-nemesis],
                               [build with NemesisII API support]),
		[case "${enableval}" in
		  yes)  enablenemesis=yes ;;
		   no)  enablenemesis=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-nemesis) ;;
		 esac],
		 [enablenemesis=yes])


		 		
  dnl The NEMESIS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablenemesis = yes); then
     NEMESIS_INCLUDE="-I\$(top_srcdir)/contrib/nemesis/Lib"
     AC_DEFINE(HAVE_NEMESIS_API, 1, [Flag indicating whether the library will be compiled with Nemesis support])
     AC_MSG_RESULT(<<< Configuring library with Nemesis support >>>)
     have_nemesis=yes
  else
     NEMESIS_INCLUDE=""
     enablenemesis=no
     have_nemesis=no
  fi

  AC_SUBST(NEMESIS_INCLUDE)
  AC_SUBST(enablenemesis)

  AM_CONDITIONAL(LIBMESH_ENABLE_NEMESIS, test x$enablenemesis = xyes)
  AC_CONFIG_FILES([contrib/nemesis/Lib/Makefile])
])
